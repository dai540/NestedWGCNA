#!/usr/bin/env python
"""Run NestedWGCNA reference analysis on a sample x gene matrix.

Inputs:
- TSV expression matrix (rows: samples, cols: genes)
- Optional TSV metadata table indexed by sample ids

Outputs:
- module assignments, enrichment tables, module score tables
- phenotype association table (if metadata is provided)
- plots and summary JSON
"""

from __future__ import annotations

import argparse
import json
import sys
import urllib.request
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import f_oneway, hypergeom, spearmanr


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
PYREF = ROOT / "inst" / "python_ref"
if PYREF.exists() and str(PYREF) not in sys.path:
    sys.path.insert(0, str(PYREF))

import GenFocus
import umap.umap_ as uu
from GenFocus import genfocus
from utils import get_clust


XCELL_URL = "https://raw.githubusercontent.com/ilyada/NestedWGCNA/main/data/xCell_signatures.csv"
BOSTONGENE_URL = (
    "https://raw.githubusercontent.com/ilyada/NestedWGCNA/main/data/bostongene_signatures.txt"
)
RANDOM_STATE = 42


def patch_runtime_compatibility() -> None:
    """Patch runtime behavior for reproducibility and modern sklearn compatibility."""
    original = uu.check_array

    def _check_array_compat(*args, force_all_finite=None, **kwargs):
        if force_all_finite is not None and "ensure_all_finite" not in kwargs:
            kwargs["ensure_all_finite"] = force_all_finite
        return original(*args, **kwargs)

    uu.check_array = _check_array_compat

    # Make set operations deterministic in GenFocus.
    GenFocus.exclusion = lambda lst1, lst2: sorted(set(lst1) - set(lst2))
    GenFocus.intersection = lambda lst1, lst2: sorted(set(lst1) & set(lst2))


def download(url: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(url, path)


def parse_xcell_csv(path: Path) -> dict[str, set[str]]:
    signatures: dict[str, set[str]] = {}
    with path.open("r", encoding="utf-8") as f:
        for i, line in enumerate(f):
            parts = [p.strip() for p in line.rstrip("\n").split(",")]
            if i == 0 or not parts or parts[0] == "":
                continue
            sig = parts[0]
            genes = [g for g in parts[2:] if g != ""]
            if genes:
                signatures[sig] = set(genes)
    return signatures


def parse_bostongene_tsv(path: Path) -> dict[str, set[str]]:
    signatures: dict[str, set[str]] = {}
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            parts = [p.strip() for p in line.rstrip("\n").split("\t")]
            if not parts or parts[0] == "":
                continue
            sig = parts[0]
            genes = [g for g in parts[2:] if g != ""]
            if genes:
                signatures[sig] = set(genes)
    return signatures


def bh_fdr(pvals: Iterable[float]) -> np.ndarray:
    pvals = np.asarray(list(pvals), dtype=float)
    if len(pvals) == 0:
        return pvals
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order] * n / (np.arange(1, n + 1))
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(n, dtype=float)
    q[order] = ranked
    return np.minimum(q, 1.0)


def normalize_input(expr: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, float]]:
    """Normalize to stable non-negative values for downstream Python reference code."""
    expr2 = expr.copy().astype(float)
    min_val = float(np.nanmin(expr2.values))
    shifted = 0.0
    if min_val < 0:
        shifted = -min_val + 1e-6
        expr2 = expr2 + shifted

    # Heuristic log-transform for count-like scales.
    q99 = float(np.nanquantile(expr2.values, 0.99))
    logged = 0
    if q99 > 50:
        expr2 = np.log1p(expr2)
        logged = 1

    meta = {
        "input_min_before_shift": min_val,
        "shift_applied": shifted,
        "log1p_applied": float(logged),
        "q99_after_shift_before_log": q99,
    }
    return expr2, meta


def top_var_genes(expr: pd.DataFrame, n: int) -> pd.DataFrame:
    if n >= expr.shape[1]:
        return expr.copy()
    selected = expr.var(axis=0).sort_values(ascending=False).head(n).index
    return expr.loc[:, selected]


def enrichment_table(
    module_genes: list[str],
    signatures: dict[str, set[str]],
    source: str,
    universe: set[str],
) -> pd.DataFrame:
    rows = []
    module = set(module_genes)
    m_total = len(universe)
    k_total = len(module)

    for sig, sig_genes in signatures.items():
        genes = sig_genes & universe
        n_sig = len(genes)
        if n_sig < 5:
            continue
        overlap = module & genes
        x = len(overlap)
        if x == 0:
            continue
        p_value = hypergeom.sf(x - 1, m_total, n_sig, k_total)
        rows.append(
            (
                source,
                sig,
                m_total,
                k_total,
                n_sig,
                x,
                p_value,
                "|".join(sorted(overlap)),
            )
        )

    if not rows:
        return pd.DataFrame(
            columns=["source", "signature", "M", "K", "n", "k", "p", "overlap_genes", "fdr"]
        )

    out = pd.DataFrame(
        rows,
        columns=["source", "signature", "M", "K", "n", "k", "p", "overlap_genes"],
    ).sort_values("p")
    out["fdr"] = bh_fdr(out["p"].to_numpy())
    return out.reset_index(drop=True)


def module_scores(expr: pd.DataFrame, assignment: pd.Series, prefix: str) -> pd.DataFrame:
    rows = {}
    for module_id in sorted(assignment[assignment != -1].unique()):
        genes = assignment[assignment == module_id].index
        if len(genes) == 0:
            continue
        rows[f"{prefix}{int(module_id)}"] = expr.loc[:, genes].mean(axis=1)
    return pd.DataFrame(rows, index=expr.index)


def phenotype_association(scores: pd.DataFrame, metadata: pd.DataFrame, layer: str) -> pd.DataFrame:
    if scores.empty or metadata.empty:
        return pd.DataFrame(
            columns=[
                "layer",
                "module",
                "phenotype",
                "test",
                "levels_or_n",
                "effect",
                "p",
                "fdr",
            ]
        )

    common = scores.index.intersection(metadata.index)
    if len(common) < 10:
        return pd.DataFrame(
            columns=[
                "layer",
                "module",
                "phenotype",
                "test",
                "levels_or_n",
                "effect",
                "p",
                "fdr",
            ]
        )

    s = scores.loc[common]
    m = metadata.loc[common]
    out = []

    for col in m.columns:
        trait = m[col]
        # numeric: Spearman correlation
        if pd.api.types.is_numeric_dtype(trait):
            if trait.notna().sum() < 10:
                continue
            for mod in s.columns:
                valid = trait.notna() & s[mod].notna()
                if valid.sum() < 10:
                    continue
                rho, p = spearmanr(s.loc[valid, mod], trait.loc[valid])
                out.append((layer, mod, col, "spearman", int(valid.sum()), float(rho), float(p)))
            continue

        # categorical: one-way ANOVA
        trait2 = trait.astype(str)
        counts = trait2.value_counts()
        keep = counts[counts >= 3].index
        if len(keep) < 2 or len(keep) > 12:
            continue
        mask = trait2.isin(keep)
        for mod in s.columns:
            groups = [s.loc[mask & (trait2 == k), mod].values for k in keep]
            try:
                stat, p = f_oneway(*groups)
                out.append(
                    (layer, mod, col, "anova", int(len(keep)), float(stat), float(p))
                )
            except Exception:
                continue

    if not out:
        return pd.DataFrame(
            columns=[
                "layer",
                "module",
                "phenotype",
                "test",
                "levels_or_n",
                "effect",
                "p",
                "fdr",
            ]
        )
    df = pd.DataFrame(
        out,
        columns=["layer", "module", "phenotype", "test", "levels_or_n", "effect", "p"],
    )
    df["fdr"] = bh_fdr(df["p"].values)
    return df.sort_values("fdr").reset_index(drop=True)


def choose_target_cgm(cgm_assign: pd.Series, cgm_enrichment: pd.DataFrame) -> tuple[int, str]:
    immune_pattern = (
        "B-cell|B_cells|T-cell|T_cells|NK|Macroph|Dendritic|Immune|MHC|Cytotoxic|MDSC|Treg|Th1|Th2"
    )
    if len(cgm_enrichment):
        tmp = cgm_enrichment.copy()
        tmp["immune_like"] = tmp["signature"].str.contains(immune_pattern, case=False, regex=True)
        cand = tmp[tmp["immune_like"]]
        if len(cand):
            best = cand.sort_values("fdr").iloc[0]
            return int(best["module"]), f"immune_signature:{best['signature']}"
    largest = int(cgm_assign[cgm_assign != -1].value_counts().idxmax())
    return largest, "largest_cgm"


def run_fgm_with_fallback(
    normalized: pd.DataFrame,
    initial_min_cluster_size: int,
) -> tuple[np.ndarray, np.ndarray, int]:
    tried = []
    for x in [initial_min_cluster_size, 10, 8, 6, 5, 4]:
        y = int(max(4, x))
        if y not in tried:
            tried.append(y)

    best = None
    for mcs in tried:
        labels, core = get_clust(normalized, min_clust_size=mcs, random_state=RANDOM_STATE)
        labels_s = pd.Series(labels)
        n_modules = int(labels_s[labels_s != -1].nunique())
        n_assigned = int((labels_s != -1).sum())
        score = (n_modules, n_assigned)
        if best is None or score > best["score"]:
            best = {"labels": labels, "core": core, "mcs": mcs, "score": score}
        if n_modules >= 2 and n_assigned >= 40:
            return labels, core, mcs
    return best["labels"], best["core"], best["mcs"]


def save_size_plot(
    assignment: pd.Series,
    title: str,
    outpath: Path,
    color: str,
) -> None:
    sizes = assignment[assignment != -1].value_counts().sort_index()
    if len(sizes) == 0:
        return
    plt.figure(figsize=(7, 4))
    plt.bar([str(i) for i in sizes.index], sizes.values, color=color)
    plt.title(title)
    plt.xlabel("module id")
    plt.ylabel("genes")
    plt.tight_layout()
    plt.savefig(outpath, dpi=160)
    plt.close()


def run_case_study(
    expr: pd.DataFrame,
    output_dir: Path,
    dataset_id: str,
    source: str,
    cgm_min_cluster_size: int,
    fgm_min_cluster_size: int,
    top_variable_genes: int,
    metadata: pd.DataFrame | None = None,
) -> dict:
    patch_runtime_compatibility()
    output_dir.mkdir(parents=True, exist_ok=True)

    meta = metadata if metadata is not None else pd.DataFrame(index=expr.index)
    expr2, input_meta = normalize_input(expr)
    expr_model = top_var_genes(expr2, n=top_variable_genes)

    sig_dir = ROOT / "inst" / "extdata" / "case_studies" / "_shared_signatures"
    sig_dir.mkdir(parents=True, exist_ok=True)
    xcell_path = sig_dir / "xCell_signatures.csv"
    boston_path = sig_dir / "bostongene_signatures.txt"
    if not xcell_path.exists():
        download(XCELL_URL, xcell_path)
    if not boston_path.exists():
        download(BOSTONGENE_URL, boston_path)
    xcell = parse_xcell_csv(xcell_path)
    boston = parse_bostongene_tsv(boston_path)

    cgm_labels, cgm_core = get_clust(
        expr_model,
        min_clust_size=cgm_min_cluster_size,
        random_state=RANDOM_STATE,
    )
    cgm_assign = pd.Series(cgm_labels, index=expr_model.columns, name="cgm")
    universe = set(expr_model.columns)

    cgm_enrichment_rows = []
    for mod, genes in cgm_assign[cgm_assign != -1].groupby(cgm_assign[cgm_assign != -1]):
        gene_list = genes.index.tolist()
        enr = pd.concat(
            [
                enrichment_table(gene_list, xcell, "xCell", universe),
                enrichment_table(gene_list, boston, "BostonGene", universe),
            ],
            ignore_index=True,
        )
        if len(enr) == 0:
            continue
        enr["module"] = int(mod)
        cgm_enrichment_rows.append(enr)
    cgm_enrichment = (
        pd.concat(cgm_enrichment_rows, ignore_index=True) if cgm_enrichment_rows else pd.DataFrame()
    )

    target_cgm, target_reason = choose_target_cgm(cgm_assign, cgm_enrichment)
    target_genes = cgm_assign[cgm_assign == target_cgm].index.tolist()
    target_expr = expr_model.loc[:, target_genes]

    genfocus_result = None
    used_threshold = None
    for threshold in [0.9, 0.85, 0.8, 0.75, 0.7]:
        try:
            result = genfocus(
                target_expr.copy(),
                focus_gene="eigengene",
                method="spearman",
                tpm=True,
                corr_thr=threshold,
                CVR_thr=0.6,
            )
        except ZeroDivisionError:
            result = None
        if result is not None:
            genfocus_result = result
            used_threshold = threshold
            break
    if genfocus_result is None:
        raise RuntimeError("GenFocus failed for tested thresholds.")

    ings, normalized = genfocus_result
    normalized = normalized.loc[target_expr.index, target_expr.columns]

    fgm_labels, fgm_core, used_fgm_mcs = run_fgm_with_fallback(normalized, fgm_min_cluster_size)
    fgm_assign = pd.Series(fgm_labels, index=normalized.columns, name="fgm")

    fgm_enrichment_rows = []
    for mod, genes in fgm_assign[fgm_assign != -1].groupby(fgm_assign[fgm_assign != -1]):
        gene_list = genes.index.tolist()
        enr = pd.concat(
            [
                enrichment_table(gene_list, xcell, "xCell", universe),
                enrichment_table(gene_list, boston, "BostonGene", universe),
            ],
            ignore_index=True,
        )
        if len(enr) == 0:
            continue
        enr["fgm"] = int(mod)
        fgm_enrichment_rows.append(enr)
    fgm_enrichment = (
        pd.concat(fgm_enrichment_rows, ignore_index=True) if fgm_enrichment_rows else pd.DataFrame()
    )

    # score tables
    cgm_scores = module_scores(expr_model, cgm_assign, "CGM_")
    fgm_scores = module_scores(normalized, fgm_assign, "FGM_")

    assoc_cgm = phenotype_association(cgm_scores, meta, "CGM")
    assoc_fgm = phenotype_association(fgm_scores, meta, "FGM")
    associations = pd.concat([assoc_cgm, assoc_fgm], ignore_index=True)

    # save tables
    pd.DataFrame(
        {"gene": expr_model.columns, "cgm": cgm_labels, "cgm_core": cgm_core}
    ).to_csv(output_dir / "cgm_assignments.tsv", sep="\t", index=False)
    pd.DataFrame(
        {"gene": normalized.columns, "fgm": fgm_labels, "fgm_core": fgm_core}
    ).to_csv(output_dir / "fgm_assignments.tsv", sep="\t", index=False)
    cgm_enrichment.to_csv(output_dir / "cgm_enrichment.tsv", sep="\t", index=False)
    fgm_enrichment.to_csv(output_dir / "fgm_enrichment.tsv", sep="\t", index=False)
    cgm_scores.to_csv(output_dir / "cgm_scores.tsv", sep="\t")
    fgm_scores.to_csv(output_dir / "fgm_scores.tsv", sep="\t")
    associations.to_csv(output_dir / "phenotype_associations.tsv", sep="\t", index=False)

    # save plots
    save_size_plot(
        cgm_assign,
        f"CGM sizes ({dataset_id})",
        output_dir / "cgm_module_sizes.png",
        color="#2a9d8f",
    )
    save_size_plot(
        fgm_assign,
        f"FGM sizes in CGM {target_cgm} ({dataset_id})",
        output_dir / "fgm_module_sizes.png",
        color="#e76f51",
    )

    summary = {
        "dataset_id": dataset_id,
        "source": source,
        "random_state": RANDOM_STATE,
        "input_meta": input_meta,
        "n_samples": int(expr_model.shape[0]),
        "n_genes_input": int(expr.shape[1]),
        "n_genes_model": int(expr_model.shape[1]),
        "cgm_min_cluster_size": int(cgm_min_cluster_size),
        "cgm_n_modules": int((cgm_assign[cgm_assign != -1]).nunique()),
        "cgm_assigned_genes": int((cgm_assign != -1).sum()),
        "cgm_noise_genes": int((cgm_assign == -1).sum()),
        "target_cgm": int(target_cgm),
        "target_cgm_reason": target_reason,
        "target_cgm_size": int(len(target_genes)),
        "genfocus_corr_thr": float(used_threshold),
        "genfocus_ings_n": int(len(ings)),
        "fgm_min_cluster_size_requested": int(fgm_min_cluster_size),
        "fgm_min_cluster_size_used": int(used_fgm_mcs),
        "fgm_n_modules": int((fgm_assign[fgm_assign != -1]).nunique()),
        "fgm_assigned_genes": int((fgm_assign != -1).sum()),
        "fgm_noise_genes": int((fgm_assign == -1).sum()),
    }

    if len(cgm_enrichment):
        top_cgm = cgm_enrichment.sort_values("fdr").groupby("module").head(3)
        summary["top_cgm_enrichment"] = [
            {
                "module": int(r.module),
                "source": r.source,
                "signature": r.signature,
                "k": int(r.k),
                "K": int(r.K),
                "fdr": float(r.fdr),
            }
            for r in top_cgm.itertuples(index=False)
        ]
    if len(fgm_enrichment):
        top_fgm = fgm_enrichment.sort_values("fdr").groupby("fgm").head(3)
        summary["top_fgm_enrichment"] = [
            {
                "fgm": int(r.fgm),
                "source": r.source,
                "signature": r.signature,
                "k": int(r.k),
                "K": int(r.K),
                "fdr": float(r.fdr),
            }
            for r in top_fgm.itertuples(index=False)
        ]
    if len(associations):
        top_assoc = associations.sort_values("fdr").head(15)
        summary["top_phenotype_associations"] = [
            {
                "layer": r.layer,
                "module": r.module,
                "phenotype": r.phenotype,
                "test": r.test,
                "effect": float(r.effect),
                "fdr": float(r.fdr),
            }
            for r in top_assoc.itertuples(index=False)
        ]

    with (output_dir / "summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="TSV matrix with samples x genes.")
    parser.add_argument("--output-dir", required=True, help="Output directory.")
    parser.add_argument("--dataset-id", required=True, help="Dataset identifier.")
    parser.add_argument("--source", required=True, help="Data source description/URL.")
    parser.add_argument("--metadata", default="", help="Optional metadata TSV (samples indexed).")
    parser.add_argument("--cgm-min-size", type=int, default=30, help="CGM min_cluster_size.")
    parser.add_argument("--fgm-min-size", type=int, default=10, help="FGM min_cluster_size.")
    parser.add_argument("--top-var-genes", type=int, default=1200, help="Number of variable genes.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    expr = pd.read_csv(args.input, sep="\t", index_col=0, low_memory=False)
    if expr.index.has_duplicates:
        expr = expr[~expr.index.duplicated(keep="first")]
    if expr.columns.has_duplicates:
        expr = expr.loc[:, ~expr.columns.duplicated(keep="first")]

    meta = None
    if args.metadata:
        md = pd.read_csv(args.metadata, sep="\t", index_col=0, low_memory=False)
        if md.index.has_duplicates:
            md = md[~md.index.duplicated(keep="first")]
        meta = md

    run_case_study(
        expr=expr,
        output_dir=Path(args.output_dir),
        dataset_id=args.dataset_id,
        source=args.source,
        cgm_min_cluster_size=args.cgm_min_size,
        fgm_min_cluster_size=args.fgm_min_size,
        top_variable_genes=args.top_var_genes,
        metadata=meta,
    )


if __name__ == "__main__":
    main()
