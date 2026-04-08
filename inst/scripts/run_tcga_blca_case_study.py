#!/usr/bin/env python
"""Run the TCGA-BLCA case study."""

from __future__ import annotations

import sys
import tempfile
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from run_matrix_case_study import run_case_study


def main() -> None:
    output_dir = ROOT / "inst" / "extdata" / "case_studies" / "tcga_blca"
    expr_url = "https://raw.githubusercontent.com/ilyada/NestedWGCNA/main/data/TCGA_BLCA.tsv"

    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td) / "tcga_blca.tsv"
        urllib.request.urlretrieve(expr_url, tmp)
        expr = pd.read_csv(tmp, sep="\t", index_col=0)

    run_case_study(
        expr=expr,
        output_dir=output_dir,
        dataset_id="tcga_blca",
        source=expr_url,
        cgm_min_cluster_size=50,
        fgm_min_cluster_size=10,
        top_variable_genes=1200,
        metadata=None,
    )


if __name__ == "__main__":
    main()
