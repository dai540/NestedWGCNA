if (!requireNamespace("testthat", quietly = TRUE)) {
  quit(save = "no", status = 0)
}

library(testthat)
library(NestedWGCNA)

test_check("NestedWGCNA")
