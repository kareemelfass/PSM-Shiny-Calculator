app_dir <- normalizePath(file.path("..", ".."))
for (f in list.files(file.path(app_dir, "R"), pattern = "^fct_|^utils", full.names = TRUE)) {
  source(f, local = FALSE)
}
library(survival)
example_ipd <- read_ipd(file.path(app_dir, "example", "example_ipd_months.csv"), "months")
example_lt <- read_lifetable(file.path(app_dir, "example", "illustrative_lifetable.csv"))
example_arms <- list(comp = "Standard care", int = "New drug")

# Exponential "pieces" with closed-form answers (no fitting needed)
exp_pieces <- function(spec, t, lam_os = 0.3, lam_pfs = 0.8, lam_os_int = 0.2, lam_pfs_int = 0.5) {
  list(
    t = t,
    OS = list(comp = exp(-lam_os * t), int = exp(-lam_os_int * t)),
    PFS = list(comp = exp(-lam_pfs * t), int = exp(-lam_pfs_int * t)),
    ToT = list(comp = NULL, int = NULL),
    gp = NULL
  )
}
