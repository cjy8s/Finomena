# install_packages.R — Installs required R packages into the app-local library.
#
# Run from the App/ directory:
#   Rscript R/install_packages.R
#
# This script detects its own location and installs packages into R/library/
# next to itself. The resulting library is platform-specific — run this once
# on each OS (Windows, Mac) to create a compatible library for that platform.

script_dir <- tryCatch(
  dirname(sys.frame(1)$ofile),
  error = function(e) {
    # Fallback: if run via Rscript, use commandArgs to find the script path
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      dirname(normalizePath(sub("^--file=", "", file_arg[1])))
    } else {
      stop("Cannot determine script directory. Run with: Rscript R/install_packages.R")
    }
  }
)

# Allow an explicit target library via first CLI argument (used by build scripts)
cli_args <- commandArgs(trailingOnly = TRUE)
if (length(cli_args) >= 1 && nzchar(cli_args[1])) {
  lib_dir <- cli_args[1]
} else {
  lib_dir <- file.path(script_dir, "library")
}
dir.create(lib_dir, showWarnings = FALSE, recursive = TRUE)

cat("Installing R packages into:", lib_dir, "\n")

required_pkgs <- c("tidyverse", "data.table", "mgcv", "emmeans")

install.packages(
  required_pkgs,
  lib = lib_dir,
  repos = "https://cloud.r-project.org",
  dependencies = TRUE
)

# Verify
missing <- required_pkgs[!sapply(required_pkgs, function(p) {
  requireNamespace(p, lib.loc = lib_dir, quietly = TRUE)
})]

if (length(missing) > 0) {
  cat("\nWARNING: Failed to install:", paste(missing, collapse = ", "), "\n")
} else {
  cat("\nAll packages installed successfully.\n")
}
