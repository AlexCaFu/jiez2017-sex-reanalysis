#!/usr/bin/env Rscript
###############################################################################
# Environment setup for the mSystems sex-analysis scripts.
# Run ONCE on a machine WITH internet:   Rscript install_env.R
#
# Installs:
#   - ggplot2, dplyr, vegan      -> CRAN packages used by the scripts
#   - curatedMetagenomicData     -> needed by 01_reanalysis_JieZ2017.R and 02_results_export.R
#     (this is the remote resource; it pulls profiles from Bioconductor's
#      ExperimentHub over the network, it is not a local dataset)
#
# Figure2_evidence_map.R needs only ggplot2, so if you just want the figure,
# the first block below is sufficient.
###############################################################################

options(timeout = 1200, repos = c(CRAN = "https://cloud.r-project.org"))

## 1. CRAN packages (Figure 2 needs only these) -----------------------------
cran <- c("ggplot2", "dplyr", "vegan")
need <- setdiff(cran, rownames(installed.packages()))
if (length(need)) install.packages(need)

## 2. Bioconductor + curatedMetagenomicData (only for the JieZ_2017 analysis)-
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
bioc <- c("curatedMetagenomicData", "SummarizedExperiment")
need <- setdiff(bioc, rownames(installed.packages()))
if (length(need)) BiocManager::install(need, update = FALSE, ask = FALSE)

## 3. Connectivity handshake (a query, not the full analysis) ---------------
check <- function(label, expr) {
  cat(sprintf("%-34s", paste0(label, " ...")))
  ok <- tryCatch({ force(expr); TRUE },
                 error = function(e) { cat(" FAILED:", conditionMessage(e), "\n"); FALSE })
  if (ok) cat(" OK\n"); invisible(ok)
}

cat("\n--- package check ---\n")
check("ggplot2 + dplyr present",
      stopifnot(all(cran %in% rownames(installed.packages()))))

cat("\n--- Bioconductor / curatedMetagenomicData ---\n")
suppressMessages({
  check("curatedMetagenomicData loads", library(curatedMetagenomicData))
  check("sampleMetadata reachable", {
    m <- sampleMetadata
    stopifnot("JieZ_2017" %in% m$study_name)
    n <- sum(m$study_name == "JieZ_2017")
    cat(sprintf("(JieZ_2017 has %d samples in the index) ", n))
  })
})

cat("\nIf every row says OK:\n")
cat("   Rscript 01_reanalysis_JieZ2017.R    # analysis as archived in v1.0.0; downloads about 12 MB on first run\n")
cat("   Rscript 02_results_export.R         # writes the result tables to results/\n")
cat("   Rscript Figure2_evidence_map.R      # Figure 2 of the original submission, not of the article (no data)\n\n")
cat(R.version.string, "\n")
if (requireNamespace("BiocManager", quietly = TRUE)) print(BiocManager::version())
