#!/usr/bin/env Rscript
###############################################################################
# 02_results_export.R  (v1.1.0)
# Writes every value of the JieZ_2017 reanalysis reported in
#   Cabrera-Fuentes HA. Source-resolved co-metabolism: locating where biological
#   sex acts on microbiome-host pathways in atherosclerotic cardiovascular
#   disease. mSystems. https://doi.org/10.1128/msystems.00445-26
# to the folder results/:
#   01 samples by group and sex, with body mass index and medication records
#   02 per-species models (observations used, estimates, P, BH FDR)
#   03 PERMANOVA, all subjects with recorded sex, no covariate
#   04 PERMANOVA, subjects with recorded body mass index, BMI as covariate
#   05 betadisper by sex
#   06 MetaCyc pathway screen (terms include cholate and 7-alpha-dehydroxylation)
#   07 sessionInfo, and console_log.txt with the full console output
# Data handling and the CLR step are those of 01_reanalysis_JieZ2017.R (v1.0.0).
# set.seed(1) is called immediately before each permutation test, so each
# P value does not depend on the tests run before it.
# Run from the repository folder in the R console: source("02_results_export.R")
###############################################################################
dir.create("results", showWarnings = FALSE)
old <- setwd("results")
sink("console_log.txt", split = TRUE)
tryCatch({
  suppressMessages({
    library(curatedMetagenomicData)
    library(SummarizedExperiment)
    library(dplyr)
    library(vegan)
  })
  getcol <- function(d, n) if (n %in% colnames(d)) as.character(d[[n]]) else rep("", nrow(d))

  ## 1. samples (same rules as 01_reanalysis_JieZ2017.R)
  studyRows <- sampleMetadata %>% dplyr::filter(study_name == "JieZ_2017")
  cat("JieZ_2017 rows in sampleMetadata:", nrow(studyRows), "\n")
  se <- returnSamples(studyRows, "relative_abundance", rownames = "short")
  md <- as.data.frame(colData(se))
  cat("profiled samples:", ncol(se), "| species rows:", nrow(se), "\n")
  md$sex  <- tolower(getcol(md, intersect(c("gender", "sex"), colnames(md))[1]))
  condstr <- tolower(paste(getcol(md, "disease"), getcol(md, "study_condition")))
  md$grp  <- ifelse(grepl("acvd|athero", condstr), "Atherosclerosis", "Healthy")
  keep <- md$sex %in% c("male", "female")
  se <- se[, keep]; md <- md[keep, , drop = FALSE]
  md$sex <- factor(md$sex, levels = c("female", "male"))
  md$grp <- factor(md$grp, levels = c("Healthy", "Atherosclerosis"))
  md$BMI <- suppressWarnings(as.numeric(getcol(md, "BMI")))
  md$medication <- getcol(md, "treatment"); md$medication[md$medication %in% c("", "NA")] <- NA
  tab <- md %>% group_by(group = grp, sex) %>%
    summarise(subjects = n(), bmi_recorded = sum(!is.na(BMI)),
              medication_recorded = sum(!is.na(medication)), .groups = "drop")
  print(as.data.frame(tab)); write.csv(tab, "01_samples_by_group_and_sex.csv", row.names = FALSE)
  cat("with recorded sex:", nrow(md), "| with BMI:", sum(!is.na(md$BMI)),
      "| medication recorded:", sum(!is.na(md$medication)), "\n")

  ## 2. prevalence filter and CLR (identical to 01_reanalysis_JieZ2017.R)
  X  <- assay(se); X <- X[rowMeans(X > 0) >= 0.10, , drop = FALSE]
  Xp <- X + (min(X[X > 0]) / 2)
  gm <- apply(Xp, 2, function(z) exp(mean(log(z))))
  clr <- t(log(sweep(Xp, 2, gm, "/")))
  cat("species after 10% prevalence filter:", ncol(clr), "\n")

  ## 3. per-species models; lm() drops subjects without BMI, so nobs is reported
  res <- lapply(colnames(clr), function(f) {
    d <- md; d$y <- clr[, f]
    fit <- lm(y ~ grp * sex + BMI, data = d)
    co <- summary(fit)$coefficients
    ix <- grep(":", rownames(co))[1]
    data.frame(feature = f, nobs = nobs(fit),
               estimate_interaction = co[ix, 1], p_interaction = co[ix, 4],
               estimate_sex = co["sexmale", 1], p_sex = co["sexmale", 4])
  }) %>% bind_rows()
  res$fdr_interaction <- p.adjust(res$p_interaction, "BH")
  res$fdr_sex         <- p.adjust(res$p_sex, "BH")
  res <- res[order(res$p_interaction), ]
  write.csv(res, "02_species_models_bmi_adjusted.csv", row.names = FALSE)
  cat("species tested:", nrow(res), "| observations per model:", paste(unique(res$nobs), collapse = ", "),
      "| lowest interaction FDR:", round(min(res$fdr_interaction), 3),
      "| lowest sex FDR:", round(min(res$fdr_sex), 3), "\n")

  ## 4. PERMANOVA (Aitchison = Euclidean distance on CLR), both models, and betadisper
  set.seed(1)
  pA <- adonis2(dist(clr) ~ grp * sex, data = md, permutations = 999, by = "terms")
  print(pA); write.csv(as.data.frame(pA), "03_permanova_all380_unadjusted.csv")
  ok <- !is.na(md$BMI)
  set.seed(1)
  pB <- adonis2(dist(clr[ok, , drop = FALSE]) ~ grp * sex + BMI, data = md[ok, ], permutations = 999, by = "terms")
  print(pB); write.csv(as.data.frame(pB), "04_permanova_bmi_subset.csv")
  bd <- anova(betadisper(dist(clr), md$sex))
  print(bd); write.csv(as.data.frame(bd), "05_betadisper_sex.csv")

  ## 5. MetaCyc pathway screen on unstratified pathways
  pw <- returnSamples(studyRows, "pathway_abundance", rownames = "short")
  pw <- pw[, colnames(pw) %in% rownames(md)]
  nm <- rownames(pw)
  unstrat <- nm[!grepl("\\|", nm)]
  terms <- c("choline", "carnitine", "betaine", "trimethylamine", "bile", "cholate", "dehydrox", "cholesterol")
  hits <- unstrat[grepl(paste(terms, collapse = "|"), tolower(unstrat))]
  pr <- data.frame(pathway = hits,
                   samples_nonzero = sapply(hits, function(h) sum(assay(pw)[h, ] > 0)),
                   samples = ncol(pw))
  print(pr); write.csv(pr, "06_pathway_screen.csv", row.names = FALSE)

  writeLines(capture.output(sessionInfo()), "07_sessionInfo.txt")
  cat("Done. Result tables written to results/\n")
}, finally = { sink(); setwd(old) })
