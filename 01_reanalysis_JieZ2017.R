#!/usr/bin/env Rscript
###############################################################################
# Reanalysis of the JieZ_2017 atherosclerotic cardiovascular disease cohort
# Companion code for:
#   Cabrera-Fuentes HA. "Source-resolved co-metabolism: locating where
#   biological sex acts on microbiome-host pathways in atherosclerosis."
#   mSystems Perspective.
#
# Reproduces every statistic reported in the reanalysis paragraph and Methods:
#   - per-species disease x sex interaction models with BH-FDR
#   - PERMANOVA on Aitchison distances (999 permutations) + betadisper
#   - screen of MetaCyc pathway abundances for the four co-metabolic pathways
#
# Data source (NOT redistributed here): the JieZ_2017 shotgun-metagenomic
# profiles are pulled at run time from Bioconductor's ExperimentHub through
# curatedMetagenomicData. First run downloads and caches the profiles.
#
# Setup once:  Rscript install_env.R
# Run:         Rscript 01_reanalysis_JieZ2017.R
###############################################################################
suppressMessages({
  library(curatedMetagenomicData)
  library(SummarizedExperiment)
  library(dplyr)
  library(vegan)
})
set.seed(1)
getcol <- function(d, n) if (n %in% colnames(d)) as.character(d[[n]]) else rep("", nrow(d))

## ---- 1. retrieve species-level relative abundance -------------------------
studyRows <- sampleMetadata %>% filter(study_name == "JieZ_2017")
se <- returnSamples(studyRows, "relative_abundance", rownames = "short")
md <- as.data.frame(colData(se))

sexcol   <- intersect(c("gender", "sex"), colnames(md))[1]
md$sex   <- tolower(getcol(md, sexcol))
condstr  <- tolower(paste(getcol(md, "disease"), getcol(md, "study_condition")))
md$grp   <- ifelse(grepl("acvd|athero", condstr), "Atherosclerosis", "Healthy")

keep <- md$sex %in% c("male", "female")
se <- se[, keep]; md <- md[keep, , drop = FALSE]
md$sex <- factor(md$sex, levels = c("female", "male"))
md$grp <- factor(md$grp, levels = c("Healthy", "Atherosclerosis"))
cat("\n== realised sample sizes (sex is partially confounded with case status) ==\n")
print(table(disease = md$grp, sex = md$sex))

## ---- 2. feature filtering (10% prevalence) + CLR --------------------------
X  <- assay(se); X <- X[rowMeans(X > 0) >= 0.10, , drop = FALSE]
Xp <- X + (min(X[X > 0]) / 2)
gm <- apply(Xp, 2, function(z) exp(mean(log(z))))
clr <- t(log(sweep(Xp, 2, gm, "/")))          # samples x features
cat(sprintf("features after 10%% prevalence filter: %d\n", ncol(clr)))

## ---- 3. per-species model: disease * sex (+ BMI), BH-FDR ------------------
md$BMI  <- suppressWarnings(as.numeric(getcol(md, "BMI")))
has_bmi <- mean(!is.na(md$BMI)) > 0.8
form <- as.formula(paste0("y ~ grp * sex", if (has_bmi) " + BMI" else ""))
cat("per-species model:", deparse(form), "\n")
res <- lapply(colnames(clr), function(f) {
  d <- md; d$y <- clr[, f]
  fit <- try(lm(form, data = d), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  co <- summary(fit)$coefficients
  ix <- grep(":", rownames(co)); sx <- grep("^sexmale", rownames(co))
  data.frame(feature = f,
             p_interaction = if (length(ix)) co[ix[1], 4] else NA,
             p_sex         = if (length(sx)) co[sx[1], 4] else NA)
}) %>% bind_rows()
res$fdr_interaction <- p.adjust(res$p_interaction, "BH")
res$fdr_sex         <- p.adjust(res$p_sex,         "BH")
write.csv(res[order(res$p_interaction), ], "JieZ2017_sex_disease_associations.csv", row.names = FALSE)
cat(sprintf("species tested: %d | min interaction FDR: %.3f | min sex-main FDR: %.3f\n",
            nrow(res), min(res$fdr_interaction, na.rm = TRUE), min(res$fdr_sex, na.rm = TRUE)))

## ---- 4. community-level: PERMANOVA (Aitchison) + betadisper ---------------
D <- dist(clr)                                  # Euclidean of CLR = Aitchison
cat("\n== PERMANOVA (full sample, unadjusted), Aitchison distance, 999 perms ==\n")
print(adonis2(D ~ grp * sex, data = md, permutations = 999, by = "terms"))
cat("\n== betadisper (sex): homogeneity of multivariate dispersion ==\n")
print(anova(betadisper(D, md$sex)))

sub <- !is.na(md$BMI)
if (sum(sub) > 20) {
  cat("\n== PERMANOVA (BMI-adjusted subset) ==\n")
  print(adonis2(dist(clr[sub, , drop = FALSE]) ~ grp * sex + BMI,
                data = md[sub, ], permutations = 999, by = "terms"))
}

## ---- 5. screen MetaCyc pathway abundances for the four pathways -----------
pw <- returnSamples(studyRows, "pathway_abundance", rownames = "short")
terms <- c("choline", "carnitine", "betaine", "trimethylamine", "bile", "cholesterol")
hits  <- grep(paste(terms, collapse = "|"), tolower(rownames(pw)), value = TRUE)
cat("\n== MetaCyc pathways matching co-metabolic terms ==\n")
if (length(hits)) print(hits) else
  cat("none: standard MetaCyc profiles do not encode microbial TMA (cutC),\n",
      "bile-acid (bsh/bai), or cholesterol (ismA) conversion as named pathways.\n")

## ---- 6. capture environment for the archive ------------------------------
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
cat("\nDone. Wrote associations CSV and sessionInfo.txt.\n")
