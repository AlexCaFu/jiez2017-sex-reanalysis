# Reanalysis code and results: source-resolved co-metabolism (JieZ_2017 cohort)

Companion deposit for: Cabrera-Fuentes HA. Source-resolved co-metabolism: locating where biological sex acts on microbiome–host pathways in atherosclerotic cardiovascular disease. mSystems. https://doi.org/10.1128/msystems.00445-26 (DOI assigned by the journal; not yet registered when this version was prepared).

Version 1.1.0, released 17 November 2026. Version 1.0.0 (2 July 2026): https://doi.org/10.5281/zenodo.21144754. All versions: https://doi.org/10.5281/zenodo.21144753.

## Changes from v1.0.0

- Final article title and citation.
- `02_results_export.R` writes every value of the reanalysis reported in the article to `results/`, with the full console output.
- PERMANOVA for all 380 subjects with recorded sex and, with body mass index as a covariate, for the 312 subjects with body mass index recorded.
- `set.seed(1)` is called immediately before each permutation test. The v1.0.0 script sets the seed once; run that way, the adjusted model gives P = 0.013 for sex and P = 0.023 for the sex-by-disease interaction instead of 0.010 and 0.021. R² values are identical; the difference is permutation error with 999 permutations.
- Per-species models report the number of observations used (312: `lm()` drops the 68 subjects without body mass index) and the coefficient estimates.
- The MetaCyc pathway screen adds the terms "cholate" and "dehydrox". The v1.0.0 terms did not match the bile acid 7α-dehydroxylation pathway (7ALPHADEHYDROX-PWY, "cholate degradation (bacteria, anaerobic)"), which is present in 62 of the 380 subjects.
- `install_env.R` names the scripts of this repository; v1.0.0 referred to `JieZ2017_sex_analysis.R`, which was never part of it.
- `CITATION.cff` and the checksum file `MANIFEST.md5` added.

## Data source

- Cohort: JieZ_2017; Jie Z et al. Nat Commun 2017;8:845, https://doi.org/10.1038/s41467-017-00900-1. Sequence data: European Nucleotide Archive PRJEB21528.
- Profiles: curatedMetagenomicData (Bioconductor; Pasolli E et al. Nat Methods 2017;14:1023–1024, https://doi.org/10.1038/nmeth.4468), resources 2021-10-14.JieZ_2017.relative_abundance (ExperimentHub EH7157) and 2021-10-14.JieZ_2017.pathway_abundance (EH7155). They are downloaded at run time (files of 0.3 MB and 12 MB, plus the ExperimentHub index) and are not redistributed here.
- Subjects: 385 profiled; 380 with recorded sex (controls: 101 women, 69 men; cases: 53 women, 157 men). Body mass index is recorded for 312 of the 380 (controls 163, cases 149). Medication is recorded for 149 of the 210 cases and for none of the controls.
- Species: 240 species present in at least 10% of the 380 samples. With `rownames = "short"`, `returnSamples()` drops one species that has no taxonomy entry in the package (Collinsella massiliensis, present in 60% of the samples) and prints a message saying so. Two species groups without a species-rank name in the package appear under genus labels in the results (genus:Enterobacter, the Enterobacter cloacae complex; genus:Streptococcus, the Streptococcus anginosus group).

## How to run (R console)

1. `source("install_env.R")` once.
2. `source("01_reanalysis_JieZ2017.R")` runs the analysis as archived in v1.0.0.
3. `source("02_results_export.R")` writes the result tables below.

## Results (`results/`)

| File | Content |
|---|---|
| 01_samples_by_group_and_sex.csv | Subjects by group and sex, with the number having body mass index and medication recorded |
| 02_species_models_bmi_adjusted.csv | 240 species: observations used, estimates, P and Benjamini–Hochberg FDR for the sex-by-disease interaction and for sex |
| 03_permanova_all380_unadjusted.csv | PERMANOVA on Aitchison distances, 999 permutations, 380 subjects, terms added sequentially |
| 04_permanova_bmi_subset.csv | The same model with body mass index, 312 subjects |
| 05_betadisper_sex.csv | Homogeneity of multivariate dispersion by sex |
| 06_pathway_screen.csv | Unstratified MetaCyc pathways matching the screen terms, with the number of samples in which each is nonzero |
| 07_sessionInfo.txt | R and package versions of the run that wrote these tables |
| console_log.txt | Console output of that run |

## Figures

Fig. 1 and Fig. 2 of the article are conceptual schematics that contain no values from this reanalysis; no code in this repository draws them. `Figure2_evidence_map.R` and `Figure3_resolution_cascade.R` are unchanged from v1.0.0. They draw two figures of the original submission (July 2026) that are not part of the published article and are kept for the record.

## Reproducibility check

The tables in `results/` were written on 17 November 2026 with R 4.6.1, vegan 2.7-5 and curatedMetagenomicData 3.20.0 (`results/07_sessionInfo.txt`). The same day the analysis was rerun in a second environment (R 4.3.3, vegan 2.6-4) from the same two ExperimentHub resources: every value agreed to within a relative difference of 1e-10. The July 2026 output (R 4.6.1, vegan 2.7-5; `sessionInfo.txt`) is reproduced exactly when `set.seed(1)` precedes each permutation test, as in `02_results_export.R`.

## License

MIT (see `LICENSE`), for the code and the result tables.
