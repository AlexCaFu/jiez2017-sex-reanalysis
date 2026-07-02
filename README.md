# Reanalysis code: source-resolved co-metabolism in atherosclerosis (JieZ_2017)

Companion code for the mSystems Perspective:

> Cabrera-Fuentes HA. *Source-resolved co-metabolism: locating where biological
> sex acts on microbiome-host pathways in atherosclerosis.* mSystems (submitted).

This repository reproduces the independent reanalysis of the **JieZ_2017**
atherosclerotic cardiovascular disease cohort reported in the manuscript.

## What it reproduces
- Per-species disease-by-sex interaction models on CLR-transformed abundances,
  with Benjamini-Hochberg false-discovery control.
- Community-level **PERMANOVA** on Aitchison distances (999 permutations) and a
  homogeneity-of-dispersion check (**betadisper**).
- A screen of MetaCyc **pathway abundances** for choline, carnitine, betaine,
  trimethylamine, bile-acid, and cholesterol terms.

## Data availability
No primary data are redistributed here. The JieZ_2017 shotgun-metagenomic
species and pathway profiles are retrieved at run time from Bioconductor's
ExperimentHub through the `curatedMetagenomicData` package. Original cohort:
Jie Z, et al. *Nat Commun* 2017;8:845. doi:10.1038/s41467-017-00900-1.
`curatedMetagenomicData`: Pasolli E, et al. *Nat Methods* 2017;14:1023-1024.
doi:10.1038/nmeth.4468.

## How to run
```r
Rscript install_env.R            # one-time: installs CRAN + Bioconductor deps
Rscript 01_reanalysis_JieZ2017.R # downloads/caches profiles, runs the analysis
```
The first run downloads several hundred MB of profiles and caches them locally.
Outputs: `JieZ2017_sex_disease_associations.csv` and `sessionInfo.txt`.

## Figures
- `Figure2_evidence_map.R` - the source-resolution evidence map (Figure 2).
- `Figure3_resolution_cascade.R` - the resolution-cascade figure (supplementary).

## Environment
Developed under R 4.6.x. A machine-specific `sessionInfo.txt` is written by the
analysis script; deposit the one produced on the machine used for the paper.

## License
Code released under the MIT License (see `LICENSE`).
