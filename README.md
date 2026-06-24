# DRD Final Project - Team 1

DNA methylation analysis of CTRL and DIS samples using Illumina HumanMethylation450K arrays.

This repository contains the Team 1 final project for the DNA & RNA Dynamics course at Alma Mater Studiorum - Universita di Bologna.

## Repository Contents

- `DRD_Team1_FinalReport.Rmd`: final report source.
- `DRD_Team1_FinalReport.R`: standalone R workflow.
- `DRD_Team1_FinalReport.html`: rendered report.
- `data/raw/`: raw input data, including IDAT files and `SampleSheet_Report_II.csv`.
- `results/rds/`: intermediate R objects.
- `results/tables/`: exported result tables.
- `results/figures/`: exported figures.
- `docs/course_scripts/`: course reference scripts.

## Workflow

The analysis follows the assigned 12-step pipeline: raw data import, fluorescence inspection, quality control, raw beta/M-value summaries, `preprocessNoob` normalization, PCA, probe-wise t-tests, multiple testing correction, volcano and Manhattan plots, and heatmap visualization.

## Team 1 Assigned Parameters

- Group ID: 1
- Probe address: `42796479`
- Detection p-value threshold: `0.05`
- Normalization method: `preprocessNoob`
- Differential methylation test: `t-test`

## Reproducibility

To run the standalone workflow from the repository root:

```r
source("DRD_Team1_FinalReport.R")
```

To render the report:

```r
rmarkdown::render("DRD_Team1_FinalReport.Rmd")
```

## Notes

Raw data and generated analysis outputs are intentionally tracked in this repository for project reproducibility.
