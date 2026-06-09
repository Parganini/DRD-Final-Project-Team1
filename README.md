# DRD Final Project - Team 1

Final project for the DNA & RNA Dynamics course.

## Team assignment

- Group ID: 1
- Address: 42796479
- Detection p-value threshold: 0.05
- Normalization: preprocessNoob
- Differential methylation test: t-test

## Deliverables

- Rendered report: HTML or PDF
- Raw R code used for each pipeline step

## Project structure

```text
drd-final-project-team1/
├── data/
│   ├── raw/
│   │   └── .gitkeep
│   └── processed/
│       └── .gitkeep
├── scripts/
│   ├── 00_setup.R
│   ├── 01_load_raw_data.R
│   ├── 02_quality_control.R
│   ├── 03_raw_beta_m_values.R
│   ├── 04_normalization_noob.R
│   ├── 05_pca.R
│   ├── 06_differential_methylation_ttest.R
│   └── 07_final_plots.R
├── report/
│   └── team1_report.Rmd
├── figures/
│   └── .gitkeep
├── results/
│   └── .gitkeep
├── docs/
│   └── .gitkeep
├── .gitignore
└── README.md
```

## Notes

Raw IDAT files and large processed objects should not be committed to GitHub.
