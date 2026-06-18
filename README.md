# DRD Team 1 Final Report

Final project for DNA & RNA Dynamics.

## Main files

- DRD_Team1_FinalReport.R: raw R workflow scaffold
- DRD_Team1_FinalReport.Rmd: final report source
- data/raw/: raw input data, including IDAT files and SampleSheet_Report_II.csv
- results/rds/: intermediate R objects
- results/tables/: exported tables
- results/figures/: exported figures

## Team 1 assignment

- Group ID: 1
- Address: 42796479
- Detection p-value threshold: 0.05
- Normalization: preprocessNoob
- Differential methylation test: t-test

## Current status

The workflow files currently contain the Team 1 constants and the 12 required section headings. The analysis code still needs to be implemented step by step.

## How to run

From the repository root, run:

`
source("DRD_Team1_FinalReport.R")
`

To render the report, open DRD_Team1_FinalReport.Rmd in RStudio and knit it.
