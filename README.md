# Analysis code for "SARS-CoV-2 variants associated with vaccine breakthrough in the Delaware Valley through summer 2021" 

To regenerate the analysis, run `source('analyzeLineages.R')` in R.

The code depends on the R package `rstan` available from CRAN (`install.packages('rstan')`) and was run in R version 4.1.1 with `rstan` version 2.21.2.

This code should generate analysis result files in the `out` directory. Example output is available in the `exampleOut` directory where:
 * countsPrediction.pdf:
 * sDropEnrichment.pdf:
 * vaccineEnrichment.pdf:
 * individualVariants.pdf:
 * strainFoldVaccineEnrichment.csv:
 * dropFoldVaccineEnrichment.csv:
 * variantProps.csv:



