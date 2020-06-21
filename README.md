Welcome to the GitHub repository for the following publication: "DiMSum: An error model and pipeline for analyzing deep mutational scanning (DMS) data and diagnosing common experimental pathologies"

Here you'll find an R package with all scripts to reproduce the figures and results from the computational analyses described in the paper. In addition we provide scripts and files to pre-process raw read datasets using DiMSum.

# Table Of Contents

* **1. [Required Software](#required-software)**
* **2. [Installation Instructions](#installation-instructions)**
* **3. [Required Data](#required-data)**
* **4. [Variant Count Data](#variant-count-data)**
* **5. [Pipeline Modes](#pipeline-modes)**
* **6. [Pipeline Stages](#pipeline-stages)**

# Required Software

To run the dimsumms pipeline you will need the following software and associated packages:

* **[_R_](https://www.r-project.org/) >=v3.5.2** (Biostrings, caTools, corpcor, cowplot, data.table, gdata, ggplot2, GGally, hexbin, lemon, optparse, parallel, pdist, plyr, ppcor, raster, reshape2, Rpdb, RColorBrewer, foreach, doMC)

The following packages are optional:

* **[_DiMSum_](https://github.com/lehner-lab/DiMSum)** (softawre for pre-processing deep mutational scanning data)

# Installation Instructions

Open R and enter:

```
# Install
if(!require(devtools)) install.packages("devtools")
devtools::install_github("lehner-lab/dimsumms")

# Load
library(dimsumms)

# Help
?dimsumms
```

# Required Data

Fitness scores, pre-processed data and required miscellaneous files should be downloaded from [here](https://www.dropbox.com/s/7yf0xcxyrklor6x/misc.zip?dl=0) to your project directory (see '_base_dir_' option) i.e. where output files should be written, and unzipped.

# Variant Count Data

If processing of variant count data is required ('_rerun_raw_' = T), files should be downloaded from [here](https://www.dropbox.com/s/b7gtqhtya6ay9ba/datasets.zip?dl=0), unzipped and saved within the [Required Data](#required-data) directory.

# Pipeline Modes

There are a number of options available for running the dimsumms pipeline depending on user requirements.

* ## Basic (default)

Default pipeline functionality ('_rerun_raw_' = F) uses fitness scores (see [Required Data](#required-data)) to reproduce all figures in the publication. The **[DiMSum](https://github.com/lehner-lab/DiMSum)** package is not required for this default functionality.

* ## Variant count processing

To preprocess variant count data (see [Variant Count Data](#variant-count-data)), the '_rerun_raw_' option should be set to TRUE ('_rerun_raw_' = T). Variant count data is optional for Stage 4 and required for Stages 6 and 7.

* ## Raw read processing

Raw read processing is not handled by the dimsumms pipeline. FastQ files from deep mutational scanning (DMS) experiments were processed using **[DiMSum](https://github.com/lehner-lab/DiMSum)** itself. DiMSum command-line arguments and Experimental design files required to obtain fitness scores and variant counts from FastQ files are included within the subdirectory 'inputfiles' of the [Required Data](#required-data).

# Pipeline Stages

The top-level function **dimsumms()** is the recommended entry point to the pipeline and reproduces the figures and results from the computational analyses described in the following publication: "DiMSum: An error model and pipeline for analyzing deep mutational scanning (DMS) data and diagnosing common experimental pathologies" (Faure AJ and Schmiedel JM et al.). See see [Required Data](#required-data) for instructions on how to obtain all required data and miscellaneous files before running the pipeline.

## Stage 1: Fitness comparisons before/after bottlenecks and filtering

This stage ('dimsumms_bottleneck_fitness_scatter') produces scatterplots of fitness scores before/after bottlenecks and soft/hard filtering based on minimum Input variant count thresholds.

## Stage 2: Biological conclusion (fitness v.s. hydrophobicity) comparisons before/after bottlenecks and filtering

This stage ('dimsumms_bottleneck_fitness_vs_hydrophobicity_scatter') produces scatterplots of fitness v.s. hydrophobicity before/after bottlenecks and soft/hard filtering based on minimum Input variant count thresholds.

## Stage 3: Error model results before/after replicate error and over-sequencing factor manipulations

This stage ('dimsumms_error_model_data_manipulations') produces plots of error model parameters before/after manipulating replicate error and over-sequencing factors for replicate 1.

## Stage 4: Error model leave-one-out cross validation benchmark

This stage ('dimsumms_errormodel_leaveoneout') performs leave-one-out cross validation on published DMS datasets to benchmark error model performance. The preprocessed datasets are supplied with the [Required Data](#required-data). To reproduce the full workflow, download the [Variant Count Data](#variant-count-data) and set the '_rerun_raw_' option to TRUE.

## Stage 5: Hierarchical abundance of DMS experiments

This stage ('dimsumms_hierarchical_abundance') produces the cartoon in Figure 1B.

## Stage 6: Scatterplot matrices of Input and Output sample variant counts for doubles from real DMS datasets

This stage ('dimsumms_real_bottleneck_scatterplot_matrices') produces the scatterplots shown in Supplementary Figure 9. To run this stage, download the [Variant Count Data](#variant-count-data) and set the '_rerun_raw_' option to TRUE.

## Stage 7: Bottleneck simulations

This stage ('dimsumms_bottleneck_simulations') simulates DNA extraction, library and replicate bottlenecks in DMS data. To run this stage, download the [Variant Count Data](#variant-count-data) and set the '_rerun_raw_' option to TRUE.


