# Overview

Welcome to the GitHub repository for the following publication: "DiMSum: A pipeline for analyzing deep mutational scanning (DMS) data and diagnosing common experimental pathologies"

Here you'll find an R package with all scripts to reproduce the figures and results from the computational analyses described in the paper.

# Required Software

To run the dimsumms pipeline you will need the following software and associated packages:

* **[R](https://www.r-project.org/) >=v3.5.2** (data.table, ggplot2, hexbin, plyr, reshape2, RColorBrewer)

The following packages are optional:

* **[DiMSum](https://github.com/lehner-lab/DiMSum)** (pipeline for pre-processing deep mutational scanning data i.e. FASTQ to counts)

# Installation and loading

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

Fitness scores, pre-processed data and required miscellaneous files should be downloaded from [here](https://www.dropbox.com/s/cjfn3c5rby5h0mb/misc.zip?dl=0) to your project directory (see 'base_dir' argument) i.e. where output files should be written, and unzipped.

# Running

There are a number of options available for running the dimsumms pipeline depending on user requirements.

* ## Basic (default)

Default pipeline functionality uses fitness scores (see 'Required Data') to reproduce all figures in the publication. The **[DiMSum](https://github.com/lehner-lab/DiMSum)** package is not required for this default functionality.

* ## Raw read processing

Raw read processing is not handled by the dimsumms pipeline. FastQ files from deep mutational scanning (DMS) experiments were processed using **[DiMSum](https://github.com/lehner-lab/DiMSum)** itself.

DiMSum command-line arguments and Experimental design files required to obtain fitness scores and variant counts from FastQ files are available [here]().

# Pipeline

The top-level function **dimsumms()** is the recommended entry point to the pipeline and reproduces the figures and results from the computational analyses described in the following publication: "DiMSum: A pipeline for analyzing deep mutational scanning (DMS) data and diagnosing common experimental pathologies" (Faure AJ and Schmiedel JM et al.). See section on "Required Data" above for instructions on how to obtain all required data and miscellaneous files before running the pipeline.

## Stage 1: Fitness comparisons before/after bottlenecks and filtering

This stage ('dimsumms_bottleneck_fitness_scatter') produces scatterplots of fitness scores before/after bottlenecks and soft/hard filtering based on minimum Input variant count thresholds.

## Stage 2: Biological conclusion (fitness v.s. hydrophobicity) comparisons before/after bottlenecks and filtering

This stage ('dimsumms_bottleneck_fitness_vs_hydrophobicity_scatter') produces scatterplots of fitness v.s. hydrophobicity before/after bottlenecks and soft/hard filtering based on minimum Input variant count thresholds.

## Stage 3: Error model results before/after replicate error and over-sequencing factor manipulations

This stage ('dimsumms_error_model_data_manipulations') produces plots of error model parameters before/after manipulating replicate error and over-sequencing factors for replicate 1.

## Stage 4: Error model leave-one-out cross validation benchmark

This stage ('dimsumms_errormodel_leaveoneout') performs leave-one-out cross validation on published DMS datasets to benchmark error model performance. The preprocessed datasets are supplied with the required data .zip-file, but the preprocessing script ('dimsumms_errormodel_leaveoneout_preprocess_datasets.R') can be run separately on the raw data to reproduce the full workflow.


