
#' dimsumms_real_bottleneck_scatterplot_matrices
#'
#' Scatterplot matrices of Input and Output sample variant counts for doubles from real DMS datasets
#'
#' @param dataset_dir directory with preprocessed DMS datasets (preprocessed by dimsumms_errormodel_leaveoneout_preprocess_datasets.R script)
#' @param outpath output path for plots and saved objects (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
dimsumms_real_bottleneck_scatterplot_matrices <- function(
  dataset_dir,
  outpath,
  execute = TRUE
){
  #Return if analysis not executed
  if(!execute){
    return()
  }

  #Display status
  message(paste("\n\n*******", "running stage: dimsumms_real_bottleneck_scatterplot_matrices", "*******\n\n"))
  
  #Create output directory
  dimsumms__create_dir(dimsumms_dir = outpath)
  
  #FOS cis data
  dataset_name <- "GD_FOSintra_2016-06-17_variant_data_merge"
  variant_data_merge <- fread(file.path(dataset_dir, "datasets", paste0(dataset_name, ".tsv")))
  #Sample names
  input_samples_pattern = "_s0_"
  output_samples_pattern = "_s1_"
  input_samples <- colnames(variant_data_merge)[grep(input_samples_pattern, colnames(variant_data_merge))]
  output_samples <- colnames(variant_data_merge)[grep(output_samples_pattern, colnames(variant_data_merge))]

  #All-vs-all sample count correlations - only doubles (in single AA mutants)
  if(length(input_samples)!=0 | length(output_samples)!=0){
    temp_dt <- variant_data_merge[,.SD,.SDcols = c(input_samples, output_samples, "Nham_nt", "Nham_aa")]
    names(temp_dt)[grep("_count", names(temp_dt))] <- names(temp_dt)[grep("_count", names(temp_dt))]
    dimsumms__ggpairs_binhex(
      input_dt = log10(temp_dt[Nham_nt==2 & Nham_aa==1,.SD,.SDcols = c(input_samples, output_samples)]+1), 
      output_file_prefix = file.path(outpath, paste0(dataset_name, "_diagnostics_report_scatterplotmatrix_doubles")),
      xlab = "log10(variant count + 1)",
      ylab = "log10(variant count + 1)",
      size = 0.1)
  }

  #FOS-JUN trans data
  dataset_name <- "GD_FOSJUN_2015-12-15_variant_data_merge"
  variant_data_merge <- fread(file.path(dataset_dir, "datasets", paste0(dataset_name, ".tsv")))
  #Sample names
  input_samples_pattern = "_s0_"
  output_samples_pattern = "_s1_"
  input_samples <- colnames(variant_data_merge)[grep(input_samples_pattern, colnames(variant_data_merge))]
  output_samples <- colnames(variant_data_merge)[grep(output_samples_pattern, colnames(variant_data_merge))]

  #All-vs-all sample count correlations - only doubles (in single AA mutants)
  if(length(input_samples)!=0 | length(output_samples)!=0){
    temp_dt <- variant_data_merge[,.SD,.SDcols = c(input_samples, output_samples, "Nham_nt", "Nham_aa")]
    names(temp_dt)[grep("_count", names(temp_dt))] <- names(temp_dt)[grep("_count", names(temp_dt))]
    dimsumms__ggpairs_binhex(
      input_dt = log10(temp_dt[Nham_nt==2 & Nham_aa==1,.SD,.SDcols = c(input_samples, output_samples)]+1), 
      output_file_prefix = file.path(outpath, paste0(dataset_name, "_diagnostics_report_scatterplotmatrix_doubles")),
      xlab = "log10(variant count + 1)",
      ylab = "log10(variant count + 1)",
      size = 0.1)
  }

  #tRNA 23C data
  dataset_name <- "Li2018_tRNA_sel23_variant_data_merge"
  variant_data_merge <- fread(file.path(dataset_dir, "datasets", paste0(dataset_name, ".tsv")))
  #Sample names
  input_samples_pattern = "_s0_"
  output_samples_pattern = "_s1_"
  input_samples <- colnames(variant_data_merge)[grep(input_samples_pattern, colnames(variant_data_merge))]
  output_samples <- colnames(variant_data_merge)[grep(output_samples_pattern, colnames(variant_data_merge))]

  #All-vs-all sample count correlations - only doubles
  if(length(input_samples)!=0 | length(output_samples)!=0){
    temp_dt <- variant_data_merge[,.SD,.SDcols = c(input_samples, output_samples, "Nham_nt")]
    names(temp_dt)[grep("_count", names(temp_dt))] <- names(temp_dt)[grep("_count", names(temp_dt))]
    dimsumms__ggpairs_binhex(
      input_dt = log10(temp_dt[Nham_nt==2,.SD,.SDcols = c(input_samples, output_samples)]+1), 
      output_file_prefix = file.path(outpath, paste0(dataset_name, "_diagnostics_report_scatterplotmatrix_doubles")),
      xlab = "log10(variant count + 1)",
      ylab = "log10(variant count + 1)",
      size = 0.1)
  }


}