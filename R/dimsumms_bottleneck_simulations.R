
#' dimsumms_bottleneck_simulations
#'
#' Simulate bottlenecks
#'
#' @param dataset_dir directory with pre-processed DMS data
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
dimsumms_bottleneck_simulations <- function(
  dataset_dir,
  execute = TRUE
){
  #Return if analysis not executed
  if(!execute){
    return()
  }

  #Display status
  message(paste("\n\n*******", "running stage: dimsumms_bottleneck_simulations", "*******\n\n"))

  #No bottleneck
  variant_data_merge <- fread(file.path(dataset_dir, "datasets", "BB_TARDBP_290_2017-06-13_variant_data_merge.tsv"))
  variant_data_merge <- variant_data_merge[apply(variant_data_merge[,.SD,.SDcols = names(variant_data_merge)[grep("_count$", names(variant_data_merge))]], 1, sum)!=0,]
  variant_data_merge <- variant_data_merge[,.SD,,.SDcols = grep("nt_seq|_count$", names(variant_data_merge))]
  names(variant_data_merge) <- sapply(strsplit(names(variant_data_merge), "_e"), '[', 1)
  write.table(variant_data_merge, file = file.path(dataset_dir, "datasets", "BB_TARDBP_290_2017-06-13_NoBottleneck_t0_variant_data_merge.tsv"), sep = "\t", quote = F, row.names = F)

  #Simulate DNA extraction bottleneck
  dimsumms__dna_extraction_bottleneck(
    input_file = file.path(dataset_dir, "datasets", "BB_TARDBP_290_2017-06-13_variant_data_merge.tsv"),
    outpath = file.path(dataset_dir, "datasets"))

  #Simulate library bottleneck
  dimsumms__library_bottleneck(
    input_file = file.path(dataset_dir, "datasets", "BB_TARDBP_290_2017-06-13_variant_data_merge.tsv"),
    outpath = file.path(dataset_dir, "datasets"))

  #Simulate replicate bottleneck
  dimsumms__replicate_bottleneck(
    input_file = file.path(dataset_dir, "datasets", "BB_TARDBP_290_2017-06-13_variant_data_merge.tsv"),
    outpath = file.path(dataset_dir, "datasets"))

  #Over-sequencing factor manipulations
  variant_data_merge <- fread(file.path(dataset_dir, "datasets", "BB_TARDBP_290_2017-06-13_NoBottleneck_t0_variant_data_merge.tsv"))
  temp <- variant_data_merge
  #Increase counts in replicate 1 input by factor of 3
  over_seq_factor <- 3
  variant_data_merge <- data.table::copy(temp)
  variant_data_merge[, input1 := input1*over_seq_factor]
  write.table(variant_data_merge, file = file.path(dataset_dir, "datasets", "BB_TARDBP_290_2017-06-13_NoBottleneck_t0_overseq3_variant_data_merge.tsv"), sep = "\t", quote = F, row.names = F)
  #Increase counts in replicate 1 input by factor of 10
  over_seq_factor <- 10
  variant_data_merge <- data.table::copy(temp)
  variant_data_merge[, input1 := input1*over_seq_factor]
  write.table(variant_data_merge, file = file.path(dataset_dir, "datasets", "BB_TARDBP_290_2017-06-13_NoBottleneck_t0_overseq10_variant_data_merge.tsv"), sep = "\t", quote = F, row.names = F)

}