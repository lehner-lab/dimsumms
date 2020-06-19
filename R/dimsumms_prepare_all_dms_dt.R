
#' dimsumms_prepare_all_dms_dt
#'
#' Load and reformat a list of fitness data.tables.
#'
#' @param fitness_dir Base directory with subfolders containing results from different DiMSum runs
#'
#' @return Nothing
#' @export
#' @import data.table
dimsumms_prepare_all_dms_dt <- function(
  fitness_dir
	){
	
	### Load data
	###########################

	dataset_names <- list(
		"NoBottleneck" = "none",
		"LibraryBottleneck_03" = "library",
		"ReplicateBottleneck_03" = "replicate",
		"DnaExtractionBottleneck_03" = "extraction")

	dt_list <- list()
	for(i in names(dataset_names)){
		#No filter
		dt_list[[i]] <- dimsumms_prepare_dms_dt(fitness_list = list("290" = file.path(fitness_dir, i)))
		dt_list[[i]][,bottleneck := dataset_names[[i]]]
		dt_list[[i]][,filter := "none"]
		#Bottleneck datasets
		if(i != "NoBottleneck"){
			#Any filter			
			dt_list[[paste0(i, "_tany")]] <- dimsumms_prepare_dms_dt(fitness_list = list("290" = file.path(fitness_dir, paste0(i, "_tany"))))
			dt_list[[paste0(i, "_tany")]][,bottleneck := dataset_names[[i]]]
			dt_list[[paste0(i, "_tany")]][,filter := "any"]
			#All filter			
			dt_list[[paste0(i, "_tall")]] <- dimsumms_prepare_dms_dt(fitness_list = list("290" = file.path(fitness_dir, paste0(i, "_tall"))))
			dt_list[[paste0(i, "_tall")]][,bottleneck := dataset_names[[i]]]
			dt_list[[paste0(i, "_tall")]][,filter := "all"]
		}
	}

	#Return data.table
	return(dt_list)

}

