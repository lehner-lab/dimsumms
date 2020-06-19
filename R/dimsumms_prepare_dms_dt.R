
#' dimsumms_prepare_dms_dt
#'
#' Load and reformat fitness data.table.
#'
#' @param fitness_list named list of folder paths with fitness estimates (required)
#'
#' @return Data.table with fitness estimates
#' @export
#' @import data.table
dimsumms_prepare_dms_dt <- function(
  fitness_list
  ){
	
	### Load data
	###########################a

	#DMS variant data 
	dms_dt1 <- dimsumms__load_toxicity(fitness_list[[1]], read_filter=F)

	#Combine regions
	dms_dt1[, region := names(fitness_list)[1]]
	dms_dt <- dms_dt1
	#Remove missing fitness values
	if("fitness1" %in% names(dms_dt) & "fitness2" %in% names(dms_dt)){
		dms_dt <- dms_dt[!is.na(fitness) | (!is.na(fitness1) & !is.na(fitness2)),]
	}
	#Absolute position (singles)
	dms_dt[, Pos_abs := Pos+as.numeric(region)-1]
	#Absolute position (doubles)
	dms_dt[, Pos_abs1 := Pos1+as.numeric(region)-1]
	dms_dt[, Pos_abs2 := Pos2+as.numeric(region)-1]
	#Mutation code (singles)
	dms_dt[, mut_code := paste0(WT_AA, Pos_abs, Mut)]
	#Mutation code (doubles)
	dms_dt[, mut_code1 := paste0(WT_AA1, Pos_abs1, Mut1)]
	dms_dt[, mut_code2 := paste0(WT_AA2, Pos_abs2, Mut2)]

	#Return data.table
	return(dms_dt)

}

