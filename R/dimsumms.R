
#' dimsumms
#'
#' Main analysis script.
#'
#' @param startStage Start at a specified analysis stage (default:1)
#' @param stopStage Stop at a specified analysis stage (default:0 i.e. no stop condition)
#' @param base_dir Base directory for all output file (default:NB private CRG server path; change accordingly)
#'
#' @return Nothing
#' @export
dimsumms <- function(
  startStage=1,
  stopStage=0,
  base_dir = "/users/blehner/afaure/DMS/Results/dimsumms_proj"
  ){

	colour_scheme <- list(
		"shade 0" = list(
			"#F4270C",
			"#F4AD0C",
			"#1B38A6",
			"#09B636"),
		"shade 1" = list(
			"#FFB0A5",
			"#FFE4A5",
			"#9DACE3",
			"#97E9AD"),
		"shade 2" = list(
			"#FF6A56",
			"#FFCB56",
			"#4C63B7",
			"#43C766"),
		"shade 3" = list(
			"#A31300",
			"#A37200",
			"#0C226F",
			"#007A20"),
		"shade 4" = list(
			"#410800",
			"#412D00",
			"#020B2C",
			"#00300D"))

  #First and last analysis stages
  first_stage <- startStage
  last_stage <- stopStage

	### Fitness distributions and fitness vs hydrophobicity scatterplots
	###########################

	#Load fitness data.tables
	fitness_list <- dimsumms_prepare_all_dms_dt(
		fitness_dir = file.path(base_dir, "misc", "DiMSum_fitness"))

	#Fitness comparisons before/after bottlenecks and filtering
	stagenum <- 1
	toxicity_aaprop_dt_nb <- dimsumms_bottleneck_fitness_scatter(
		fitness_list = fitness_list,
		outpath = dimsumms__format_dir(dir_suffix="_dimsumms_bottlenck_fitness_scatter", stagenum=stagenum, base_dir=base_dir),
		execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)))

	#Biological conclusion (fitness v.s. hydrophobicity) comparisons before/after bottlenecks and filtering
	stagenum <- 2
	toxicity_aaprop_dt_nb <- dimsumms_bottleneck_fitness_vs_hydrophobicity_scatter(
		fitness_list = fitness_list,
		outpath = dimsumms__format_dir(dir_suffix="_dimsumms_bottleneck_fitness_vs_hydrophobicity_scatter", stagenum=stagenum, base_dir=base_dir),
		aaprop_file = file.path(base_dir, "misc", "amino_acid_properties", "amino_acid_properties_annotated_supplementary.txt"),
		aaprop_file_selected = file.path(base_dir, "misc", "amino_acid_properties", "selected.amino_acid_properties.txt"),
		colour_scheme = colour_scheme,
		execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)))

	#Error model results before/after replicate error and over-sequencing factor manipulations
	stagenum <- 3
	dimsumms_error_model_data_manipulations(
		error_model_dir = file.path(base_dir, "misc", "DiMSum_errormodel"),
		outpath = dimsumms__format_dir(dir_suffix="_dimsumms_error_model_data_manipulations", stagenum=stagenum, base_dir=base_dir),
		execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)))

	#Error model performance plots
	stagenum <- 4
	dimsumms_error_model_performance(
		performance_dir = file.path(base_dir, "misc", "Performance"),
		outpath = dimsumms__format_dir(dir_suffix="_dimsumms_error_model_performance", stagenum=stagenum, base_dir=base_dir),
		execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)))

}



