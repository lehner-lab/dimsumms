
#' dimsumms_errormodel_data_manipulations
#'
#' Plot error model parameters before/after replicate error and over-sequencing factor manipulations in TDP43 290 library
#'
#' @param error_model_dir Base directory with subfolders containing results from different DiMSum runs
#' @param outpath output path for plots and saved objects (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
dimsumms_errormodel_data_manipulations <- function(
  error_model_dir,
  outpath,
  execute = TRUE
  ){

	#Return if analysis not executed
	if(!execute){
		return()
	}

	#Display status
	message(paste("\n\n*******", "running stage: dimsumms_error_model_data_manipulations", "*******\n\n"))

	#Create output directory
	dimsumms__create_dir(dimsumms_dir = outpath)

	### Error model results before/after replicate error and over-sequencing factor manipulations
	###########################

	#Error model plots
	error_model_list <- list()
	error_model_list[["original"]] <- fread(file.path(error_model_dir, "NoBottleneck", "errormodel.txt"))
	for(i in c("03", "10")){
		error_model_list[[paste0("overseq", i)]] <- fread(file.path(error_model_dir, paste0("NoBottleneck_overseq", i), "errormodel.txt"))
		error_model_list[[paste0("reperror", i)]] <- fread(file.path(error_model_dir, paste0("NoBottleneck_reperror", i), "errormodel.txt"))
	}
	#Combine
	plot_error_model <- do.call("rbind", error_model_list)
	#Dataset name
	plot_error_model[, dataset := rep(names(error_model_list), each = dim(error_model_list[[1]])[1])]

  #Add columns for upper/lower bounds of parameter estimates (for plotting)
  plot_error_model[,upper := mean_value + sd_value]
  plot_error_model[,lower := mean_value - sd_value]
  plot_error_model[lower < 0,lower := mean_value]
  # print(plot_error_model)
  
  #Plot1: input and output over-sequencing factor parameters +- sd
  a <- ggplot2::ggplot(plot_error_model[parameter %in% c("input","output") & dataset %in% c("original", "overseq03", "overseq10")],
    ggplot2::aes(x=interaction(parameter, rep), mean_value, ymin = lower, ymax = upper, color = factor(rep))) +
    ggplot2::geom_pointrange() +
    # ggplot2::scale_y_log10(limits = c(min(c(1, plot_error_model[parameter %in% c("input","output"), mean_value], plot_error_model[parameter %in% c("input", "output"), lower])),
    #   max(c(2.5, plot_error_model[parameter %in% c("input","output"), mean_value], plot_error_model[parameter %in% c("input", "output"), upper])))) +
    ggplot2::scale_y_log10() + ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::labs(y = "Over-sequencing factor", x = "Replicate (input or output)") +
    ggplot2::facet_grid(~dataset)
  ggplot2::ggsave(file.path(outpath, "errormodel_overseq.pdf"), a, width = 6, height = 4)

  #Plot2: replicate error parameters +- sd
  b <- ggplot2::ggplot(plot_error_model[parameter == "reperror" & !dataset %in% c("overseq03", "overseq10")], 
  	ggplot2::aes(x = as.factor(rep), y = sqrt(mean_value), ymin = sqrt(lower), ymax = sqrt(upper), color = factor(rep))) +
    ggplot2::geom_pointrange() +
    # ggplot2::scale_y_log10(limits = c(min(plot_error_model[parameter == "reperror", sqrt(lower)]),
    #   max(c(0.1, plot_error_model[parameter == "reperror", sqrt(upper)])))) +
    ggplot2::scale_y_log10() + ggplot2::theme_bw() + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::labs(y = "Replicate error", x = "Replicate") +
    ggplot2::facet_grid(~dataset)
  ggplot2::ggsave(file.path(outpath, "errormodel_reperror.pdf"), b, width = 5, height = 4)

}

