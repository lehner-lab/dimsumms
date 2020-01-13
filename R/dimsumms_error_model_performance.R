
#' dimsumms_error_model_performance
#'
#' Plot error model performance results
#'
#' @param performance_dir Base directory with error model performance results
#' @param outpath output path for plots and saved objects (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
dimsumms_error_model_performance <- function(
  performance_dir,
  outpath,
  execute = TRUE
  ){

	#Return if analysis not executed
	if(!execute){
		return()
	}

	#Display status
	message(paste("\n\n*******", "running stage: dimsumms_error_model_performance", "*******\n\n"))

	#Create output directory
	dimsumms__create_dir(dimsumms_dir = outpath)

	### Error model performance plots
	###########################

	#Load data
	load(file.path(performance_dir, "all_datasets.RData"))

	#look at errors across datasets
	for (i in c(1,2,3,4,5,6,9,12,13,14)) {
	  X <- data.table::fread(file.path(performance_dir, paste0("validation_",dataset_label[i],"_error_estimates.txt")))
	  X[,dataset := dataset_label[i]]
	  
	  if (i == 1) {
	    all_datasets <- X
	  } else {
	    all_datasets <- rbind(all_datasets,X)
	  }
	}
	adm <- all_datasets[test_rep != "avg" & method != "MioAtrue",.(value = 2^mean((log2(sd_val)))),.(Nmut,method,dataset)]

	dataset_dict <- list(
		"TDP43_290" = "TDP-43 (290-331); Bolognesi et al. 2019",
		"TDP43_290_134" = "TDP-43 (290-331) only reps 1,3,4",
		"TDP43_332" = "TDP-43 (332-373); Bolognesi et al. 2019",
		"FOSJUN" = "FOS-JUN deepPCA; Diss et al. 2018",
		"FOScis" = "FOS deepPCA; Diss et al. 2018",
		"GRB2_GPD_epPCR" = "GRB2 deepPCA; in prep.",
		"GRB2_CYC_epPCR" = "GRB2 stabilityPCA; in prep.",
		"GB1" = "GB1; Olson et al. 2014",
		"tRNA_Li2018_s23" = "tRNA 23C; Li & Zhang 2018",
		"tRNA_Li2018_s30" = "tRNA 30C; Li & Zhang 2018",
		"tRNA_Li2018_s37" = "tRNA 37C; Li & Zhang 2018",
		"tRNA_Li2018_sDMSO" = "tRNA DMSO; Li & Zhang 2018")
	adm[,dataset := unlist(dataset_dict[dataset])]
	adm[,method := factor(method,levels=c("naive","cbe","ire","MioA","A","Mio"))]
	levels(adm$method) <- c("s.d.-based","count-based","CB + variant corr.","DiMSum full","DiMSum rep.","DiMSum mult.")
	adm[,Nmut := factor(Nmut,levels=c("1","2","all"))]
	levels(adm$Nmut) <- c("single mutants","double mutants","all variants")
	p <- ggplot2::ggplot(adm[!method %in% c("DiMSum rep.", "DiMSum mult.") & Nmut=="all variants",],ggplot2::aes(x = method, y = 1/(value))) +
	  ggplot2::geom_boxplot(ggplot2::aes(group = method), outlier.shape = NA) +
	  ggplot2::geom_point(ggplot2::aes(color = dataset), position = ggplot2::position_dodge(width = 0.2), size = 3) +
	  # scale_y_continuous(breaks = c(0,1,2,4,8)) +
	  # ggplot2::facet_wrap(~Nmut) +
	  # scale_color_brewer(palette="Set1") +
	  ggplot2::geom_hline(yintercept = 1,lty = 2) +
	  # ggplot2::geom_vline(xintercept = 4.5,lty = 3) +
	  	  ggplot2::labs(title = "average deviation of error estimates", y = "estimated / observed error", x = "") +
	  ggplot2::theme_bw() +
	  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
	ggplot2::ggsave(plot = p, filename = file.path(outpath, "validation_error_estimates_across_datasets_reps_fig3.pdf"), width = 5, height = 5)

	p <- ggplot2::ggplot(adm,ggplot2::aes(x = method, y = 1/(value))) +
	  ggplot2::geom_boxplot(ggplot2::aes(group = method), outlier.shape = NA) +
	  ggplot2::geom_point(ggplot2::aes(color = dataset), position = ggplot2::position_dodge(width = 0.2), size = 3) +
	  # scale_y_continuous(breaks = c(0,1,2,4,8)) +
	  ggplot2::facet_wrap(~Nmut) +
	  # scale_color_brewer(palette="Set1") +
	  ggplot2::geom_hline(yintercept = 1,lty = 2) +
	  ggplot2::geom_vline(xintercept = 4.5,lty = 3) +
	  	  ggplot2::labs(title = "average deviation of error estimates", y = "estimated / observed error", x = "") +
	  ggplot2::theme_bw() +
	  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
	ggplot2::ggsave(plot = p, filename = file.path(outpath, "validation_error_estimates_across_datasets_reps_supp.pdf"), width = 10, height = 5)

}

