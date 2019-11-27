
#' dimsumms_bottleneck_fitness_vs_hydrophobicity_scatter
#'
#' Plot fitness before/after bottleneck and different minimum input count thresholds
#'
#' @param fitness_list list of data.tables (required)
#' @param outpath output path for plots and saved objects (required)
#' @param aaprop_file path to amino acid properties file (required)
#' @param aaprop_file_selected path to file with selected subset of identifiers
#' @param colour_scheme colour scheme file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
dimsumms_bottleneck_fitness_vs_hydrophobicity_scatter <- function(
  fitness_list,
  outpath,
  aaprop_file,
  aaprop_file_selected,
  colour_scheme,
  execute = TRUE
  ){

	#Return if analysis not executed
	if(!execute){
		return()
	}

	#Display status
	message(paste("\n\n*******", "running stage: dimsumms_bottleneck_fitness_vs_hydrophobicity_scatter", "*******\n\n"))

	#Create output directory
	dimsumms__create_dir(dimsumms_dir = outpath)

	### Bottleneck biological conclusions comparisons
	###########################

	#Calculate single and double mutant effects from AA PCA (using single mutants)
	for(i in seq_along(names(fitness_list))){
		fitness_list[[i]] <- dimsumms_aa_properties_mutant_effects(
			toxicity_dt = fitness_list[[i]],
			outpath = outpath,
			aaprop_file = aaprop_file,
			aaprop_file_selected = aaprop_file_selected,
			colour_scheme = colour_scheme,
			report = (i==1))
		#Set fitness = fitness_uncorr if fitness_uncorr!=NA
		fitness_list[[i]][!is.na(fitness_uncorr), fitness := fitness_uncorr]
		#Set sigma = sigma_uncorr if fitness_uncorr!=NA
		fitness_list[[i]][!is.na(fitness_uncorr), sigma := sigma_uncorr]
		#Remove missing fitness values and STOPs
		fitness_list[[i]] <- fitness_list[[i]][!is.na(fitness) & Nham_aa %in% c(1,2) & STOP==FALSE,]
	}

	# #Combine
	# toxicity_aaprop_dt_nb_tany <- data.table::copy(toxicity_aaprop_dt_nb)
	# toxicity_aaprop_dt_nb_tall <- data.table::copy(toxicity_aaprop_dt_nb)

	# #No bottleneck
	# toxicity_aaprop_dt_nb[, filter := "none"]
	# toxicity_aaprop_dt_nb[, bottleneck := "none"]
	# toxicity_aaprop_dt_nb_tany[, filter := "any"]
	# toxicity_aaprop_dt_nb_tany[, bottleneck := "none"]
	# toxicity_aaprop_dt_nb_tall[, filter := "all"]
	# toxicity_aaprop_dt_nb_tall[, bottleneck := "none"]
	# #Library bottleneck
	# toxicity_aaprop_dt_lb[, filter := "none"]
	# toxicity_aaprop_dt_lb[, bottleneck := "library"]
	# toxicity_aaprop_dt_lb_tany[, filter := "any"]
	# toxicity_aaprop_dt_lb_tany[, bottleneck := "library"]
	# toxicity_aaprop_dt_lb_tall[, filter := "all"]
	# toxicity_aaprop_dt_lb_tall[, bottleneck := "library"]
	# #Replicate bottleneck
	# toxicity_aaprop_dt_rb[, filter := "none"]
	# toxicity_aaprop_dt_rb[, bottleneck := "replicate"]
	# toxicity_aaprop_dt_rb_tany[, filter := "any"]
	# toxicity_aaprop_dt_rb_tany[, bottleneck := "replicate"]
	# toxicity_aaprop_dt_rb_tall[, filter := "all"]
	# toxicity_aaprop_dt_rb_tall[, bottleneck := "replicate"]
	# #DNA extraction bottleneck
	# toxicity_aaprop_dt_db[, filter := "none"]
	# toxicity_aaprop_dt_db[, bottleneck := "extraction"]
	# toxicity_aaprop_dt_db_tany[, filter := "any"]
	# toxicity_aaprop_dt_db_tany[, bottleneck := "extraction"]
	
	#Merged data.table
	tox_dt <- do.call("rbind", fitness_list)

	#Scatter plots - singles and doubles in hotspot
	feat_dt <- tox_dt[Nham_aa %in% c(1,2) & STOP==FALSE,.SD,,.SDcols = c("fitness", "PC1 (Hydrophobicity)", "filter", "bottleneck")]
	plot_df <- as.data.frame(feat_dt)
	colnames(plot_df)[2] <- "hydrophobicity"
	plot_df[,"bottleneck"] <- factor(plot_df[,"bottleneck"], levels = c("none", "library", "replicate", "extraction"))
	plot_df[,"filter"] <- factor(plot_df[,"filter"], levels = c("none", "all", "any"))
	#Remove NAs
	plot_df <- plot_df[!is.na(plot_df[,"fitness"]),]
  #plot_colours = c(colour_scheme[["shade 0"]][[1]], colour_scheme[["shade 0"]][[2]], colour_scheme[["shade 0"]][[3]])
  plot_colours = dimsumms__gg_color_hue(2)
  temp_cor <- plyr::ddply(plot_df, c("filter", "bottleneck"), plyr::summarize, cor = round(cor(hydrophobicity, fitness, use = "pairwise.complete"), 2), n = length(hydrophobicity))
	d <- ggplot2::ggplot(plot_df, ggplot2::aes(hydrophobicity, fitness)) +
    ggplot2::stat_binhex(bins=50) +
	  ggplot2::xlab("Hydrophobicity") +
	  ggplot2::ylab("Fitness") +
	  ggplot2::theme_bw() +
	  ggplot2::facet_grid(bottleneck ~ filter, scales="free") +
	  ggplot2::geom_smooth(method = "lm", linetype = 2, se = F, color = "black") +
	  # ggplot2::coord_cartesian(ylim=c(-0.4, 0.4)) +
	  ggplot2::geom_text(data = temp_cor, ggplot2::aes(label = paste0("R = ", cor, "\nn = ", n)), x = 0.25, y = -0.3, size = 2) +
	  ggplot2::scale_fill_gradientn(colours=plot_colours,name = "Frequency", na.value=plot_colours[length(plot_colours)], limits=c(0, 20))
	ggplot2::ggsave(file=file.path(outpath, 'fitness_vs_hydrophobicity_scatter_all.pdf'), width=5, height=6)

	#Scatter plots - singles and doubles in hotspot
	feat_dt <- tox_dt[Nham_aa %in% c(1,2) & STOP==FALSE & (Pos>=23 | (Pos1>=23 & Pos2>=23)),.SD,,.SDcols = c("fitness", "PC1 (Hydrophobicity)", "filter", "bottleneck")]
	plot_df <- as.data.frame(feat_dt)
	colnames(plot_df)[2] <- "hydrophobicity"
	plot_df[,"bottleneck"] <- factor(plot_df[,"bottleneck"], levels = c("none", "library", "replicate", "extraction"))
	plot_df[,"filter"] <- factor(plot_df[,"filter"], levels = c("none", "all", "any"))
	#Remove NAs
	plot_df <- plot_df[!is.na(plot_df[,"fitness"]),]
  #plot_colours = c(colour_scheme[["shade 0"]][[1]], colour_scheme[["shade 0"]][[2]], colour_scheme[["shade 0"]][[3]])
  plot_colours = dimsumms__gg_color_hue(2)
  temp_cor <- plyr::ddply(plot_df, c("filter", "bottleneck"), plyr::summarize, cor = round(cor(hydrophobicity, fitness, use = "pairwise.complete"), 2), n = length(hydrophobicity))
	d <- ggplot2::ggplot(plot_df, ggplot2::aes(hydrophobicity, fitness)) +
    ggplot2::stat_binhex(bins=50) +
	  ggplot2::xlab("Hydrophobicity") +
	  ggplot2::ylab("Fitness") +
	  ggplot2::theme_bw() +
	  ggplot2::facet_grid(bottleneck ~ filter, scales="free") +
	  ggplot2::geom_smooth(method = "lm", linetype = 2, se = F, color = "black") +
	  # ggplot2::coord_cartesian(ylim=c(-0.4, 0.4)) +
	  ggplot2::geom_text(data = temp_cor, ggplot2::aes(label = paste0("R = ", cor, "\nn = ", n)), x = 0.25, y = -0.3, size = 2) +
	  ggplot2::scale_fill_gradientn(colours=plot_colours,name = "Frequency", na.value=plot_colours[length(plot_colours)], limits=c(0, 20))
	ggplot2::ggsave(file=file.path(outpath, 'fitness_vs_hydrophobicity_scatter_all_hotspot.pdf'), width=6, height=6)

}

