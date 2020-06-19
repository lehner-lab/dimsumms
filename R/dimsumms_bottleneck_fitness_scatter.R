
#' dimsumms_bottleneck_fitness_scatter
#'
#' Plot fitness before/after bottleneck and different minimum input count thresholds
#'
#' @param fitness_list list of data.tables (required)
#' @param outpath output path for plots and saved objects (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
dimsumms_bottleneck_fitness_scatter <- function(
  fitness_list,
  outpath,
  execute = TRUE
  ){

	#Return if analysis not executed
	if(!execute){
		return()
	}

	#Display status
	message(paste("\n\n*******", "running stage: dimsumms_bottleneck_fitness_scatter", "*******\n\n"))

	#Create output directory
	dimsumms__create_dir(dimsumms_dir = outpath)

	### Fitness scatterplots 
	###########################

	#Merged data.table
	tox_dt <- do.call("rbind", fitness_list)

	#Subset to nonsynonymous variants
	tox_dt <- tox_dt[Nham_aa!=0,]

	#Variant codes
	tox_dt[, var_code := apply(.SD, 1, paste0, collapse = ":"),,.SDcols = c("mut_code", "mut_code1", "mut_code2")]

	# #Variant class
	# tox_dt[, var_class := 3]
	# for(b in c("library", "replicate", "extraction")){
	# 	temp_list <- as.list(tox_dt[bottleneck==b,table(var_code)])
	# 	tox_dt[bottleneck==b, var_class := as.numeric(unlist(temp_list[var_code]))]
	# }
	# #Recode variant class
	# tox_dt[var_class==1, filter_status := "fail"]
	# tox_dt[var_class==2, filter_status := "conditional"]
	# tox_dt[var_class==3, filter_status := "pass"]
	# tox_dt[, filter_status := factor(filter_status, levels = c("pass", "conditional", "fail"))]
	# #Remove duplicate variants for each bottleneck
	# tox_dt <- tox_dt[filter=="none",]

	#Fitness and error of double amino acid mutants
	tox_dt[Nham_aa==2,fitness := fitness_uncorr]
	tox_dt[Nham_aa==2,sigma := sigma_uncorr]

	#Convert bottleneck to factor
	tox_dt[, bottleneck := factor(bottleneck, levels = c("none", "library", "replicate", "extraction"))]

	#Comparative fitness scatterplots
	plot_dt_lb <- merge(tox_dt[bottleneck=="none"], tox_dt[bottleneck=="library"], by = "var_code")
	plot_dt_rb <- merge(tox_dt[bottleneck=="none"], tox_dt[bottleneck=="replicate"], by = "var_code")
	plot_dt_db <- merge(tox_dt[bottleneck=="none"], tox_dt[bottleneck=="extraction"], by = "var_code")
	plot_dt <- rbind(plot_dt_lb, plot_dt_rb, plot_dt_db)
  temp_cor <- plyr::ddply(plot_dt, c("bottleneck.y", "filter.y"), plyr::summarize, cor = round(cor(fitness.x, fitness.y, use = "pairwise.complete"), 2), n = length(fitness.x))
	d <- ggplot2::ggplot(plot_dt, ggplot2::aes(fitness.x, fitness.y)) +
    ggplot2::geom_hex() +
    ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
	  ggplot2::xlab("Fitness (no bottleneck)") +
	  ggplot2::ylab("Fitness (bottleneck)") +
	  ggplot2::geom_abline(linetype = 2) +
	  ggplot2::theme_bw() +
	  ggplot2::facet_grid(filter.y ~ bottleneck.y, scales="free") +
	  ggplot2::geom_text(data = temp_cor, ggplot2::aes(label = paste0("R = ", cor, "\nn = ", n)), x = 1, y = -3, size = 2)
	ggplot2::ggsave(file=file.path(outpath, 'fitness_scatterplots_bottleneck.pdf'), width=6, height=5)

}

