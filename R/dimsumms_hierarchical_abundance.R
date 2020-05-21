
#' dimsumms_hierarchical_abundance
#'
#' Plots showing hierarchical abundance of DMS experiments
#'
#' @param outpath output path for plots and saved objects (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
dimsumms_hierarchical_abundance <- function(
  outpath,
  execute = TRUE
  ){

	#Return if analysis not executed
	if(!execute){
		return()
	}

	#Display status
	message(paste("\n\n*******", "running stage: dimsumms_hierarchical_abundance", "*******\n\n"))

	#Create output directory
	dimsumms__create_dir(dimsumms_dir = outpath)

	### 
	###########################


	nmut <- seq(0,2)
	seq_length <- 100
	construct_mutations_pervar <- 1
	mut_per_pos <- 3
	seqerror_rate <- 10^-3
	total_reads <- 1e7

	# whats the frequency with which variants with a specific number of mutations appear in the dataset?
#	mut_freq <- dpois(x = nmut, lambda = construct_mutations_pervar)
	mut_freq <- c(1/3, 1/3, 1/3)

	# how many possible variants exist with the same number of mutations?
	# poss_vars_permut = factorial(seq_length) / factorial(seq_length - nmut) / factorial(nmut) * (mut_per_pos ^ nmut)
	poss_vars_permut <- c(1, seq_length * mut_per_pos, seq_length * (seq_length - 1) / 2 * mut_per_pos ^ 2)
	# obs_vars_permut = poss_vars_permut * c(1, 0.5, 0.25)


	# how many reads does each individual variant have?
	reads_permut <- mut_freq / poss_vars_permut
	# ...relative to wildtype?
	reads_permut_relwt = reads_permut / reads_permut[1]


	#likelihood of mutations per sequencing read (0 to 2 mutations)
	x_down <- dpois(x = nmut, lambda = seq_length * seqerror_rate* mut_per_pos / (mut_per_pos + 1))

	# how many variants flow from wild-type to singles or doubles?
	wildtype_to_others = reads_permut_relwt[1] * (1 - x_down[1]) / poss_vars_permut[1] #how many reads flow away from wild-type
	wildtype_to_others[2:3] = reads_permut_relwt[1] * x_down[2:3] / poss_vars_permut[2:3] #to each single or each double variant

	# how many variants flow from singles to wild-type or doubles?
	single_to_others = mut_freq[2] / mut_freq[1] * seqerror_rate / (mut_per_pos + 1) #how many reads are added to wild-type by ALL singles
	single_to_others[2] = reads_permut_relwt[2] * (1 - x_down[1]) / poss_vars_permut[1] #how many reads flow away PER single variant
	single_to_others[3] = reads_permut_relwt[2] * x_down[2] / poss_vars_permut[2] * 2 #how many reads are added to each double variant from singles

	doubles_to_others = mut_freq[3] / mut_freq[1] * (seqerror_rate / (mut_per_pos + 1))^2 #how many reads are added to wild-type by ALL doubles
	doubles_to_others[2] = mut_freq[3] / mut_freq[1] * 2 * seqerror_rate / (mut_per_pos + 1) / poss_vars_permut[2] #how many reads are added to each single by ALL doubles
	doubles_to_others[3] = reads_permut_relwt[3] * (1 - x_down[1]) / poss_vars_permut[1] #how many reads flow away PER double variant

	DT <- data.table(nmut, mut_freq, 
	                poss_vars_permut, reads_permut_relwt, 
	                wildtype_to_others,
	                single_to_others,
	                doubles_to_others)
	DT[, nmut := factor(nmut)]
	#first plot showing frequency per mutation category
	ggplot2::ggplot(DT, ggplot2::aes(nmut, mut_freq, fill = nmut)) + 
	  ggplot2::geom_bar(stat = "identity") +
	#  labs(x = "# mutations", y = "frequency of all variant with X mutations")
	  ggplot2::labs(x = "# mutations", y = "frequency") +
	  ggplot2::theme_classic()
	ggplot2::ggsave(file=file.path(outpath, 'hierarchical_abundance_1.pdf'), width=4, height=3)

	#second plot showing number of possible variants
	ggplot2::ggplot(DT, ggplot2::aes(nmut, poss_vars_permut+1, fill = nmut)) + 
	  ggplot2::geom_bar(stat="identity") + 
	  ggplot2::scale_y_log10() +
	  ggplot2::labs(x = "# mutations", y = "log10(possible variants+1)") +
	  ggplot2::theme_classic()
	ggplot2::ggsave(file=file.path(outpath, 'hierarchical_abundance_2.pdf'), width=4, height=3)

	#third plot showing reads per individual variant
	ggplot2::ggplot(DT, ggplot2::aes(nmut, reads_permut_relwt*1e7*DT[nmut==0,mut_freq]/DT[,sum(mut_freq)], fill = nmut)) + 
	  ggplot2::geom_bar(stat="identity") + 
	  ggplot2::scale_y_log10() +
	  ggplot2::labs(x = "# mutations", y = "log10(reads per variant)") +
	  ggplot2::theme_classic()
	ggplot2::ggsave(file=file.path(outpath, 'hierarchical_abundance_3.pdf'), width=4, height=3)

}

