#' dimsumms_errormodel_fit
#'
#' estimate multiplicative input/output and additive error terms from DMS dataset
#'
#' @param DT DMS dataset data.table
#' @param reps replicates for which the error model parameters should be estimated
#' @param Ncores how many cores to use for parallelization
#' @param Nbootstraps # bootstraps for error model calculation
#' @param maxN how many variants are used per bootstrap
#' @param max_tries_per_fit how many restarts of the fitting routine in case nls doesn't yield results
#' @param lower_bound_add lower limited for additive error terms
#' @param Fcorr fitness normalization parameters
#'
#' @return Nothing
#' @export
#' @import data.table
#' @import foreach
dimsumms_errormodel_fit = function(
  DT,
  reps,
  Ncores=1,
  Nbootstraps=100,
  maxN = 10000,
  max_tries_per_fit = 20,
  lower_bound_add = 10^-4,
  Fcorr= c()) {
  
  reps_num = as.numeric(strsplit(reps,"")[[1]])
  
  work_data = copy(DT)
  ######## estimate replicate & over-sequencing error model
  #calcualte density of data along mean count based error (for density dependent weighting)
  # bins = 50
  # work_data[,mean_cbe := rowMeans(.SD),.SDcols = grep(paste0("^cbe[",reps,"]*$"),names(work_data))]
  # error_range = seq(work_data[input_above_threshold == T & all_reads ==T,log10(quantile(mean_cbe^2,probs=0.001))],0,length.out = bins)
  # D = density(x = work_data[input_above_threshold == T & all_reads ==T,log10(mean_cbe^2)],
  #             from = work_data[,log10(quantile(mean_cbe^2,probs=0.001))],to=0,n=bins)
  # work_data[,bin_error := findInterval(mean_cbe^2,vec = 10^error_range)]
  # work_data[,bin_error_density := D$y[bin_error],bin_error]
  # work_data[bin_error == 0,bin_error_density := work_data[bin_error >0][which.min(bin_error),unique(bin_error_density)]]
  
  work_data[, bin_error_density := sqrt(max(.N, sqrt(nrow(work_data)))), Nmut]

  #set up fitting parameters
  Nreps = nchar(reps)
  doMC::registerDoMC(cores=Ncores)
  
  # calculate all combinations of replicates of length >= 2
  idx = list()
  for (i in (nchar(reps)):2) {
    idx = c(idx,combn(nchar(reps),i,function(x) list(x)))
  }
  
  # fitting routine
  t=proc.time()
  parameters = foreach(m=1:Nbootstraps,.combine = rbind) %dopar% {
    y = NA
    counter=0
    while (!is.list(y) & counter < max_tries_per_fit) { #this is a routine to catch failed fitting attempts
      tryCatch({
        #create data concatenation for all combinations of replicates
        bs_data = work_data[input_above_threshold == T & all_reads ==T & Nmut > 0][sample(.N,min(.N,maxN),replace = T)]
        F_data_list = list() #fitness data
        E_data_list = list() #count based error data (for weighting of datapoints)
        C_data_list = list() #count data
        NR = list() # simple counter for how many replicates are in each data row
        for (i in 1:length(idx)) {
          F_data_list[[i]] = bs_data[,.SD,.SDcols = grep(paste0("^fitness[",paste0(reps_num[idx[[i]]],collapse=""),"]$"),names(work_data))]
          E_data_list[[i]] = bs_data[,.SD,.SDcols = grep(paste0("^cbe[",paste0(reps_num[idx[[i]]],collapse=""),"]*$"),names(work_data))]
          C_data_list[[i]] = bs_data[,1/.SD,.SDcols = grep(paste0("put[",paste0(reps_num[idx[[i]]],collapse=""),"]*$"),names(work_data))]
          NR[[i]] = data.table(rep(length(idx[[i]]),nrow(bs_data)))
        }
        F_data = as.matrix(rbindlist(F_data_list,fill=T)) #make matrices with #replicates = #columns
        E_data = as.matrix(rbindlist(E_data_list,fill=T))
        C_data = as.matrix(rbindlist(C_data_list,fill=T))
        NRT = as.matrix(rbindlist(NR))
        C_data[is.na(C_data)] = 0 #avoid NA for fitting
        if (length(Fcorr) > 0) { #correct count terms for fitness normalization factors
          C_data = C_data * matrix(rep(Fcorr,2),nrow = nrow(C_data),ncol=Nreps*2,byrow=T)
        }
        V_data = apply(F_data,1,var,na.rm=T) #calculate variance from fitness data
        Dw = (bs_data$bin_error_density + quantile(bs_data$bin_error_density,na.rm=T,probs=0.001)) #density based weighting of data points
        Ew = rowMeans(E_data,na.rm=T)^2 #count based error weighting of data points (needed because deviation is proportional to expcectation)
        
        #binary variable to avoid replicate error fitting for replicates that are not present in a data row
        BV = F_data
        BV[!is.na(BV)] = 1
        BV[is.na(BV)] = 0
        
        #fit formula
        fit_formula = as.formula(paste0("V_data ~ (",paste0("BV[,",1:Nreps,"] * (p[",1:Nreps,"] * C_data[,",1:Nreps,"] +  p[",Nreps + 1:Nreps,"] * C_data[,",Nreps + 1:Nreps,"] + p[",2*Nreps+1:Nreps,"])",collapse = " + "),") / NRT"))
        #fitting
        y=nls(formula = fit_formula,
              start = list(p = c(rep(1,each=Nreps*2)+10^rnorm(2*Nreps),rep(0,each=Nreps)+0.01*10^rnorm(1*Nreps))),
              lower = c(rep(1,each=Nreps*2),rep(lower_bound_add,each=Nreps)),
              weights = c((NRT-1)/ Ew / Dw),
              algorithm = "port")
      },error=function(cond){})
      counter = counter + 1}
    #if fit was succesful after less than max_tries_per_fit attempts, transfer parameters to p
    if (counter < max_tries_per_fit) {
      p=summary(y)$parameter[,1]
    } else { #otherwise move on to next fit
      p=rep(NA,3*Nreps)
    }
    return(p)
  }
  print(proc.time() -t)
  return(parameters)
}
