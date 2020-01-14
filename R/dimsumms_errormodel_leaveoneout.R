#' dimsumms_errormodel_leaveoneout
#'
#' Perform leave one out cross validation to benchmark error models
#'
#' @param dataset_dir directory with preprocessed DMS datasets (preprocessed by dimsumms_errormodel_leaveoneout_preprocess_datasets.R script)
#' @param Ncores number of cores to use to parallelize up analyes
#' @param outpath output path for plots and saved objects (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
dimsumms_errormodel_leaveoneout <- function(
  dataset_dir,
  Ncores = 1,
  outpath,
  execute = TRUE
){
  #Return if analysis not executed
  if(!execute){
    return()
  }
  
  #Display status
  message(paste("\n\n*******", "running stage: dimsumms_errormodel_leaveoneout", "*******\n\n"))
  
  #Create output directory
  dimsumms__create_dir(dimsumms_dir = outpath)
  
  #fetch files
  files = list.files(file.path(dataset_dir))
  dataset_label = sapply(1:length(files),function(X){strsplit(files[X],"\\.")[[1]][1]})
  
  ###############################
  ###  calculate error model ####
  ###############################
  
  for (dataset_idx in seq_along(files)) {
    
    print(paste0("running leave one out cross validation for ",dataset_label[dataset_idx]))
    
    #load dataset
    work_data = fread(file.path(dataset_dir,files[dataset_idx]))
    
    #how many replicates?
    reps = paste0(gsub("input","",grep("input",names(work_data),value=T)),collapse="")
    reps_num = as.numeric(strsplit(reps,"")[[1]])
    
    #####################
    ###### fitness ######
    #####################
    ## calculate fitness
    for (j in reps_num) {
      wt_corr = as.numeric(work_data[WT==T,log(.SD[,2]/.SD[,1]),,
                                     .SDcols = c(grep(paste0("^input",j,"$"),names(work_data)),grep(paste0("^output",j,"$"),names(work_data)))])
      
      work_data[,paste0("fitness",j) := log(.SD[,2]/.SD[,1]) - wt_corr,,
                .SDcols = c(grep(paste0("^input",j,"$"),names(work_data)),grep(paste0("^output",j,"$"),names(work_data)))]
    }
    
    # flag variants that don't have reads in all input/output replicates
    work_data[,all_reads := rowSums(.SD > 0) == (2*nchar(reps)),,.SDcols = grep(paste0("put[",reps,"]$"),names(work_data))]
    
    # find input read threshold for full fitness range
    input_count_threshold = work_data[all_reads == T,exp(-quantile(.SD,probs = 0.01,na.rm=T)),,.SDcols = grep("fitness",names(work_data))]
    
    #define variants above threshold for later use
    work_data[,input_above_threshold := rowSums(.SD > input_count_threshold) == nchar(reps),,.SDcols = grep(paste0("input[",reps,"]$"),names(work_data))]
    
    # normalize fitness by scaling and shifting
    set.seed(1603)
    minF = function(p) {
      F_norm = (F_data + matrix(p[(nchar(reps)+1):(2*nchar(reps))],nrow=nrow(F_data),ncol=nchar(reps),byrow = T)) * matrix(p[1:nchar(reps)],nrow=nrow(F_data),ncol=nchar(reps),byrow = T)
      F_avg = rowMeans(F_data + matrix(p[(nchar(reps)+1):(2*nchar(reps))],nrow=nrow(F_data),ncol=nchar(reps),byrow = T))
      diffF = sqrt(rowSums((F_norm - F_avg)^2))
      return(sum(diffF))
    }
    F_data = work_data[input_above_threshold == T & all_reads ==T,as.matrix(.SD),.SDcols = grep(paste0("fitness[",reps,"]$"),names(work_data))]
    W_data = work_data[input_above_threshold == T & all_reads ==T,as.matrix(1/.SD^2),.SDcols = grep(paste0("cbe[",reps,"]$"),names(work_data))]
    x=nlm(minF,rep(c(1,0),each=nchar(reps)))
    # print(x)
    p = x$estimate
    p[1:nchar(reps)] = p[1:nchar(reps)]/p[1]
    
    fitness_norm_model = data.table(t(p))
    names(fitness_norm_model) = c(paste0("scale_",reps_num),paste0("shift_",reps_num))
    
    #wild-type correction such that mean(wild-type) = 0
    wt_corr = work_data[WT ==T,rowMeans((.SD + unlist(fitness_norm_model[,.SD,,.SDcols = grep(paste0("shift_[",reps,"]"),names(fitness_norm_model))])) * 
                                          unlist(fitness_norm_model[,.SD,,.SDcols = grep(paste0("scale_[",reps,"]"),names(fitness_norm_model))])),
                        ,.SDcols = grep(paste0("fitness[",reps,"]"),names(work_data))]
    #normalize fitness values
    for (j in seq_along(reps_num)) {
      work_data[all_reads ==T,paste0("fitness",reps_num[j]) := (.SD + unlist(fitness_norm_model[,.SD,,.SDcols = paste0("shift_",reps_num[j])])) * unlist(fitness_norm_model[,.SD,,.SDcols = paste0("scale_",reps_num[j])]) - wt_corr,
                ,.SDcols = paste0("fitness",reps_num[j])]
    }
    
    #calculate count-based error for each variant and each replicate
    for (j in as.numeric(strsplit(reps,"")[[1]])) {
      wt_corr = as.numeric(work_data[WT==T,1/.SD[,2] + 1/.SD[,1],,
                                     .SDcols = c(grep(paste0("^input",j,"$"),names(work_data)),grep(paste0("^output",j,"$"),names(work_data)))])
      
      work_data[,paste0("cbe",j) := sqrt(unlist(fitness_norm_model[,.SD,.SDcols=paste0("scale_",j)])) * sqrt(1/.SD[,2] + 1/.SD[,1] + wt_corr),,
                .SDcols = c(grep(paste0("^input",j,"$"),names(work_data)),grep(paste0("^output",j,"$"),names(work_data)))]
    }
    
    ######################################
    ### leave one out cross validation ###
    ######################################
    # calculate error model on training replicates to estimate fitness in leftout replicate
    for (i in seq_along(reps_num)) {
      
      training_reps = paste0(strsplit(reps,"")[[1]][-i],collapse="")
      training_reps_num = as.numeric(strsplit(training_reps,"")[[1]])
      NTreps = nchar(training_reps)
      test_rep = strsplit(reps,"")[[1]][i]
      training_data = copy(work_data)
      
      #########################
      ###### error models #####
      #########################
      
      ##############################
      ### naive SD based error model
      training_data[,fitness_naive := rowMeans(.SD[,1:NTreps],na.rm=T),
                    ,.SDcols = grep(paste0("^fitness[", training_reps, "]$"),names(training_data))]
      training_data[,error_naive := sqrt(rowMeans((rowMeans(.SD[,1:NTreps],na.rm=T) - .SD[,1:NTreps])^2,na.rm=T)) / sqrt(NTreps),
                    ,.SDcols = grep(paste0("^fitness[", training_reps, "]$"),names(training_data))]
      
      ###########################
      ### count based error model
      training_data[,fitness_cbe := rowSums(.SD[,1:NTreps]/(.SD[,(NTreps+1):(2*NTreps)]^2),na.rm=T) / 
                      rowSums(1/(.SD[,(NTreps+1):(2*NTreps)]^2),na.rm=T),
                    ,.SDcols = c(grep(paste0("^fitness[", training_reps, "]$"),names(training_data)),
                                 grep(paste0("^cbe[", training_reps, "]$"),names(training_data)))]
      
      training_data[,error_cbe := sqrt(1 / rowSums(1/(.SD[,1:NTreps]^2 ),na.rm=T)),
                    ,.SDcols = grep(paste0("^cbe[", training_reps, "]$"),names(training_data))]
      
      ###########################################################################
      ### error model with replicate specific additive and multiplicative errors
      
      #note fitness replication factor
      Fcorr = unlist(fitness_norm_model[,.SD,.SDcols=grep(paste0("scale_[",training_reps,"]"),names(fitness_norm_model))])
      
      #estimate additve and multiplicative factors
      print(paste0("with replicate ",test_rep," held out"))
      parameters = dimsumms_errormodel_fit(DT = training_data,
                                                     reps = training_reps,
                                                     Ncores = Ncores,
                                                     Fcorr = Fcorr)
      #  table with error model parameters
      error_model_cycle = data.table(parameter = rep(c("input","output","reperror"),each=NTreps),
                                     rep = rep(training_reps_num,3),
                                     mean_value = colMeans(parameters,na.rm=T),
                                     sd_value = apply(parameters,2,sd,na.rm=T),
                                     ensemble = sum(!is.na(parameters[,1])),
                                     which_reps = training_reps)
      
      if (i == 1) {
        error_model = copy(error_model_cycle)
      } else {
        error_model = rbind(error_model,copy(error_model_cycle))
      }
      
      # calculate upper and lower error for parameters
      error_model_cycle[,upper := mean_value + sd_value]
      error_model_cycle[,lower := mean_value - sd_value]
      error_model_cycle[lower < 0,lower := mean_value]
      
      #######
      #fitness and estimated error from full error model (MioA - Multiplicative Input/Output terms + additive error terms)
      for (j in seq_along(training_reps_num)) {
        Corr = matrix(unlist(fitness_norm_model[,.SD,.SDcols=paste0("scale_",training_reps_num[j])]),ncol = 1,nrow=nrow(training_data))
        training_data[,paste0("error_MioA",training_reps_num[j]) := sqrt(Corr * rowSums(matrix(unlist(error_model_cycle[parameter %in% c("input","output") & rep == training_reps_num[j],mean_value]),nrow=.N,ncol=2,byrow = T)/.SD) + 
                                                                           matrix(error_model_cycle[parameter %in% c("reperror") & rep == training_reps_num[j],mean_value],nrow=.N,ncol=1,byrow = T)),,
                      .SDcols = c(grep(paste0("^input",training_reps_num[j],"$"),names(training_data)),
                                  grep(paste0("^output",training_reps_num[j],"$"),names(training_data)))]
      }
      #merged fitness values
      training_data[,fitness_MioA := rowSums(.SD[,1:NTreps]/(.SD[,(NTreps+1):(2*NTreps)]^2),na.rm=T) / 
                      rowSums(1/(.SD[,(NTreps+1):(2*NTreps)]^2),na.rm=T),
                    ,.SDcols = c(grep(paste0("^fitness[", training_reps, "]$"),names(training_data)),
                                 grep(paste0("^error_MioA[", training_reps, "]$"),names(training_data)))]
      #and merged error
      training_data[,error_MioA := sqrt(1/rowSums(1/.SD^2)),,
                    .SDcols = grep(paste0("^error_MioA[", training_reps, "]$"),names(training_data))]
      
      ########
      #fitness and estimated error when ONLY using multiplicative input/output terms 
      for (j in seq_along(training_reps_num)) {
        Corr = matrix(unlist(fitness_norm_model[,.SD,.SDcols=paste0("scale_",training_reps_num[j])]),ncol = 1,nrow=nrow(training_data))
        
        training_data[,paste0("error_Mio",training_reps_num[j]) := sqrt(Corr * rowSums(matrix(unlist(error_model_cycle[parameter %in% c("input","output") & rep == training_reps_num[j],mean_value]),nrow=.N,ncol=2,byrow = T)/.SD)),,
                      .SDcols = c(grep(paste0("^input",training_reps_num[j],"$"),names(training_data)),
                                  grep(paste0("^output",training_reps_num[j],"$"),names(training_data)))]
      }
      #merged fitness values
      training_data[,fitness_Mio := rowSums(.SD[,1:NTreps]/(.SD[,(NTreps+1):(2*NTreps)]^2),na.rm=T) / 
                      rowSums(1/(.SD[,(NTreps+1):(2*NTreps)]^2),na.rm=T),
                    ,.SDcols = c(grep(paste0("^fitness[", training_reps, "]$"),names(training_data)),
                                 grep(paste0("^error_Mio[", training_reps, "]$"),names(training_data)))]
      #and merged error
      training_data[,error_Mio := sqrt(1/rowSums(1/.SD^2)),,
                    .SDcols = grep(paste0("^error_Mio[", training_reps, "]$"),names(training_data))]
      
      ######
      ##fitness and estimated error when ONLY using additive terms 
      for (j in seq_along(training_reps_num)) {
        Corr = matrix(unlist(fitness_norm_model[,.SD,.SDcols=paste0("scale_",training_reps_num[j])]),ncol = 1,nrow=nrow(training_data))
        training_data[,paste0("error_A",training_reps_num[j]) := sqrt(Corr * rowSums(1/.SD) + 
                                                                        matrix(error_model_cycle[parameter %in% c("reperror") & rep == training_reps_num[j],mean_value],nrow=.N,ncol=1,byrow = T)),,
                      .SDcols = c(grep(paste0("^input",training_reps_num[j],"$"),names(training_data)),
                                  grep(paste0("^output",training_reps_num[j],"$"),names(training_data)))]
      }
      #merged fitness values
      training_data[,fitness_A := rowSums(.SD[,1:NTreps]/(.SD[,(NTreps+1):(2*NTreps)]^2),na.rm=T) / 
                      rowSums(1/(.SD[,(NTreps+1):(2*NTreps)]^2),na.rm=T),
                    ,.SDcols = c(grep(paste0("^fitness[", training_reps, "]$"),names(training_data)),
                                 grep(paste0("^error_A[", training_reps, "]$"),names(training_data)))]
      #and merged error
      training_data[,error_A := sqrt(1/rowSums(1/.SD^2)),,
                    .SDcols = grep(paste0("^error_A[", training_reps, "]$"),names(training_data))]
      
      #######################################
      ### variant-specific random effect error model (Rubin2017)
      training_data[input_above_threshold == T & all_reads ==T,random_error := dimsumms_errormodel_random_effect(training_data[input_above_threshold == T & all_reads ==T,.SD,,.SDcols = grep(paste0("^fitness[", training_reps, "]$"),names(training_data))],
                                                                                                                 training_data[input_above_threshold == T & all_reads ==T,.SD,,.SDcols = grep(paste0("^cbe[", training_reps, "]$"),names(training_data))])]
      training_data[input_above_threshold == T & all_reads ==T,fitness_ire := rowSums(.SD[,1:NTreps]/(.SD[,(NTreps+1):(2*NTreps)]^2 + random_error),na.rm=T) / 
                      rowSums(1/(.SD[,(NTreps+1):(2*NTreps)]^2 + random_error),na.rm=T),
                    ,.SDcols = c(grep(paste0("^fitness[", training_reps, "]"),names(training_data)),
                                 grep(paste0("^cbe[", training_reps, "]$"),names(training_data)))]
      training_data[input_above_threshold == T & all_reads ==T,error_ire := sqrt(1 / rowSums(1/(.SD[,1:NTreps]^2 + random_error),na.rm=T)),
                    ,.SDcols = grep(paste0("^cbe[", training_reps, "]$"),names(training_data))]
      
      
      ##############################
      ### predict test replicate ###
      ##############################
      
      #only use variants that have non-NA fitness and error values in all error models
      training_data[,test_rep_ok := !is.na(.SD) & is.finite(unlist(.SD)) & all_reads == T & input_above_threshold == T &
                      !is.na(fitness_naive) & is.finite(fitness_naive) & is.finite(error_naive) & error_naive != 0 &
                      !is.na(fitness_cbe) & is.finite(fitness_cbe) & is.finite(error_cbe) &
                      !is.na(fitness_ire) & is.finite(fitness_ire) & is.finite(error_ire) &
                      !is.na(fitness_Mio) & is.finite(fitness_Mio) & is.finite(error_Mio) &
                      !is.na(fitness_A) & is.finite(fitness_A) & is.finite(error_A) &
                      !is.na(fitness_MioA) & is.finite(fitness_MioA) & is.finite(error_MioA),,
                    .SDcols = paste0("fitness",test_rep)]
      
      #estimate model parameters for test rep as mean of those derived from training reps
      Mi = error_model_cycle[parameter == "input",mean(mean_value)]
      Mo = error_model_cycle[parameter == "output",mean(mean_value)]
      A = error_model_cycle[parameter == "reperror",mean(sqrt(mean_value))^2]
      #scaling factor from fitness normalization
      fitness_scale = fitness_norm_model[,unlist(.SD),.SDcols = paste0("scale_",test_rep)]
      
      #estimate ratio of fitness deviation in test vs training replicates versus estimated error
      report_moments = rbind(training_data[test_rep_ok==T,.(training_reps,test_rep,Nmut = "all",
                                                            # are error estimates of the right order of magnitude vs unseen replicate?
                                                            sd_naive = sd((fitness_naive-unlist(.SD[,1]))/(error_naive*sqrt(NTreps+1)),na.rm=T),
                                                            sd_cbe = sd((fitness_cbe-unlist(.SD[,1]))/sqrt(error_cbe^2 + unlist(.SD[,2])^2),na.rm=T),
                                                            sd_ire = sd((fitness_ire-unlist(.SD[,1]))/sqrt(error_ire^2 + unlist(.SD[,2])^2 + random_error),na.rm=T),
                                                            sd_A = sd((fitness_A-unlist(.SD[,1]))/sqrt(error_A^2 + unlist(.SD[,2])^2 + A),na.rm=T),
                                                            sd_Mio = sd((fitness_Mio-unlist(.SD[,1]))/sqrt(error_Mio^2 + fitness_scale*(Mi/unlist(.SD[,3]) + Mo/unlist(.SD[,4]))),na.rm=T),
                                                            sd_MioA = sd((fitness_MioA-unlist(.SD[,1]))/sqrt(error_MioA^2 + fitness_scale*(Mi/unlist(.SD[,3]) + Mo/unlist(.SD[,4])) + A),na.rm=T),
                                                            #robust estimates of sd
                                                            robust_naive = abs(diff(quantile((fitness_naive-unlist(.SD[,1]))/(error_naive*sqrt(NTreps+1)),probs = c(0.159,0.841),na.rm=T)))/2,
                                                            robust_cbe = abs(diff(quantile((fitness_cbe-unlist(.SD[,1]))/sqrt(error_cbe^2 + unlist(.SD[,2])^2),probs = c(0.159,0.841),na.rm=T)))/2, 
                                                            robust_ire = abs(diff(quantile((fitness_ire-unlist(.SD[,1]))/sqrt(error_ire^2 + unlist(.SD[,2])^2 + random_error),probs = c(0.159,0.841),na.rm=T)))/2, 
                                                            robust_A = abs(diff(quantile((fitness_A-unlist(.SD[,1]))/sqrt(error_A^2 + unlist(.SD[,2])^2 + A),probs = c(0.159,0.841),na.rm=T)))/2, 
                                                            robust_Mio = abs(diff(quantile((fitness_Mio-unlist(.SD[,1]))/sqrt(error_Mio^2 + fitness_scale*(Mi/unlist(.SD[,3]) + Mo/unlist(.SD[,4]))),probs = c(0.159,0.841),na.rm=T)))/2,
                                                            robust_MioA = abs(diff(quantile((fitness_MioA-unlist(.SD[,1]))/sqrt(error_MioA^2 + fitness_scale*(Mi/unlist(.SD[,3]) + Mo/unlist(.SD[,4])) + A),probs = c(0.159,0.841),na.rm=T)))/2,
                                                            
                                                            # are merged fitness values unbiased compared to unseen replicate?
                                                            mean_naive = (mean((fitness_naive-unlist(.SD[,1])),na.rm=T)), 
                                                            mean_cbe = (mean((fitness_cbe-unlist(.SD[,1])),na.rm=T)), 
                                                            mean_ire = (mean((fitness_ire-unlist(.SD[,1])),na.rm=T)), 
                                                            mean_A = (mean((fitness_A-unlist(.SD[,1])),na.rm=T)), 
                                                            mean_Mio = (mean((fitness_Mio-unlist(.SD[,1])),na.rm=T)), 
                                                            mean_MioA = (mean((fitness_MioA-unlist(.SD[,1])),na.rm=T))),
                                           ,.SDcols = c(paste0("fitness",test_rep),
                                                        paste0("cbe",test_rep),
                                                        paste0("input",test_rep),
                                                        paste0("output",test_rep))],
                             
                             ## same calculation but stratified by variants' number of mutations
                             training_data[test_rep_ok==T & Nmut > 0,.(training_reps,test_rep,
                                                                       # are error estimates of the right order of magnitude vs unseen replicate?
                                                                       sd_naive = sd((fitness_naive-unlist(.SD[,1]))/(error_naive*sqrt(NTreps+1)),na.rm=T), 
                                                                       sd_cbe = sd((fitness_cbe-unlist(.SD[,1]))/sqrt(error_cbe^2 + unlist(.SD[,2])^2),na.rm=T), 
                                                                       sd_ire = sd((fitness_ire-unlist(.SD[,1]))/sqrt(error_ire^2 + unlist(.SD[,2])^2 + random_error),na.rm=T), 
                                                                       sd_A = sd((fitness_A-unlist(.SD[,1]))/sqrt(error_A^2 + unlist(.SD[,2])^2 + A),na.rm=T), 
                                                                       sd_Mio = sd((fitness_Mio-unlist(.SD[,1]))/sqrt(error_Mio^2 + fitness_scale*(Mi/unlist(.SD[,3]) + Mo/unlist(.SD[,4]))),na.rm=T),
                                                                       sd_MioA = sd((fitness_MioA-unlist(.SD[,1]))/sqrt(error_MioA^2 + fitness_scale*(Mi/unlist(.SD[,3]) + Mo/unlist(.SD[,4])) + A),na.rm=T),
                                                                       #robust estimates of sd
                                                                       robust_naive = abs(diff(quantile((fitness_naive-unlist(.SD[,1]))/(error_naive*sqrt(NTreps+1)),probs = c(0.159,0.841),na.rm=T)))/2,
                                                                       robust_cbe = abs(diff(quantile((fitness_cbe-unlist(.SD[,1]))/sqrt(error_cbe^2 + unlist(.SD[,2])^2),probs = c(0.159,0.841),na.rm=T)))/2, 
                                                                       robust_ire = abs(diff(quantile((fitness_ire-unlist(.SD[,1]))/sqrt(error_ire^2 + unlist(.SD[,2])^2 + random_error),probs = c(0.159,0.841),na.rm=T)))/2, 
                                                                       robust_A = abs(diff(quantile((fitness_A-unlist(.SD[,1]))/sqrt(error_A^2 + unlist(.SD[,2])^2 + A),probs = c(0.159,0.841),na.rm=T)))/2, 
                                                                       robust_Mio = abs(diff(quantile((fitness_Mio-unlist(.SD[,1]))/sqrt(error_Mio^2 + fitness_scale*(Mi/unlist(.SD[,3]) + Mo/unlist(.SD[,4]))),probs = c(0.159,0.841),na.rm=T)))/2,
                                                                       robust_MioA = abs(diff(quantile((fitness_MioA-unlist(.SD[,1]))/sqrt(error_MioA^2 + fitness_scale*(Mi/unlist(.SD[,3]) + Mo/unlist(.SD[,4])) + A),probs = c(0.159,0.841),na.rm=T)))/2,
                                                                       # are merged fitness values unbiased compared to unseen replicate?
                                                                       mean_naive = (mean((fitness_naive-unlist(.SD[,1])),na.rm=T)), 
                                                                       mean_cbe = (mean((fitness_cbe-unlist(.SD[,1])),na.rm=T)), 
                                                                       mean_ire = (mean((fitness_ire-unlist(.SD[,1])),na.rm=T)), 
                                                                       mean_A = (mean((fitness_A-unlist(.SD[,1])),na.rm=T)), 
                                                                       mean_Mio = (mean((fitness_Mio-unlist(.SD[,1])),na.rm=T)), 
                                                                       mean_MioA = (mean((fitness_MioA-unlist(.SD[,1])),na.rm=T))),
                                           .(Nmut),.SDcols = c(paste0("fitness",test_rep),
                                                               paste0("cbe",test_rep),
                                                               paste0("input",test_rep),
                                                               paste0("output",test_rep))])
      
      
      if (i == 1) {
        report_moments_all = report_moments
      } else {
        report_moments_all = rbind(report_moments_all,report_moments)
      }
    }
    
    #save results
    X = melt(report_moments_all,id.vars = c("Nmut","training_reps","test_rep"))
    setkey(X,Nmut,test_rep,variable)
    Y = cbind(X[grepl("sd",variable),.(Nmut,test_rep,sd=variable,sd_val = value)],
              X[grepl("robust",variable),.(wq=variable,wq_val = value)],
              X[grepl("^mean",variable),.(mean = variable,mean_val = value)])
    Y[,method := paste0(strsplit(x=as.character(sd),"_")[[1]][-1],collapse="_"),sd]
    write.table(x = rbind(Y[,.(Nmut,test_rep = as.character(test_rep),method,
                               sd_val = round(sd_val,digits=2),
                               wq_val = round(wq_val,digits=2),
                               mean_val=round(mean_val,digits=2))],
                          Y[,.(test_rep = "avg",
                               sd_val = round(mean(sd_val),digits=2),
                               wq_val = round(wq_val,digits=2),
                               mean_val=round(mean(mean_val),digits=2)),.(Nmut,method)]),
                file = gsub("//","/",file.path(outpath,paste0(dataset_label[dataset_idx],"_leaveoneout_performance.txt"))),quote=F,row.names=F)
  }
  
  #####################################
  ### Error model performance plots ###
  #####################################
  
  #Load data
  files = list.files(outpath)
  dataset_label = gsub("_leaveoneout_performance.txt","",files)
  
  #look at errors across datasets
  for (i in seq_along(files)) {
    X <- fread(file.path(outpath , paste0(dataset_label[i],"_leaveoneout_performance.txt")))
    X[,dataset := dataset_label[i]]
    
    if (i == 1) {
      all_datasets <- X
    } else {
      all_datasets <- rbind(all_datasets,X)
    }
  }
  performance <- all_datasets[test_rep != "avg",.(value = 2^mean((log2(sd_val)))),.(Nmut,method,test_rep,dataset)]
  
  dataset_dict <- list(
    "TDP43_290" = "TDP-43 (290-331); Bolognesi et al. 2019",
    "TDP43_332" = "TDP-43 (332-373); Bolognesi et al. 2019",
    "FOSJUN" = "FOS-JUN deepPCA; Diss et al. 2018",
    "FOScis" = "FOS deepPCA; Diss et al. 2018",
    "GRB2_GPD" = "GRB2 deepPCA; in prep.",
    "GRB2_CYC" = "GRB2 stabilityPCA; in prep.",
    "GB1" = "GB1; Olson et al. 2014",
    "tRNA_sel23" = "tRNA 23C; Li & Zhang 2018",
    "tRNA_sel30" = "tRNA 30C; Li & Zhang 2018",
    "tRNA_sel37" = "tRNA 37C; Li & Zhang 2018",
    "tRNA_selDMSO" = "tRNA 30C, 3% DMSO; Li & Zhang 2018")
  performance[,dataset := unlist(dataset_dict[dataset])]
  performance[,method := factor(method,levels=c("naive","cbe","ire","MioA","A","Mio"))]
  levels(performance$method) <- c("s.d.-based","count-based","CB + variant corr.","DiMSum full","DiMSum rep.","DiMSum mult.")
  performance[,Nmut := factor(Nmut,levels=c("1","2","all"))]
  levels(performance$Nmut) <- c("single mutants","double mutants","all variants")
  p <- ggplot2::ggplot(performance[!method %in% c("DiMSum rep.", "DiMSum mult.") & Nmut=="all variants",],ggplot2::aes(x = method, y = 1/(value))) +
    ggplot2::geom_boxplot(ggplot2::aes(group = method), outlier.shape = NA) +
    ggplot2::geom_point(ggplot2::aes(color = dataset), position = ggplot2::position_dodge(width = 0.6), size = 2) +
    ggplot2::geom_hline(yintercept = 1,lty = 2) +
    ggplot2::labs(y = "error model performance", x = "") +
    ggplot2::scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1,1.25,1.5),expand=c(0,0.015)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
  ggplot2::ggsave(plot = p, filename = file.path(outpath, "leaveoneout_crossvalidation_fig3.pdf"),width=6,height=5)
  
  p <- ggplot2::ggplot(performance,ggplot2::aes(x = method, y = 1/(value))) +
    ggplot2::geom_boxplot(ggplot2::aes(group = method), outlier.shape = NA) +
    ggplot2::geom_point(ggplot2::aes(color = dataset), position = ggplot2::position_dodge(width = 0.6), size = 2) +
    ggplot2::facet_wrap(~Nmut) +
    ggplot2::geom_hline(yintercept = 1,lty = 2) +
    ggplot2::geom_vline(xintercept = 4.5,lty = 3) +
    ggplot2::labs(y = "error model performance", x = "") +
    ggplot2::scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1,1.5,2,2.5,3),expand=c(0,0.015)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
  ggplot2::ggsave(plot = p, filename = file.path(outpath, "leaveoneout_crossvalidation_suppFigure10.pdf"),width=12,height=9)
  
  
}