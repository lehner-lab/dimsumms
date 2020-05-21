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
  
    #Preprocess datasets
    dimsumms_errormodel_leaveoneout_preprocess_datasets(dataset_dir)

    #fetch files
    files = list.files(file.path(dataset_dir, "processed_data"))
    dataset_label = sapply(1:length(files),function(X){strsplit(files[X],"\\.")[[1]][1]})
    
    ###############################
    ###  calculate error model ####
    ###############################
    
    for (dataset_idx in seq_along(files)) {
      
      print(paste0("running leave one out cross validation for ",dataset_label[dataset_idx]))
      
      #load dataset
      work_data = fread(file.path(dataset_dir, "processed_data", files[dataset_idx]))
      
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
      ### leave-one-out cross validation ###
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
        training_data[,error_naive := apply(.SD, 1, sd) / sqrt(nchar(training_reps)),
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
        # table with error model parameters
        error_model = data.table(parameter = rep(c("input", "output", "reperror"), each = NTreps),
                                       rep = c(rep(training_reps_num,3)),
                                       mean_value = colMeans(parameters, na.rm = T),
                                       sd_value = apply(parameters, 2, sd, na.rm = T),
                                       ensemble = sum(!is.na(parameters[, 1])),
                                       which_reps = training_reps)

        if (i == 1) {
           parameters_full = dimsumms_errormodel_fit(DT = training_data,
                                                       reps = reps,
                                                       Ncores = Ncores,
                                                       Fcorr = Fcorr)
          error_model_full = data.table(parameter = c(rep(c("input","output", "reperror"), each = nchar(reps))),
                                       rep = c(rep(reps_num,3)),
                                       mean_value = round(colMeans(parameters_full, na.rm = T), 4),
                                       sd_value = round(apply(parameters_full, 2, sd, na.rm = T), 4))

          error_model_full = rbind(error_model_full, 
            error_model_full[parameter %in% c("input", "output"), .(rep = "average",
              mean_value = round(mean(mean_value),3), 
              sd_value = round(sd(mean_value),3)), parameter],
            error_model_full[parameter %in% c("reperror"), .(parameter = "sqrt-reperror", 
              rep = "average",
              mean_value = round(mean(sqrt(mean_value)),3), 
              sd_value = round(sd(sqrt(mean_value)),3))])

          write.table(error_model_full,
            file = gsub("//","/",file.path(outpath,paste0(dataset_label[dataset_idx],"_full_errormodel.txt"))),
            quote = F,
            row.names = F)
        }
        

        #######
        #fitness and estimated error from full error model (MioA - Multiplicative Input/Output terms + replicate-specific additive error terms)
        for (j in seq_along(training_reps_num)) {
          Corr = matrix(unlist(fitness_norm_model[,.SD,.SDcols=paste0("scale_",training_reps_num[j])]),ncol = 1,nrow=nrow(training_data))
          training_data[,paste0("error_MioA",training_reps_num[j]) := sqrt(Corr * rowSums(matrix(unlist(error_model[parameter %in% c("input","output") & rep == training_reps_num[j],mean_value]),nrow=.N,ncol=2,byrow = T)/.SD) + 
                                                                             matrix(error_model[parameter %in% c("reperror") & rep == training_reps_num[j],mean_value],nrow = .N, ncol = 1, byrow = T)),,
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
        

        #######
        # from fitting using all replicates and replicate-specific additive terms
        for (j in seq_along(training_reps_num)) {
          Corr = matrix(unlist(fitness_norm_model[,.SD,.SDcols=paste0("scale_",training_reps_num[j])]),ncol = 1,nrow=nrow(training_data))
          training_data[,paste0("error_MioA_full",training_reps_num[j]) := sqrt(Corr * rowSums(matrix(unlist(error_model_full[parameter %in% c("input","output") & rep == training_reps_num[j],mean_value]),nrow=.N,ncol=2,byrow = T)/.SD) + 
                                                                             matrix(error_model_full[parameter %in% c("reperror") & rep == training_reps_num[j],mean_value],nrow = .N, ncol = 1, byrow = T)),,
                        .SDcols = c(grep(paste0("^input",training_reps_num[j],"$"),names(training_data)),
                                    grep(paste0("^output",training_reps_num[j],"$"),names(training_data)))]
        }

        #merged fitness values
        training_data[,fitness_MioA_full := rowSums(.SD[,1:NTreps]/(.SD[,(NTreps+1):(2*NTreps)]^2),na.rm=T) / 
                        rowSums(1/(.SD[,(NTreps+1):(2*NTreps)]^2),na.rm=T),
                      ,.SDcols = c(grep(paste0("^fitness[", training_reps, "]$"),names(training_data)),
                                   grep(paste0("^error_MioA_full[", training_reps, "]$"),names(training_data)))]
        #and merged error
        training_data[,error_MioA_full := sqrt(1/rowSums(1/.SD^2)),,
                      .SDcols = grep(paste0("^error_MioA_full[", training_reps, "]$"),names(training_data))]



        ########
        #fitness and estimated error when ONLY using multiplicative input/output terms 
        for (j in seq_along(training_reps_num)) {
          Corr = matrix(unlist(fitness_norm_model[,.SD,.SDcols=paste0("scale_",training_reps_num[j])]),ncol = 1,nrow=nrow(training_data))
          
          training_data[,paste0("error_Mio",training_reps_num[j]) := sqrt(Corr * rowSums(matrix(unlist(error_model[parameter %in% c("input","output") & rep == training_reps_num[j],mean_value]),nrow=.N,ncol=2,byrow = T)/.SD)),,
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
                                                                          matrix(error_model[parameter %in% c("reperror") & rep == training_reps_num[j],mean_value],nrow=.N,ncol=1,byrow = T)),,
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
        training_data[input_above_threshold == T & all_reads ==T,
          random_error := dimsumms_errormodel_random_effect(training_data[input_above_threshold == T & all_reads ==T,.SD,,.SDcols = grep(paste0("^fitness[", training_reps, "]$"),names(training_data))],
                                                            training_data[input_above_threshold == T & all_reads ==T,.SD,,.SDcols = grep(paste0("^cbe[", training_reps, "]$"),names(training_data))])]

        training_data[input_above_threshold == T & all_reads ==T,fitness_ire := rowSums(.SD[,1:NTreps]/(.SD[,(NTreps+1):(2*NTreps)]^2 + random_error),na.rm=T) / 
                        rowSums(1/(.SD[,(NTreps+1):(2*NTreps)]^2 + random_error),na.rm=T),
                      ,.SDcols = c(grep(paste0("^fitness[", training_reps, "]"),names(training_data)),
                                   grep(paste0("^cbe[", training_reps, "]$"),names(training_data)))]
        training_data[input_above_threshold == T & all_reads ==T,error_ire := sqrt(1 / rowSums(1/(.SD[,1:NTreps]^2 + random_error),na.rm=T)),
                      ,.SDcols = grep(paste0("^cbe[", training_reps, "]$"),names(training_data))]

        #######################################
        ### bayesian estimate of variance (Weile2019)
        training_data[, input_min_count := apply(.SD, 1, min), 
          .SDcols = grep(paste0("input[", training_reps, "]"), names(training_data))]
        splinemat <- training_data[input_above_threshold == T & all_reads == T, .(logmin = log10(input_min_count),
                                                                    fitness = fitness_naive,
                                                                    logsd = log10(error_naive))]

        z <- lm(logsd ~ fitness + logmin, data = splinemat)
        training_data[, error_naive_fit := 10^predict(z, newdata = .(logmin = log10(input_min_count),fitness = fitness_naive))]
        training_data[, error_br := sqrt((3 * error_naive_fit^2 + (nchar(training_reps) - 1) * error_naive^2) / (1 + nchar(training_reps)))]
        
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
        Mi = error_model[parameter == "input",mean(mean_value)]
        Mo = error_model[parameter == "output",mean(mean_value)]
        A = error_model[parameter == "reperror", mean(sqrt(mean_value))^2]

        Mi_full = error_model_full[parameter == "input" & rep == test_rep,mean(mean_value)]
        Mo_full = error_model_full[parameter == "output" & rep == test_rep,mean(mean_value)]
        A_full = error_model_full[parameter == "reperror" & rep == test_rep, mean_value]

        #scaling factor from fitness normalization
        fitness_scale = fitness_norm_model[,unlist(.SD),.SDcols = paste0("scale_",test_rep)]
        
        dt_zscore <- training_data[test_rep_ok==T,
            .(test_rep, Nmut,
                naive = (fitness_naive-unlist(.SD[,1]))/(error_naive*sqrt(NTreps+1)),
                br = (fitness_naive-unlist(.SD[,1]))/(error_br*sqrt(NTreps+1)),
                cbe = (fitness_cbe-unlist(.SD[,1]))/sqrt(error_cbe^2 + unlist(.SD[,2])^2),
                ire = (fitness_ire-unlist(.SD[,1]))/sqrt(error_ire^2 + unlist(.SD[,2])^2 + random_error),
                MioA = (fitness_MioA-unlist(.SD[,1]))/sqrt(error_MioA^2 + fitness_scale*(Mi/unlist(.SD[,3]) + Mo/unlist(.SD[,4])) + A),
                A = (fitness_A-unlist(.SD[,1]))/sqrt(error_A^2 + A),
                Mio = (fitness_Mio-unlist(.SD[,1]))/sqrt(error_Mio^2 + fitness_scale*(Mi/unlist(.SD[,3]) + Mo/unlist(.SD[,4]))),
                MioA_full = (fitness_MioA_full-unlist(.SD[,1]))/sqrt(error_MioA_full^2 + fitness_scale*(Mi_full/unlist(.SD[,3]) + Mo_full/unlist(.SD[,4])) + A_full)),
           ,.SDcols = c(paste0("fitness",test_rep),
                        paste0("cbe",test_rep),
                        paste0("input",test_rep),
                        paste0("output",test_rep))]

          if (i == 1) {
            dt_zscore_all = dt_zscore
          } else {
            dt_zscore_all = rbind(dt_zscore_all, dt_zscore)
          }
        }

        # save results for dataset
        write.table(x = dt_zscore_all,
          file = gsub("//","/",file.path(outpath,paste0(dataset_label[dataset_idx],"_leaveoneout_zscore.txt"))),quote=F,row.names=F)

        # plot results for dataset
        dt_zscore_melt <- melt(dt_zscore_all[sample(.N, min(1e4, .N))], id.vars = c("test_rep", "Nmut"))
        dt_zscore_melt[, type := "leave-one-out"]
        dt_zscore_melt[variable == "MioA_full", type := "all replicates"]
        dt_zscore_melt[, type := factor(type, levels = c("leave-one-out", "all replicates"))]

        method_dict <- list(
          "naive" = "s.d.-based", 
          "br" = "bayes-reg s.d., Weile et al. 2019", 
          "cbe" = "count-based", 
          "ire" = "Enrich2, Rubin et al. 2017", 
          "MioA" = "DiMSum",
          "MioA_full" = "DiMSum")
        dt_zscore_melt[, variable := factor(variable, levels = names(method_dict))]
        levels(dt_zscore_melt$variable) <- unlist(method_dict)

        zrange = 3.5
        p1 <- ggplot2::ggplot(dt_zscore_melt, ggplot2::aes(sample = value, color = variable, lty = type)) +
          ggplot2::geom_abline(lty = 2) + 
          ggplot2::geom_qq(geom = "line") +
          ggplot2::coord_cartesian(xlim = c(-zrange, zrange), ylim = c(-zrange, zrange)) +
          ggplot2::labs(x = 'expected', y = 'observed', color = "error model", lty = "") +
          ggplot2::theme_bw()

        dt_zscore_melt[, pvalue := 2*pnorm(-abs(value))]

        p2 <- ggplot2::ggplot(dt_zscore_melt, ggplot2::aes(pvalue, color = variable, lty = type)) +
          ggplot2::geom_abline(lty = 2) + 
          ggplot2::stat_ecdf() +
          ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0.01)) +
          ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0.01)) +
          ggplot2::labs(x = 'p value', y = 'empirical CDF', color = "error model", lty = "") +
          ggplot2::theme_bw()
        p3 <- ggplot2::ggplot(dt_zscore_melt, ggplot2::aes(pvalue, color = variable, lty = type)) +
          ggplot2::geom_abline(lty = 2) + 
          ggplot2::stat_ecdf() +
          ggplot2::coord_cartesian(xlim = c(0, .1), ylim = c(0, .1)) +
          ggplot2::scale_x_continuous(expand = c(0, 0.01)) +
          ggplot2::scale_y_continuous(expand = c(0, 0.01)) +
          ggplot2::labs(x = 'p value', y = 'empirical CDF', color = "error model", lty = "") +
          ggplot2::theme_bw()

        p <- gridExtra::grid.arrange(p1, p2, p3,  nrow = 2)

        ggplot2::ggsave(plot = p, file = gsub("//","/",file.path(outpath,paste0(dataset_label[dataset_idx],"_leaveoneout_pvalue.pdf"))),
              width = 12, height = 6)

      }
    
    #####################################
    ### Error model performance plots ###
    #####################################
    

    ####################
    #Load data
    files = list.files(outpath)
    files = grep("_leaveoneout_zscore.txt", files, value = T)
    dataset_label = gsub("_leaveoneout_zscore.txt","",files)
    
    #look at errors across datasets
    for (i in seq_along(files)) {
      dt_zscore <- fread(file.path(outpath , paste0(dataset_label[i],"_leaveoneout_zscore.txt")))
      dt_zscore_melt <- melt(dt_zscore, id.vars = c("test_rep", "Nmut"))
      dt_zscore_melt[, type := "leave-one-out"]
      dt_zscore_melt[variable == "MioA_full", type := "all replicates"]
      dt_zscore_melt[, type := factor(type, levels = c("leave-one-out", "all replicates"))]

      method_dict <- list(
        "naive" = "s.d.-based", 
        "br" = "bayes-reg s.d., Weile et al. 2019", 
        "cbe" = "count-based", 
        "ire" = "Enrich2, Rubin et al. 2017", 
        "MioA" = "DiMSum",
        "A" = "DiMSum additive only",
        "Mio" = "DiMSum multiplicative only",
        "MioA_full" = "DiMSum")
      dt_zscore_melt[, variable := factor(variable, levels = names(method_dict))]
      levels(dt_zscore_melt$variable) <- unlist(method_dict)

      error_model <- fread(file.path(outpath , paste0(dataset_label[i],"_full_errormodel.txt")))

      zscore_moments <- rbind(dt_zscore_melt[,.(dataset = dataset_label[i], 
          Nmut = "all",
          mult = mean(error_model[rep == "average" & parameter %in% c("input", "output"), mean_value]),
          sd_value = sd(value),
          q_value = diff(quantile(value, probs = c(0.159, 0.841)))/2), 
        .(variable, type)],
        dt_zscore_melt[between(Nmut, 1, 2),.(dataset = dataset_label[i], 
          mult = mean(error_model[rep == "average" & parameter %in% c("input", "output"), mean_value]),
          sd_value = sd(value),
          q_value = diff(quantile(value, probs = c(0.159, 0.841)))/2), 
        .(variable, type, Nmut)])

      if (i == 1) {
        zscore_moments_all <- zscore_moments
      } else {
        zscore_moments_all <- rbind(zscore_moments_all,zscore_moments)
      }
    }
    
    Nreps_dict <- list(
      "TDP43_290" = 3,
      "TDP43_332" = 4,
      "FOSJUN" = 3,
      "FOScis" = 3,
      "GRB2_GPD" = 3,
      "GRB2_CYC" = 3,
      "GB1" = 3,
      "tRNA_Phylogeny" = 6,
      "tRNA_sel23" = 5,
      "tRNA_sel30" = 5,
      "tRNA_sel37" = 3,
      "tRNA_selDMSO" = 3)
    zscore_moments_all[, Nreps := unlist(Nreps_dict[dataset])]

    dataset_dict <- list(
      "TDP43_290" = "TDP-43 (290-331); Bolognesi et al. 2019",
      "TDP43_332" = "TDP-43 (332-373); Bolognesi et al. 2019",
      "FOSJUN" = "FOS-JUN; Diss et al. 2018",
      "FOScis" = "FOS; Diss et al. 2018",
      "GRB2_GPD" = "GRB2 unpublished dataset 1; Domingo et al. in prep.",
      "GRB2_CYC" = "GRB2 unpublished dataset 2; Domingo et al. in prep.",
      "GB1" = "GB1; Olson et al. 2014",
      "tRNA_Phylogeny" = "tRNA 37C, NaCl; Domingo et al. 2018",
      "tRNA_sel23" = "tRNA 23C; Li & Zhang 2018",
      "tRNA_sel30" = "tRNA 30C; Li & Zhang 2018",
      "tRNA_sel37" = "tRNA 37C; Li & Zhang 2018",
      "tRNA_selDMSO" = "tRNA 30C, 3% DMSO; Li & Zhang 2018")
    zscore_moments_all[, dataset := factor(dataset, levels = names(dataset_dict))]
    levels(zscore_moments_all$dataset) <- unlist(dataset_dict)

    

    zscore_moments_all[,Nmut := factor(Nmut, levels=c("1", "2", "all"))]
    levels(zscore_moments_all$Nmut) <- c("single mutants", "double mutants", "all variants")

    p <- ggplot2::ggplot(zscore_moments_all[Nmut=="all variants" & !is.na(dataset) & type == "leave-one-out"],
      ggplot2::aes(x = variable, y = 1/sd_value)) +
      ggplot2::geom_boxplot(ggplot2::aes(group = variable), outlier.shape = NA) +
      ggplot2::geom_point(ggplot2::aes(color = dataset), position = ggplot2::position_dodge(width = 0.6), size = 2) +
      ggplot2::geom_hline(yintercept = 1,lty = 2) +
      ggplot2::geom_vline(xintercept = 5.5,lty = 3) +
      ggplot2::labs(y = "error model performance", x = "") +
      ggplot2::scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1,1.25,1.5),expand=c(0,0.015)) +
      # ggplot2::scale_y_log10(expand=c(0,0.015)) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
    ggplot2::ggsave(plot = p, 
      filename = file.path(outpath, "leaveoneout_crossvalidation_zscore_mainFigure.pdf"),
      width=6,
      height=7)
      
    p <- ggplot2::ggplot(zscore_moments_all[!is.na(dataset)],
         ggplot2::aes(x = variable, y = sd_value)) +
      ggplot2::geom_boxplot(ggplot2::aes(group = variable), outlier.shape = NA) +
      ggplot2::geom_point(ggplot2::aes(color = log10(mult), shape = factor(Nreps)), 
        position = ggplot2::position_dodge(width = .6), size = 2) +
      ggplot2::facet_wrap(~Nmut) +
      ggplot2::geom_hline(yintercept = 1,lty = 2) +
      ggplot2::geom_vline(xintercept = 5.5,lty = 3) +
      ggplot2::labs(y = "error model performance", x = "", color = "log10(over-dispersion)", shape = "# replicates") +
      ggplot2::scale_color_gradient(low = "#fbc98e", high = "#c8b0d3") +
      ggplot2::scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1,1.5,2,2.5,3),expand=c(0,0.015)) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
    ggplot2::ggsave(plot = p, filename = file.path(outpath, "leaveoneout_crossvalidation_zscore_supplementaryFigure.pdf"),width=12,height=9)

  }