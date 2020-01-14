#' dimsumms_errormodel_leaveoneout_preprocess_datasets
#'
#' Preprocess datasets for leave one out cross validation benchmarking
#'
#' @param dataset_dir base directory with downloaded DiMSum manuscript data (folder 'misc')
#' @param outpath output path for processed datasets
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
dimsumms_errormodel_leaveoneout_preprocess_datasets <- function(
  dataset_dir,
  execute = TRUE
){
  #Return if analysis not executed
  if(!execute){
    return()
  }

  # this script preprocesses (cuts down) DiMSum processed datasets to feed into the leave one out cross validation script
  # the DiMSum datasets themselves are not provided, but the shell scripts need to run DMS to recreate them from the raw data
  
############ TDP43
#290-331 library
all_data = fread(file.path(dataset_dir,"misc/DiMSum_errormodel/datasets/TDP43_290_2017-06-13_variant_data_merge.tsv"))[Nham_nt < 3 & (input1_e1_s0_bNA_count + input3_e3_s0_bNA_count+ input4_e4_s0_bNA_count) > 10,
                                   .(nt_seq,aa_seq,WT,STOP,STOP_readthrough,Nham_nt,Nham_aa,barcode_valid,indel,constant_region,permitted,too_many_substitutions,
                                     Nmut = Nham_nt,
                                     input1 = input1_e1_s0_bNA_count,
                                     input3 = input3_e3_s0_bNA_count,
                                     input4 = input4_e4_s0_bNA_count,
                                     output1 = output1A_e1_s1_b1_count + output1B_e1_s1_b2_count,
                                     output3 = output3A_e3_s1_b1_count + output3B_e3_s1_b2_count,
                                     output4 = output4A_e4_s1_b1_count + output4B_e4_s1_b2_count)]
write.table(all_data,file.path(dataset_dir,"misc/DiMSum_errormodel/processed_data/TDP43_290.txt"),row.names=F,quote=F)

#332-373 library
all_data = fread(file.path(dataset_dir,'misc/DiMSum_errormodel/datasets/TDP43_332_2018-04-06_variant_data_merge.txt'))[Nmut_nt <= 2 & Nins_nt == 0 & Ndel_nt == 0,
                                   .(Nmut = Nsub_nt,STOP,WT,
                                     input1 = input1_e1_s0_bNA_count,
                                     input2 = input2_e2_s0_bNA_count,
                                     input3 = input3_e3_s0_bNA_count,
                                     input4 = input4_e4_s0_bNA_count,
                                     output1 = output1B_e1_s1_b2_count,
                                     output2 = output2A_e2_s1_b1_count + output2B_e2_s1_b2_count,
                                     output3 = output3A_e3_s1_b1_count + output3B_e3_s1_b2_count,
                                     output4 = output4A_e4_s1_b1_count + output4B_e4_s1_b2_count)]
all_data[is.na(WT),WT := FALSE]
write.table(all_data,file.path(dataset_dir,"misc/DiMSum_errormodel/processed_data/TDP43_332.txt"),row.names=F,quote=F)


###tRNA Li&Zhang 2018
#23C
all_data = fread(file.path(dataset_dir,"misc/DiMSum_errormodel/datasets/Li2018_tRNA_sel23_variant_data_merge.tsv"))[Nham_nt < 3 & inputA_e1_s0_bNA_count > 2000 & sel23A_e1_s1_b1_count > 200 & sel23B_e2_s1_b1_count > 200 & sel23C_e3_s1_b1_count > 200 & sel23D_e4_s1_b1_count > 200 & sel23E_e5_s1_b1_count > 200,
                                                                                             .(Nmut = Nham_nt,nt_seq,
                                                                                               input1 = inputA_e1_s0_bNA_count,
                                                                                               input2 = inputA_e1_s0_bNA_count,
                                                                                               input3 = inputA_e1_s0_bNA_count,
                                                                                               input4 = inputA_e1_s0_bNA_count,
                                                                                               input5 = inputA_e1_s0_bNA_count,
                                                                                               output1 = sel23A_e1_s1_b1_count,
                                                                                               output2 = sel23B_e2_s1_b1_count,
                                                                                               output3 = sel23C_e3_s1_b1_count,
                                                                                               output4 = sel23D_e4_s1_b1_count,
                                                                                               output5 = sel23E_e5_s1_b1_count)]
all_data[,WT := nt_seq == "ttccgttggcgtaatggtaacgcgtctccctcctaaggagaagactgcgggttcgagtcccgtacggaa"]
write.table(all_data,file.path(dataset_dir,"misc/DiMSum_errormodel/processed_data/tRNA_sel23.txt"),row.names=F,quote=F)

#30C
all_data = fread(file.path(dataset_dir,"misc/DiMSum_errormodel/datasets/Li2018_tRNA_sel30_variant_data_merge.tsv"))[Nham_nt < 3 & inputA_e1_s0_bNA_count > 2000 & sel30A_e1_s1_b1_count > 200 & sel30B_e2_s1_b1_count > 200 & sel30C_e3_s1_b1_count > 200 & sel30D_e4_s1_b1_count > 200 & sel30E_e5_s1_b1_count > 200,
                                                                                              .(Nmut = Nham_nt,nt_seq,
                                                                                                input1 = inputA_e1_s0_bNA_count,
                                                                                                input2 = inputA_e1_s0_bNA_count,
                                                                                                input3 = inputA_e1_s0_bNA_count,
                                                                                                input4 = inputA_e1_s0_bNA_count,
                                                                                                input5 = inputA_e1_s0_bNA_count,
                                                                                                output1 = sel30A_e1_s1_b1_count,
                                                                                                output2 = sel30B_e2_s1_b1_count,
                                                                                                output3 = sel30C_e3_s1_b1_count,
                                                                                                output4 = sel30D_e4_s1_b1_count,
                                                                                                output5 = sel30E_e5_s1_b1_count)]
all_data[,WT := nt_seq == "ttccgttggcgtaatggtaacgcgtctccctcctaaggagaagactgcgggttcgagtcccgtacggaa"]
write.table(all_data,file.path(dataset_dir,"misc/DiMSum_errormodel/processed_data/tRNA_sel30.txt"),row.names=F,quote=F)

#37C
all_data = fread(file.path(dataset_dir,"misc/DiMSum_errormodel/datasets/Li2018_tRNA_sel37_variant_data_merge.tsv"))[Nham_nt < 3 & inputA_e1_s0_bNA_count > 2000 & sel37A_e1_s1_b1_count > 200 & sel37B_e2_s1_b1_count > 200 & sel37C_e3_s1_b1_count > 200,
                                                                                             .(Nmut = Nham_nt,nt_seq,
                                                                                               input1 = inputA_e1_s0_bNA_count,
                                                                                               input2 = inputA_e1_s0_bNA_count,
                                                                                               input3 = inputA_e1_s0_bNA_count,
                                                                                               output1 = sel37A_e1_s1_b1_count,
                                                                                               output2 = sel37B_e2_s1_b1_count,
                                                                                               output3 = sel37C_e3_s1_b1_count)]
all_data[,WT := nt_seq == "ttccgttggcgtaatggtaacgcgtctccctcctaaggagaagactgcgggttcgagtcccgtacggaa"]
write.table(all_data,file.path(dataset_dir,"misc/DiMSum_errormodel/processed_data/tRNA_sel37.txt"),row.names=F,quote=F)

#30C, 3%DMSO
all_data = fread(file.path(dataset_dir,"misc/DiMSum_errormodel/datasets/Li2018_tRNA_selDMSO_variant_data_merge.tsv"))[Nham_nt < 3 & inputA_e1_s0_bNA_count > 2000 & selDMSOA_e1_s1_b1_count > 200 & selDMSOB_e2_s1_b1_count > 200 & selDMSOC_e3_s1_b1_count > 200,
                                                                                               .(Nmut = Nham_nt,nt_seq,
                                                                                                 input1 = inputA_e1_s0_bNA_count,
                                                                                                 input2 = inputA_e1_s0_bNA_count,
                                                                                                 input3 = inputA_e1_s0_bNA_count,
                                                                                                 output1 = selDMSOA_e1_s1_b1_count,
                                                                                                 output2 = selDMSOB_e2_s1_b1_count,
                                                                                                 output3 = selDMSOC_e3_s1_b1_count)]
all_data[,WT := nt_seq == "ttccgttggcgtaatggtaacgcgtctccctcctaaggagaagactgcgggttcgagtcccgtacggaa"]
write.table(all_data,file.path(dataset_dir,"misc/DiMSum_errormodel/processed_data/tRNA_selDMSO.txt"),row.names=F,quote=F)


########### GRB2 Domingo et al. in preparation
# GPD stabilityPCA assay
all_data = fread(file.path(dataset_dir,"misc/DiMSum_errormodel/datasets/GRB2_GPD_variant_data_merge_Q25.txt"))[Nmut_nt <= 2 & Nins_nt == 0 & Ndel_nt == 0,
                                                                                               .(Nmut = Nsub_nt,STOP,WT,
                                                                                                 input1 = input1_e1_s0_bNA_count,
                                                                                                 input2 = input2_e2_s0_bNA_count,
                                                                                                 input3 = input3_e3_s0_bNA_count,
                                                                                                 output1 = output1_e1_s1_b1_count,
                                                                                                 output2 = output2_e2_s1_b1_count,
                                                                                                 output3 = output3_e3_s1_b1_count)]
all_data[is.na(WT),WT:=F]
write.table(all_data,file.path(dataset_dir,"misc/DiMSum_errormodel/processed_data/GRB2_GPD.txt"),row.names=F,quote=F)

# CYC bindingPCA assay
all_data = fread(file.path(dataset_dir,"misc/DiMSum_errormodel/datasets/GRB2_CYC_variant_data_merge_Q25.txt"))[Nmut_nt <= 2 & Nins_nt == 0 & Ndel_nt == 0,
                                                                                               .(Nmut = Nsub_nt,STOP,WT,
                                                                                                 input1 = input1_e1_s0_bNA_count,
                                                                                                 input2 = input2_e2_s0_bNA_count,
                                                                                                 input3 = input3_e3_s0_bNA_count,
                                                                                                 output1 = output1_e1_s1_b1_count,
                                                                                                 output2 = output2_e2_s1_b1_count,
                                                                                                 output3 = output3_e3_s1_b1_count)]
all_data[is.na(WT),WT:=F]
write.table(all_data,file.path(dataset_dir,"misc/DiMSum_errormodel/processed_data/GRB2_CYC.txt"),row.names=F,quote=F)

############ FOS-JUN Diss&Lehner 2018
# FOS-JUN trans deepPCA (mutations in both proteins)
output = fread(file.path(dataset_dir,"misc/DiMSum_errormodel/datasets/FOSJUN_trans_Q20.txt"))
output[wt_aa_F == mut_aa_F,':=' (aa_pos_F = NA,mut_aa_F = NA)]
output[wt_aa_J == mut_aa_J,':=' (mut_aa_J = NA,aa_pos_J = NA)]

output = output[,.(input1 = sum(i1),input2 = sum(i2),input3 = sum(i3),
                   output1 = sum(o1),output2 = sum(o2),output3 = sum(o3)),
                .(aa_pos_F,aa_pos_J,mut_aa_F,mut_aa_J)]
names(output)[1:4] = c("Pos1","Pos2","Mut1","Mut2")
output[,Pos1 := Pos1 + 1] #correct positions starting from 0
output[,Pos2 := Pos2 + 1]
output[,FOSid := paste0(Pos1,Mut1),.(Pos1,Mut1)]
output[,JUNid := paste0(Pos2,Mut2),.(Pos2,Mut2)]
output[,Nmut := as.numeric(!is.na(Mut1)) + as.numeric(!is.na(Mut2)),.(FOSid,JUNid)]
output[,WT := Nmut == 0]
write.table(output[,.(Nmut,WT,input1,input2,input3,output1,output2,output3)],
            file.path(dataset_dir,"misc/DiMSum_errormodel/processed_data/FOSJUN.txt"),row.names=F,quote=F)

# FOS-JUN cis deepPCA (mutations only in FOS)
all_data = fread(file.path(dataset_dir,"misc/DiMSum_errormodel/datasets/FOScis_2016-06-17_variant_data_merge.tsv"))[Nmut_nt <= 2,
                                                                                             .(Nmut = Nmut_nt,STOP,WT,
                                                                                               input1 = input1_e1_s0_bNA_count,
                                                                                               input2 = input2_e2_s0_bNA_count,
                                                                                               input3 = input3_e3_s0_bNA_count,
                                                                                               output1 = output1_e1_s1_b1_count,
                                                                                               output2 = output2_e2_s1_b1_count,
                                                                                               output3 = output3_e3_s1_b1_count)]
all_data[is.na(WT),WT:=F]
write.table(all_data,file.path(dataset_dir,"misc/DiMSum_errormodel//processed_data/FOScis.txt"),row.names=F,quote=F)


######## GB1 Olson et al. 2014
all_data = fread(file.path(dataset_dir,"misc/DiMSum_errormodel/datasets/Olson2014_GB1_rawcounts.txt"))[,.(Nmut,WT = Nmut==0,Mut,
                                                                                   input1 = input,input2=input,input3=input,
                                                                                   output1,output2,output3)]
write.table(all_data,file.path(dataset_dir,"misc/DiMSum_errormodel/processed_data/GB1.txt"),row.names=F,quote=F)
}