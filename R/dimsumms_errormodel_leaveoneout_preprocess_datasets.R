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
  
############ TDP43 Bolognesi, Faure et al. 2019
#290-331 library
all_data = fread(file.path(dataset_dir,"datasets/BB_TARDBP_290_2017-06-13_variant_data_merge.tsv"))[
  Nham_nt < 3 & 
    (input1_e1_s0_bNA_count + input3_e3_s0_bNA_count+ input4_e4_s0_bNA_count) > 10,
 .(WT,
   Nmut = Nham_nt,
   input1 = input1_e1_s0_bNA_count,
   input3 = input3_e3_s0_bNA_count,
   input4 = input4_e4_s0_bNA_count,
   output1 = output1A_e1_s1_b1_count + output1B_e1_s1_b2_count,
   output3 = output3A_e3_s1_b1_count + output3B_e3_s1_b2_count,
   output4 = output4A_e4_s1_b1_count + output4B_e4_s1_b2_count)]
write.table(all_data,file.path(dataset_dir,"processed_data/TDP43_290.txt"),row.names=F,quote=F)

#332-373 library
all_data = fread(file.path(dataset_dir,'datasets/BB_TARDBP_332_2018-04-06_Q25_variant_data_merge.tsv'))[
  Nham_nt <= 2,
 .(Nmut = Nham_nt,
   WT,
   input1 = input1_e1_s0_bNA_count,
   input2 = input2_e2_s0_bNA_count,
   input3 = input3_e3_s0_bNA_count,
   input4 = input4_e4_s0_bNA_count,
   output1 = output1B_e1_s1_b2_count,
   output2 = output2A_e2_s1_b1_count + output2B_e2_s1_b2_count,
   output3 = output3A_e3_s1_b1_count + output3B_e3_s1_b2_count,
   output4 = output4A_e4_s1_b1_count + output4B_e4_s1_b2_count)]
all_data[is.na(WT),WT := FALSE]
write.table(all_data,file.path(dataset_dir,"processed_data/TDP43_332.txt"),row.names=F,quote=F)

############ tRNA
### tRNA phylogeny Domingo et al. 2018
#Heat & Salt
all_data = fread(file.path(dataset_dir,"datasets/JD_Phylogeny_tR-R-CCU_variant_data_merge.tsv"))[,
   .(Nmut = Nham_nt,
     WT,
     input1 = input1A_e1_s0_bNA_count,
     input2 = input1B_e2_s0_bNA_count,
     input3 = input1C_e3_s0_bNA_count,
     input4 = input2A_e4_s0_bNA_count,
     input5 = input2B_e5_s0_bNA_count,
     input6 = input2C_e6_s0_bNA_count,
     output1 = output1A_e1_s1_b1_count,
     output2 = output1B_e2_s1_b1_count,
     output3 = output1C_e3_s1_b1_count,
     output4 = output2A_e4_s1_b1_count,
     output5 = output2B_e5_s1_b1_count,
     output6 = output2C_e6_s1_b1_count)]
all_data[is.na(WT),WT := FALSE]
write.table(all_data,file.path(dataset_dir,"processed_data/tRNA_Phylogeny.txt"),row.names=F,quote=F)

### tRNA Li&Zhang 2018
#23C
all_data = fread(file.path(dataset_dir,"datasets/Li2018_tRNA_sel23_variant_data_merge.tsv"))[
  Nham_nt <= 2 & 
    inputA_e1_s0_bNA_count > 2000 & 
    sel23A_e1_s1_b1_count > 200 & 
    sel23B_e2_s1_b1_count > 200 & 
    sel23C_e3_s1_b1_count > 200 & 
    sel23D_e4_s1_b1_count > 200 & 
    sel23E_e5_s1_b1_count > 200,
 .(Nmut = Nham_nt,
   nt_seq,
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
write.table(all_data,file.path(dataset_dir,"processed_data/tRNA_sel23.txt"),row.names=F,quote=F)

#30C
all_data = fread(file.path(dataset_dir,"datasets/Li2018_tRNA_sel30_variant_data_merge.tsv"))[
  Nham_nt <= 2 & 
    inputA_e1_s0_bNA_count > 2000 & 
    sel30A_e1_s1_b1_count > 200 & 
    sel30B_e2_s1_b1_count > 200 & 
    sel30C_e3_s1_b1_count > 200 & 
    sel30D_e4_s1_b1_count > 200 & 
    sel30E_e5_s1_b1_count > 200,
  .(Nmut = Nham_nt,
    nt_seq,
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
write.table(all_data,file.path(dataset_dir,"processed_data/tRNA_sel30.txt"),row.names=F,quote=F)

#37C
all_data = fread(file.path(dataset_dir,"datasets/Li2018_tRNA_sel37_variant_data_merge.tsv"))[
  Nham_nt <= 2 & 
    inputA_e1_s0_bNA_count > 2000 & 
    sel37A_e1_s1_b1_count > 200 & 
    sel37B_e2_s1_b1_count > 200 & 
    sel37C_e3_s1_b1_count > 200,
 .(Nmut = Nham_nt,
   nt_seq,
   input1 = inputA_e1_s0_bNA_count,
   input2 = inputA_e1_s0_bNA_count,
   input3 = inputA_e1_s0_bNA_count,
   output1 = sel37A_e1_s1_b1_count,
   output2 = sel37B_e2_s1_b1_count,
   output3 = sel37C_e3_s1_b1_count)]
all_data[,WT := nt_seq == "ttccgttggcgtaatggtaacgcgtctccctcctaaggagaagactgcgggttcgagtcccgtacggaa"]
write.table(all_data,file.path(dataset_dir,"processed_data/tRNA_sel37.txt"),row.names=F,quote=F)

#30C, 3%DMSO
all_data = fread(file.path(dataset_dir,"datasets/Li2018_tRNA_selDMSO_variant_data_merge.tsv"))[
  Nham_nt <= 2 & 
    inputA_e1_s0_bNA_count > 2000 & 
    selDMSOA_e1_s1_b1_count > 200 & 
    selDMSOB_e2_s1_b1_count > 200 & 
    selDMSOC_e3_s1_b1_count > 200,
  .(Nmut = Nham_nt,
    nt_seq,
   input1 = inputA_e1_s0_bNA_count,
   input2 = inputA_e1_s0_bNA_count,
   input3 = inputA_e1_s0_bNA_count,
   output1 = selDMSOA_e1_s1_b1_count,
   output2 = selDMSOB_e2_s1_b1_count,
   output3 = selDMSOC_e3_s1_b1_count)]
all_data[,WT := nt_seq == "ttccgttggcgtaatggtaacgcgtctccctcctaaggagaagactgcgggttcgagtcccgtacggaa"]
write.table(all_data,file.path(dataset_dir,"processed_data/tRNA_selDMSO.txt"),row.names=F,quote=F)


########### GRB2 Domingo et al. in preparation
# GPD
all_data = fread(file.path(dataset_dir,"datasets/JD_GRB2_epPCA_stabilityPCA_variant_data_merge.tsv"))[
  Nham_nt <= 2,
 .(Nmut = Nham_nt,
   WT,
   input1 = input1_e1_s0_bNA_count,
   input2 = input2_e2_s0_bNA_count,
   input3 = input3_e3_s0_bNA_count,
   output1 = output1_e1_s1_b1_count,
   output2 = output2_e2_s1_b1_count,
   output3 = output3_e3_s1_b1_count)]
all_data[is.na(WT),WT:=F]
write.table(all_data,file.path(dataset_dir,"processed_data/GRB2_GPD.txt"),row.names=F,quote=F)

# CYC errorpronePCR library bindingPCA assay
all_data = fread(file.path(dataset_dir,"datasets/JD_GRB2_epPCA_bindingPCA_variant_data_merge.tsv"))[
  Nham_nt <= 2,
 .(Nmut = Nham_nt,
   WT,
   input1 = input1_e1_s0_bNA_count,
   input2 = input2_e2_s0_bNA_count,
   input3 = input3_e3_s0_bNA_count,
   output1 = output1_e1_s1_b1_count,
   output2 = output2_e2_s1_b1_count,
   output3 = output3_e3_s1_b1_count)]
all_data[is.na(WT),WT:=F]
write.table(all_data,file.path(dataset_dir,"processed_data/GRB2_CYC.txt"),row.names=F,quote=F)

############ FOS-JUN Diss&Lehner 2018
# FOS-JUN trans deepPCA (mutations in both proteins)
all_data = fread(file.path(dataset_dir,"datasets/GD_FOSJUN_2015-12-15_variant_data_merge.tsv"))[
  WT == TRUE | (between(Nham_aa,1,2) & Nham_aa == Nmut_codons),
 .(Nmut = Nham_aa,
   WT,
   input1 = input1_e1_s0_bNA_count,
   input2 = input2_e2_s0_bNA_count,
   input3 = input3_e3_s0_bNA_count,
   output1 = output1_e1_s1_b1_count,
   output2 = output2_e2_s1_b1_count,
   output3 = output3_e3_s1_b1_count)]
all_data[is.na(WT),WT:=F]
write.table(all_data,file.path(dataset_dir,"processed_data/FOSJUN.txt"),row.names=F,quote=F)

# FOS cis deepPCA (mutations only in FOS)
all_data = fread(file.path(dataset_dir,"datasets/GD_FOSintra_2016-06-17_variant_data_merge.tsv"))[
  Nham_nt <= 2,
 .(Nmut = Nham_nt,
   WT,
   input1 = input1_e1_s0_bNA_count,
   input2 = input2_e2_s0_bNA_count,
   input3 = input3_e3_s0_bNA_count,
   output1 = output1_e1_s1_b1_count,
   output2 = output2_e2_s1_b1_count,
   output3 = output3_e3_s1_b1_count)]
all_data[is.na(WT),WT:=F]
write.table(all_data,file.path(dataset_dir,"processed_data/FOScis.txt"),row.names=F,quote=F)

######## GB1 Olson et al. 2014
all_data = fread(file.path(dataset_dir,"datasets/Olson2014_GB1_rawcounts.txt"))[,.(Nmut,WT = Nmut==0,Mut,
                                                                                   input1 = input,input2=input,input3=input,
                                                                                   output1,output2,output3)]
write.table(all_data,file.path(dataset_dir,"processed_data/GB1.txt"),row.names=F,quote=F)

}