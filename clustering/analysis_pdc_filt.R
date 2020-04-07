rm(list=ls())
library(scDissector)

#####################
# CHANGE TO APPROPRIATE PATHS

setwd("/users/andrew leader/Dropbox/Andrew&Barbara/figure_R_scripts/data/") 

# where do the compiled mouse samples live on google drive?
compiled_dir <- "/users/andrew leader/google drive/merad/scRNAseq_analysis/compiled/Mouse_compiled/"

scDissector_dir <- "/users/andrew leader/Documents/GitHub/scDissector/"
sample_annots_path <- "/users/andrew leader/Google Drive/merad/scRNAseq_analysis/sc_data_main/sample_annots.csv"

#####################


#   defining load_dataset_and_model() inline by sourcing the script, which lets me use a  
#   modified version of the update_alpha_single_batch() function with max alpha very small,
#   defined within the projector.R script in the local directory
source(paste(scDissector_dir,"R/load_dataset_and_model.R",sep=""))
source("scripts/projector.R")


sample_ids <- c(209,210,213,214,816,815)
samples_fn <- paste(compiled_dir,"Mousedata_",sample_ids,".rd",sep="")
names(samples_fn) <- c("Naive","KP","Naive_CCR7KO","KP_CCR7KO","Naive_B16","B16")

#K=10
# model_fn <- "models/model_AL_C160_no_cycle_10.rd"
# ldm <- load_dataset_and_model(model_fn,samples_fn,min_umis=800)
# run_scDissector(ldm)

#K=5
model_fn <- "models/model_AL_C160_no_cycle_pdc_filt_5.rd"
ldm <- load_dataset_and_model(model_fn,samples_fn,min_umis=800)
#run_scDissector(ldm)
annots <- read.delim("/users/andrew leader/Dropbox/Andrew&Barbara/figure_R_scripts/data/models/model_AL_C160_no_cycle_pdc_filt_5_cluster_sets.txt",sep="\t")
