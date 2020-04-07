scDissector_dir="~/Documents/GitHub/scDissector/"
source(paste(scDissector_dir,"load_dataset.R",sep=""))
source(paste(scDissector_dir,"projector.R",sep=""))

scDissector_datadir="~/GoogleDrive/work/shared/clustering_data/clustering_data_public/"
model_fn=paste(scDissector_datadir,"model_pbmc_v1_30.rd",sep="/")
samples=read.csv(paste(scDissector_datadir,"/samples.csv",sep="/"),stringsAsFactors = F,row.names = 1)
sample_names=rownames(samples)
samples_fn=paste(scDissector_datadir,samples[sample_names,"path"],sep="/")
names(samples_fn)=sample_names                   

default_model_dataset=load_dataset_and_model(model_fn,samples_fn,500,"pbmc_v1")
library(shiny)
runApp(scDissector_dir)