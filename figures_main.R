
library(Matrix)

rm(list=ls())

############################
#Specify a working directory. This is the path where the data will download to, and where figures will be produced.
# It should also be the path to the scripts from the github.
wd <- getwd()
############################

setwd(wd)
scClustering_dir <- "clustering/scClustering/"
source(file.path(scClustering_dir,"DE.r"))
source(file.path(scClustering_dir,"clustering3.r"))
source("clustering/adt_list_to_matrix.R")


# set up the directory
if(!dir.exists(wd)){
  dir.create(wd)
}
data_dir <- file.path(wd,"data")
if(!dir.exists(data_dir)){
  dir.create(data_dir)
}
figure_dir <- file.path(wd,"figures_out")
if(!dir.exists(figure_dir)){
  dir.create(figure_dir)
}

# to download & load mouse DC data
data_url <- "https://www.dropbox.com/s/3ffq3z75af37rgr/mouse_dc.rd?dl=1"
if(Sys.info()["sysname"]=="Windows"){
  download.file(url = data_url,destfile = file.path(data_dir,"mouse_dc.rd"),mode="wb")
}else{
  download.file(url = data_url,destfile = file.path(data_dir,"mouse_dc.rd"))
}
download_check <- try(load(file.path(data_dir,"mouse_dc.rd")))
if(download_check!="ldm"){stop("Error using download.file(), due to operating system?-- please manually download human DC dataset from https://www.dropbox.com/s/3ffq3z75af37rgr/mouse_dc.rd?dl=1 and place in working directory and remove this chunk of code")}


gene_list <- rev(strsplit("Ccr7,Fscn1,Cd40,Cd274,Pdcd1lg2,Cd200,Fas,Aldh1a2,Il4ra,Il4i1,Socs2,Relb,Xcr1,Clec9a,Cadm1,Naaa,Cd209a,Sirpa,H2-DMb2,Itgam",",")[[1]])
adt_mat <- adt_list_to_matrix(ldm$dataset$adt_by_sample)
rownames(adt_mat) <- c("CD103","CD11b","CD11c","I-A/I-E","XCR1")
adt_mat <- adt_mat[c(4,3,2,1,5),]
annots <- data.frame("node"=c("5","2","4","1","3"),"parent"=c("mregDC","mac&contam","DC1","mac&contam","DC2"))
clusts <- annots[annots$parent%in%c("mregDC","DC1","DC2"),"node"]
ds <- ldm$dataset$ds[[4]][,ldm$dataset$cell_to_cluster[colnames(ldm$dataset$ds[[4]])]%in%clusts]
ds <- ds[,ldm$dataset$cell_to_sample[colnames(ds)]%in%c("Naive","KP")]
clusts <- clusts[order(factor(annots[annots$node%in%clusts,"parent"],c("mregDC","DC1","DC2")))]

#plotting 

source("plot_figures/plot_1beg.R")
plot_1beg()

source("plot_figures/plot_1c.R")
plot_1c()

source("plot_figures/plot_1d.R")
plot_1d()

source("plot_figures/plot_1f.R")
plot_1f()

source("plot_figures/plot_1h.R")
plot_1h()

source("plot_figures/plot_s1a.R")
plot_s1a()

source("plot_figures/plot_2a.R")
plot_2a()

source("plot_figures/plot_2e.R")
plot_2e()

source("plot_figures/plot_s2a.R")
plot_s2a()

source("plot_figures/plot_s2b.R")
plot_s2b()

source("plot_figures/plot_s2i.R")
plot_s2i()


source("plot_figures/plot_3h.R")
plot_3h()

source("plot_figures/plot_s3e.R")
plot_s3e()

source("plot_figures/plot_s3f.R")
plot_s3f()



#################
# # to download & load human DC counts and metadata
data_url <- "https://www.dropbox.com/s/e0cvefqzoy4hz0r/human_dc.rd?dl=1"
if(Sys.info()["sysname"]=="Windows"){
  download.file(url = data_url,destfile = file.path(wd,"human_dc.rd"),mode="wb")
}else{
  download.file(url = data_url,destfile = file.path(wd,"human_dc.rd"))
}
download_check <- try(load(file.path(wd,"human_dc.rd")))
if(download_check!="human_dc"){stop("Error using download.file(), due to operating system?-- please manually download human DC dataset from https://www.dropbox.com/s/e0cvefqzoy4hz0r/human_dc.rd?dl=1 and place in working directory and remove this chunk of code")}
###############

source("plot_figures/plot_4a.R")
plot_4a()

#4b & 4c
source("plot_figures/plot_4b.R")
plot_4b()

source("plot_figures/plot_4e_s4a.R")
plot_4e_s4a()

source("plot_figures/plot_4f.R")
plot_4f()

#################
# # to download & load human DC counts and metadata
data_url <- "https://www.dropbox.com/s/owhksv8jsbtci3o/lambrechts_dc.rd?dl=1"
if(Sys.info()["sysname"]=="Windows"){
  download.file(url = data_url,destfile = file.path(wd,"lambrechts_dc.rd"),mode="wb")
}else{
  download.file(url = data_url,destfile = file.path(wd,"lambrechts_dc.rd"))
}
download_check <- try(load(file.path(wd,"lambrechts_dc.rd")))
if(download_check!="lambrechts_dc"){stop("Error using download.file(), due to operating system?-- please manually download human DC dataset from https://www.dropbox.com/s/owhksv8jsbtci3o/lambrechts_dc.rd?dl=1 and place in working directory and remove this chunk of code")}
###############

source("plot_figures/plot_s4b.R")
plot_s4b()
