
rm(list=ls())

# set up
sample_annots <- read.csv('/users/andrew leader/google drive/merad/scRNAseq_analysis/sc_data_main/sample_annots.csv',
                          r=1,h=1,stringsAsFactors=F)

scClustering_dir <- "/users/andrew leader/Documents/GitHub/scClustering/"
source("/users/andrew leader/Documents/GitHub/scClustering/DE.r")
source("/users/andrew leader/Documents/GitHub/scClustering/clustering3.r")

source("/users/andrew leader/Dropbox/Andrew&Barbara/figure_R_scripts/data/analysis_pdc_filt.R")

setwd("/users/andrew leader/Dropbox/Andrew&Barbara/figure_R_scripts/make_figures/")
gene_list <- rev(strsplit("Ccr7,Fscn1,Cd40,Cd274,Pdcd1lg2,Cd200,Fas,Aldh1a2,Il4ra,Il4i1,Socs2,Relb,Xcr1,Clec9a,Cadm1,Naaa,Cd209a,Sirpa,H2-DMb2,Itgam",",")[[1]])
source("/users/andrew leader/Documents/GitHub/scDissector/R/adt_list_to_matrix.R")
adt_mat <- adt_list_to_matrix(ldm$dataset$adt_by_sample)
rownames(adt_mat) <- c("CD103","CD11b","CD11c","I-A/I-E","XCR1")
adt_mat <- adt_mat[c(4,3,2,1,5),]
clusts <- annots[annots$parent%in%c("mregDC","DC1","DC2"),"node"]
ds <- ldm$dataset$ds[[4]][,ldm$dataset$cell_to_cluster[colnames(ldm$dataset$ds[[4]])]%in%clusts]
ds <- ds[,ldm$dataset$cell_to_sample[colnames(ds)]%in%c("Naive","KP")]
clusts <- clusts[order(factor(annots[annots$node%in%clusts,"parent"],c("mregDC","DC1","DC2")))]

#plotting 

source("plot_1beg.R")
plot_1beg()

source("plot_1c.R")
plot_1c()

source("plot_1d.R")
plot_1d()

source("plot_1f.R")
plot_1f()

source("plot_1h.R")
plot_1h()

source("plot_s1a.R")
plot_s1a()

source("plot_2a.R")
plot_2a()

source("plot_2e.R")
plot_2e()

source("plot_s2a.R")
plot_s2a()

source("plot_s2b.R")
plot_s2b()

source("plot_s2i.R")
plot_s2i()


source("plot_3h.R")
plot_3h()

source("plot_s3e.R")
plot_s3e()

source("plot_s3f.R")
plot_s3f()
## load human

hum <- new.env()
load(envir=hum,"/users/andrew leader/google drive/merad/scRNAseq_analysis/clustering_metadata/ldm_lung_190813_dc.rd")

source("plot_4a.R")
plot_4a()

#4b & 4c
source("plot_4b.R")
plot_4b()

source("plot_4e_s4a.R")
plot_4e_s4a()

source("plot_4f.R")
plot_4f()

source("plot_s4b.R")
plot_s4b()
