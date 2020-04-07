
setwd("~/GoogleDrive/work/scClustering/example/")
source("~/GoogleDrive/work//scClustering/clustering2.r")

rps=unique(read.csv(file="~/GoogleDrive/work/scClustering/gene_lists/ribsomial_genes.csv",stringsAsFactors = F)[,1])
malat=c("MALAT1")
mts=c('MT-ND1','MT-ND2','MT-CO1','MT-CO2','MT-ATP8','MT-ATP6','MT-CO3','MT-ND3','MT-ND4L','MT-ND4','MT-ND5','MT-ND6','MT-CYB')

insilico_gating=list()

insilico_gating$MITO=list()
insilico_gating$MITO$genes=mts
insilico_gating$MITO$interval=c(0,0.1)



libnames=c("Fresh_pbmc_68k","pbmc_33k_v1")
data_l=read_multiple_mtx("~/GoogleDrive/work/shared/data/public/human/10x_data/",libnames,cell_interval = c(500,25000),noise_interval = c(100,25000))
cluster(data_l,train_set_size  = 500,seed_test_set_size = 1000,model_name = "pbmc_v1",k=15,insilico_gating = insilico_gating,genes_excluded = c(mts,malat),var_mean_seeding_thresh=.5,min_umis_per_seeding_gene=30,genes_excluded_from_seeding = c(rps),n_seeding_attempts=100,seeding_min_umis_quantile_range=c(0,.2),load_seed = F,reg=1e-7)


libnames=c("pbmc_33k_v1","cd_34")
data_l=read_multiple_mtx("~/GoogleDrive/work/shared/data/public/human/10x_data/",libnames,cell_interval = c(500,25000),noise_interval = c(100,25000))
cluster(data_l,train_set_size  = 500,seed_test_set_size = 1000,model_name = "pbmc_cd34_v1",k=15,insilico_gating = insilico_gating,genes_excluded = c(mts,malat),var_mean_seeding_thresh=.5,min_umis_per_seeding_gene=30,genes_excluded_from_seeding = c(rps),n_seeding_attempts=100,seeding_min_umis_quantile_range=c(0,.2),load_seed = F,reg=1e-7)
