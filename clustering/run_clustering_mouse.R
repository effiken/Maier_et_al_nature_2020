
#Andrew Leader
#AL_C160
# cluster CITEseq samples

library(methods)

###  change on minerva
scClustering_dir <- "scClustering/"
setwd("/sc/orga/work/leadea02/AL_C160_citeseq_mouse_DC_cluster")
###

source(paste(scClustering_dir,"clustering3.r",sep=""))

amp_batches <- as.character(c(158:159))

#load data structure containing count matrices
data_l_fn <- "tmp_data_l.rd"

#load gene lists
mts <- read.csv("mt.csv",stringsAsFactors = F,h=1,r=1)$x
cycle <- read.csv("cycle.csv",stringsAsFactors=F,h=1,r=1)$x

#set gene lists for sorting out cell barcodes we're not interested in
insilico_gating=list()
insilico_gating$MITO=list()
insilico_gating$MITO$genes=mts
insilico_gating$MITO$interval=c(0,0.25)
insilico_gating$PDC=list()
insilico_gating$PDC$genes=strsplit("Sell,Atp1b1,Fcrla,Hsd11b1,Slpi,Gm43291,Cd8b1,Cd8a,Klk1,Siglech,Cox6a2,Bst2,Rpgrip1,Ccr9,Ly6d,Ly6c2,Mzb1",",")[[1]]
insilico_gating$PDC$interval=c(0,0.01)

#clustering set up & function call, for execution on lsf system
k_iter <- 5
model_name=paste("AL_C160_no_cycle_pdc_filt",k_iter,sep="_")
pdf(paste(model_name,".pdf",sep=""))
cluster(data_l_path=data_l_fn,model_name = model_name,k=k_iter,load_seed=F,
        running_mode="LSF_seeding",
        params=list(train_set_size  = 1000,test_set_size = 1000,
                    insilico_gating = insilico_gating,
                    genes_excluded = c(mts),
                    seeding_varmean_quantile=.92,
                    min_umis_per_seeding_gene=20,
                    genes_excluded_from_seeding = c(cycle),
                    init_min_umis_quantile_range=c(0,.2),
                    fixed_downsampling_min_umis=NA,
                    reg=5e-6,
                    n_init_seeds=1000,
                    max_n_cores=1,
                    init_method="TGL_kmeans",
                    km_reg=.1,
                    init_alpha_noise=0.04),
        parallel_params=list(
          account="acc_Meradm01a",
          queue="premium",
          memk=16000,
          wall_time="01:00",
          wall_time_EM="04:00",
          headnode="mothra"
        ))
dev.off()



