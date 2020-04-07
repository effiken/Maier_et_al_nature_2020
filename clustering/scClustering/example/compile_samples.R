setwd("~/GoogleDrive/work/scClustering/example/")
source("../compile_data.r")

datain_path="~/GoogleDrive/work/shared/data/public/human/10x_data/"
libs=list.dirs(datain_path,recursive = F,full.names = F)
output_path="~/GoogleDrive/work/shared/clustering_data/clustering_data_public/"
compile_data(input_path = datain_path,output_path = output_path,libnames =libs)

