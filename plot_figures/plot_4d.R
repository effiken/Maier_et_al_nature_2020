
#load("/users/andrew leader/google drive/merad/scRNAseq_analysis/clustering_metadata/data_l_lungDC_190819.rd")
load("/users/andrew leader/google drive/merad/scRNAseq_analysis/clustering_metadata/ldm_lung_190813_dc.rd")


annots_list <- read.csv("/users/andrew leader/google drive/merad/scRNAseq_analysis/compiled/clustering_data_lung4/model_lung_190606_annot_lists.csv",
                        r=1,h=1,stringsAsFactors=F)
cell_to_annot <- annots_list[as.character(lung_ldm$dataset$cell_to_cluster),]$sub_lineage
names(cell_to_annot) <- names(lung_ldm$dataset$cell_to_cluster)
s <- split(names(cell_to_annot),cell_to_annot)
exprs <- lapply(s,function(x){rs <- rowSums(lung_ldm$dataset$umitab[,x]); return(rs/sum(rs))})

reg <- 1e-6
dc2_vs_dc1 <- log2((reg+exprs$cDC2)/(reg+exprs$cDC1))
mreg_vs_dc1 <- log2((reg+exprs$mregDC)/(reg+exprs$cDC1))
mreg_vs_dc2 <- log2((reg+exprs$mregDC)/(reg+exprs$cDC2))

gene_sets <- list()
gene_sets$dc2 <- names(dc2_vs_dc1)[dc2_vs_dc1 > 1 & exprs$cDC2>1e-5]
gene_sets$dc1 <- names(dc2_vs_dc1)[dc2_vs_dc1 < -1 & exprs$cDC1>1e-5]
gene_sets$mreg <- names(dc2_vs_dc1)[mreg_vs_dc1 > 1.5 & mreg_vs_dc2 > 1.5 & exprs$mregDC>1e-5]

numi <- colSums(lung_ldm$dataset$umitab)

gene_scores <- lapply(gene_sets,function(x){log10(1e-4+colSums(lung_ldm$dataset$umitab[x,])/numi)})

cells <- names(cell_to_annot)[cell_to_annot%in%c("cDC1","cDC2","mregDC")]
annot2col <- rgb(t(col2rgb(c(7,4,5))*.7/255)); names(annot2col) <- c("mregDC","cDC1","cDC2")

colvec <- log10(1e-4+lung_ldm$dataset$umitab["IL12B",cells]/numi[cells])
colvec <- round((colvec-min(colvec))/diff(range(colvec))*49+1)
colvec <- rev(magma(50))[colvec]

xlim=c(-0.9,1.4); ylim=c(-2.8,-0.9)

png("/users/andrew leader/Dropbox/Andrew&Barbara/figure_R_scripts/fig4/fig4d_191028.png",height=1.4,width=4.49,
    units="in",res=1000,pointsize=5)
layout(matrix(1:2,nrow=1))
par(bg=NA,bty="n",pin=c(1.1,0.9))
plot(gene_scores$dc2[cells]-gene_scores$dc1[cells],gene_scores$mreg[cells],pch=16,col=alpha(annot2col[cell_to_annot[cells]],0.6),
     xlim=xlim,ylim=ylim,xlab="",ylab="")

plot(gene_scores$dc2[cells]-gene_scores$dc1[cells],gene_scores$mreg[cells],pch=16,col=alpha(colvec,0.6),
     xlim=xlim,ylim=ylim,xlab="",ylab="")
dev.off()


