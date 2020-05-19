plot_4f <- function(){

#get hvg across DCs in each species
mus_mask <- colnames(ldm$dataset$umitab)[ldm$dataset$cell_to_cluster%in%c("5","4","3") & ldm$dataset$cell_to_sample%in%c("KP","Naive")]
mus_mask <- intersect(mus_mask,colnames(ldm$dataset$ds[[4]]))
hvg_mus <- get_highly_variable_genes(ds = ldm$dataset$ds[[4]][,mus_mask],
                                     seeding_genes = rownames(ldm$dataset$ds[[4]]),
                                     min_umis_per_var_gene = 300,varmean_quantile=.8)

hum_mask <- colnames(human_dc$filtered_ds)
hvg_hum <- get_highly_variable_genes(ds=human_dc$filtered_ds,
                                     seeding_genes = rownames(human_dc$filtered_ds),
                                     min_umis_per_var_gene = 300,varmean_quantile=.8)


overlap <- intersect(hvg_hum,toupper(hvg_mus))

#generate avg human signatures
s <- split(colnames(human_dc$filtered_umitab),human_dc$cell_to_annot)
exprs <- lapply(s,function(x){rs <- rowSums(human_dc$filtered_umitab[,x]); return(rs/sum(rs))})
exprs <- do.call(cbind,exprs)
exprs <- exprs[,c("mregDC","DC1","DC2")]

reg <- 1e-4
exprs <- log2((exprs+reg)/rowMeans(exprs+reg))

s <- split(colnames(ldm$dataset$umitab),ldm$dataset$cell_to_cluster)
mouse_exprs <- lapply(s,function(x){rs <- rowSums(ldm$dataset$umitab[,x]); return(rs/sum(rs))})
mouse_exprs <- do.call(cbind,mouse_exprs)
mouse_exprs <- mouse_exprs[,c("5","4","3")]
colnames(mouse_exprs) <- c("mDC","DC1","DC2")

mouse_exprs <- log2((mouse_exprs+reg)/rowMeans(mouse_exprs+reg))
rownames(mouse_exprs) <- toupper(rownames(mouse_exprs))

cross_species_mat <- cbind(exprs[overlap,],mouse_exprs[overlap,])
colnames(cross_species_mat) <- c("h-mregDC","h-DC1","h-DC2","m-mregDC","m-DC1","m-DC2")

samp_h <- hclust(as.dist((1-cor(cross_species_mat))/2))
#genes_h <- hclust(as.dist((1-cor(t(cross_species_mat)))/2))


##############
#seed 4; k=5 gave great results.
set.seed(9)
genes_k <- kmeans(cross_species_mat,5)



celltype_score <- matrix(nrow=nrow(cross_species_mat),ncol=3,dimnames=list(rownames(cross_species_mat),c("mDC","DC1","DC2")))
for(iter in 1:3){
  celltype_score[,iter] <- rowMeans(cross_species_mat[,c(iter,iter+3)])
}
for(iter in 1:nrow(celltype_score)){
  celltype_score[iter,] <- sort(celltype_score[iter,])
}

diff_ord <- rownames(celltype_score)[order(celltype_score[,3])]

gene_ord <- diff_ord[order(factor(genes_k$cluster[diff_ord],c("4","2","1","3","5")))]




pdf(file.path(figure_dir,"fig4f.pdf"),width=6.63,height=2.37)
layout(matrix(1:9,nrow=3),widths=c(2,7,1.7),heights=c(1,.25,4))

par(mar=c(0,0,0,0),oma=c(4,4,.1,.1))

#col 1
plot.new()
plot.new()
par(mar=c(0,0,0,6.5))
plot(as.dendrogram(samp_h),horiz=TRUE,xlab="sqrt((1-cor)/2)",type="triangle",yaxs="i",axes=F,leaflab="none")
axis(side=1,at=c(0,.5,1),cex.axis=0.6,xpd="",mgp=c(1.5,.5,0))
title(xlab=expression("(1-"*italic("cor")*")/2"),cex=.5,xpd="",line=1.5)


#col 2
par(mar=c(0,0,0,0))
plot.new()
image(as.matrix(as.numeric(genes_k$cluster[gene_ord])),col=rgb(.7*t(col2rgb(c(7,4,7,5,6))),max=255),xaxt="n",yaxt="n")
mtext(side=2,las=2,"Kmeans cluster:",line=.5,font=2,cex=.75)
breaks <- which(diff(genes_k$cluster[gene_ord])!=0)/length(gene_ord)
box()
abline(v=breaks)
mtext("226 conserved, variable DC genes",line=.5)


mat <- cross_species_mat[gene_ord,samp_h$order]
thresh <- 1
mat[mat > thresh] <- thresh; mat[mat < -thresh] <- -thresh
image(mat,col=bluered(50),xaxt="n",yaxt="n")
box()
abline(v=breaks)
mtext(c(expression("mregDC"[italic("H. sap.")]),expression("mregDC"[italic("M. mus.")]),
        expression("DC1"[italic("H. sap.")]),expression("DC1"[italic("M. mus.")]),
      expression("DC2"[italic("H. sap.")]),expression("DC2"[italic("M. mus.")])),side=2,las=2,at=seq(0,1,1/(5)),cex=.75,adj=.5,line=3)


plot.new()
plot.new()

par(mar=c(3,3.2,3,3.2))
image(as.matrix(t(1:50)),col=bluered(50),xaxt="n",yaxt="n")
box()
mtext(c("\u2264 -1","\u2265 1"),at=c(0,1),side=4,las=2,line=.25,cex=.7)
mtext("Log2(Expression\n/Species Avg)",side=3,cex=.7,line=.5)



dev.off()



}
