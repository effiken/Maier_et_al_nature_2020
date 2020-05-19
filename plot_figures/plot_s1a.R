
library(scDissector)
library(matrixStats)

plot_s1a <- function(){
  
ds <- ldm$dataset$ds[[4]]
ds <- ds[,ldm$dataset$cell_to_sample[colnames(ds)]%in%c("Naive","KP")]
  
clust_ord <- c(5,4,3,2,1)

cell_ord <- order(factor(ldm$dataset$cell_to_cluster[colnames(ds)],clust_ord),
                  factor(ldm$dataset$cell_to_sample[colnames(ds)],c("Naive","KP")))



ds <- ds[,cell_ord]

gene_list <- rev(strsplit("Csf2ra,Csf2rb,Flt3,Cd24a,Itgax,Zbtb46,Csf1r,Cd14,C1qa,C1qb,C1qc,Mertk,Fcgr1,Ly6c2,Sftpc,Pecam1,Dcn,Cd3d,Cd19,Ncr1,Ccr7,Fscn1,Cd40,Cd274,Pdcd1lg2,Cd200,Fas,Aldh1a2,Il4ra,Il4i1,Socs2,Relb",",")[[1]])

mat <- ds[gene_list,]
mat <- log2(1+mat)
mat[mat > 4] <- 4
breaks <- which(diff(as.integer(ldm$dataset$cell_to_cluster[colnames(mat)]))!=0)/ncol(mat)


cell2annot <- annots[match(ldm$dataset$cell_to_cluster[colnames(mat)],annots$node),"parent"]

cell2annot <- factor(cell2annot,c("mregDC","DC1","DC2","mac&contam"))

s <- seq(1,0,-1/49)

mid_clusts <- (c(0,breaks)+c(breaks,1))/2

cell2sample <- ldm$dataset$cell_to_sample[colnames(mat)]



#RNA truthplot
png(file.path(figure_dir,"fig_s1a.png"),height=2.36,width=3.24,units="in",res=1000)
layout(matrix(1:6,nrow=3,ncol=2),widths=c(10,4),heights=c(.5,10,.5))

par(oma=c(1,5,1.5,.1))

#column 1
par(mar=c(0,0,0,0))

lin_col <- rgb(t(col2rgb(c(7,4,5,6)))*.7,max=255)
image(as.matrix(as.numeric(cell2annot)),xaxt="n",yaxt="n",col=lin_col)
box()
abline(v=breaks,col="red")
mtext("Annotation",font=2,side=2,las=2,line=.25,cex=.6)

image(t(log2(1+as.matrix(mat))),col=rgb(s,s,s),xaxt="n",yaxt="n")
box()
abline(v=breaks,col="red")
mtext(gene_list,side=2,at=seq(0,1,1/(length(gene_list)-1)),las=2,line=.25,font=3,cex=.4)
abline(h=seq(-0.5,31.5,1)[13]/31,col="red",lty=2,lwd=0.75)
 
par(mar=c(0,0,0,0))
image(as.matrix(as.numeric(factor(cell2sample,c("Naive","KP")))),col=c(rgb(255,192,203,max=255),rgb(165,75,42,max=255)),xaxt="n",yaxt="n")
abline(v=breaks,col="red")
box()
mtext("Sample type",font=2,side=2,las=2,line=.25,cex=.6)


#column 2
par(mar=c(0,0,0,0))
plot.new()
par(mar=c(4,2.5,4,2.5))
image(t(matrix(1:50)),col=rgb(s,s,s),xaxt="n",yaxt="n")
box()
mtext(c("0","\u22654"),at=c(0,1),las=2,side=4,cex=.7,line=.25)
mtext("Log2\n1+#UMI per 2000",side=3,line=.5,cex=.4)
par(mar=c(0,0,0,0))
plot.new()
legend(x = 0,y=2,legend = c("Naive","KP"),fill = c(rgb(255,192,203,max=255),rgb(165,75,42,max=255)),
      bty="n",cex=.8,xpd=NA)

dev.off()
}
