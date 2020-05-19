
plot_s2a <- function(){
  
ds <- ldm$dataset$ds[[4]][,ldm$dataset$cell_to_cluster[colnames(ldm$dataset$ds[[4]])]%in%clusts]
ds <- ds[,ldm$dataset$cell_to_sample[colnames(ds)]%in%c("Naive","Naive_CCR7KO")]

clusts <- clusts[order(factor(annots[annots$node%in%clusts,"parent"],c("mregDC","DC1","DC2")))]


cell2cluster <- factor(ldm$dataset$cell_to_cluster,clusts)

cell_ord <- order(cell2cluster[colnames(ds)])
mat <- ds[gene_list,cell_ord]
mat <- log2(1+mat)
mat[mat > 4] <- 4
breaks <- cumsum(table(cell2cluster[colnames(mat)]))/ncol(mat)
breaks <- breaks[breaks<1]
cell2annot <- annots[cell2cluster[colnames(ds)][cell_ord],"parent"]

barcode_2_plot <- colnames(ds)[cell_ord]
s <- seq(1,0,-1/49)

mid_clusts <- (c(0,breaks)+c(breaks,1))/2

cell2sample <- ldm$dataset$cell_to_sample[colnames(ds)][cell_ord]



#RNA truthplot
png(file.path(figure_dir,"figS2a.png"),height=3.18,width=2.53,units="in",res=1000)
layout(matrix(1:6,nrow=3,ncol=2),widths=c(10,4),heights=c(.5,10,.5))

par(oma=c(1,5,1.5,.1))

#column 1
par(mar=c(0,0,0,0))

lin_col <- rgb(t(col2rgb(c(4,5,7)))*.7,max=255)
image(as.matrix(as.numeric(factor(cell2annot))),xaxt="n",yaxt="n",col=lin_col)
box()
abline(v=breaks,col="red")
mtext(at=mid_clusts,c("mregDC","DC1","DC2"),cex=.7,line=.25)
mtext("Annotation",font=2,side=2,las=2,line=.25,cex=.6)
axis(side=3,at=mid_clusts,labels=F,tck=-.2)

image(t(log2(1+as.matrix(mat))),col=rgb(s,s,s),xaxt="n",yaxt="n")
box()
abline(v=breaks,col="red")
mtext(gene_list,side=2,at=seq(0,1,1/(length(gene_list)-1)),las=2,line=.25,font=3,cex=.6)

par(mar=c(0,0,0,0))
image(as.matrix(as.numeric(factor(cell2sample,c("Naive","Naive_CCR7KO")))),col=c(rgb(255,192,203,max=255),rgb(0,178.5,0,max=255)),xaxt="n",yaxt="n")
abline(v=breaks,col="red")
box()
mtext("Sample type",font=2,side=2,las=2,line=.25,cex=.6)


#column 2
par(mar=c(0,0,0,0))
plot.new()
par(mar=c(5,1.5,5,1.5))
image(t(matrix(1:50)),col=rgb(s,s,s),xaxt="n",yaxt="n")
box()
mtext(c("0","\u22654"),at=c(0,1),las=2,side=4,cex=.7,line=.25)
mtext("Log2\n1+#UMI per 2000",side=3,line=.5,cex=.4)
#mtext("1+#UMI/2000",line=.5,cex=.5)
par(mar=c(0,0,0,0))
plot.new()
legend(x = 0,y=2,legend = c("WT",expression(italic("Ccr7"^"-/-"))),fill = c(rgb(255,192,203,max=255),rgb(0,178.5,0,max=255)),
       bty="n",cex=.8,xpd=NA)

dev.off()

}
