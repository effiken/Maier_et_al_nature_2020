
plot_1beg <- function(){
  

cell2cluster <- factor(ldm$dataset$cell_to_cluster,clusts)

cell_ord <- order(cell2cluster[colnames(ds)])
mat <- ds[gene_list,cell_ord]
mat <- log2(1+mat)
mat[mat > 4] <- 4
breaks <- cumsum(table(cell2cluster[colnames(mat)]))/ncol(mat)
breaks <- breaks[breaks<1]
cell2annot <- annots[match(cell2cluster[colnames(ds)][cell_ord],annots$node),"parent"]
cell2annot <- factor(cell2annot,c("mregDC","DC1","DC2","mac&contam"))

barcode_2_plot <- colnames(ds)[cell_ord]
s <- seq(1,0,-1/49)

mid_clusts <- (c(0,breaks)+c(breaks,1))/2

cell2sample <- ldm$dataset$cell_to_sample[colnames(ds)][cell_ord]

lin_col <- rgb(t(col2rgb(c(7,4,5)))*.7,max=255)


#RNA truthplot
png(file.path(figure_dir,"fig1b.png"),height=3.18,width=2.53,units="in",res=1000)
layout(matrix(1:6,nrow=3,ncol=2),widths=c(10,4),heights=c(.5,10,.5))

par(oma=c(1,5,1.5,.1))

#column 1
par(mar=c(0,0,0,0))

image(as.matrix(as.numeric(cell2annot)),xaxt="n",yaxt="n",col=lin_col)
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
image(as.matrix(as.numeric(factor(cell2sample,c("Naive","KP")))),col=c(rgb(255,192,203,max=255),rgb(165,75,42,max=255)),xaxt="n",yaxt="n")
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
legend(x = 0,y=2,legend = c("Naive","KP"),fill = c(rgb(255,192,203,max=255),rgb(165,75,42,max=255)),
       bty="n",cex=.8,xpd=NA)

dev.off()

#CITEseq dot plot
pdf(file.path(figure_dir,"fig1e.pdf"),height=2.21,width=3.78,pointsize = 12)
#png("fig1e.png",height=2.21,width=3.78,units="in",res=3000)
layout(matrix(1:5,nrow=5),heights=c(array(1,5)))

par(oma=c(.2,4,2,3))
par(mar=c(0,0,0,0))

for(iter in 1:nrow(adt_mat)){
  plot(log2(1+adt_mat[iter,barcode_2_plot]),
       xaxs="i",pch=".",cex=2,
       col=lin_col[as.matrix(as.numeric(cell2annot))],
       xaxt="n",tcl=-.15,yaxt="n")
  axis(side=4,at=floor(range(log2(1+adt_mat[iter,barcode_2_plot])))+c(1,-1),las=2,cex=.5,tcl=-.2,hadj=.75)
  abline(v=breaks*length(barcode_2_plot),col="red")
  

  mtext(rownames(adt_mat)[iter],side=2,line=.5,las=1,cex=.7)
 # mtext(floor(range(log2(1+adt_mat[iter,barcode_2_plot])))+c(1,-1),side=4,
#        at=floor(range(log2(1+adt_mat[iter,barcode_2_plot])))+c(1,-1),las=2,cex=.5,line=.25)
  
  if(iter==1){
    mtext(at=mid_clusts*length(cell_ord),c("mregDC","DC1","DC2"),cex=.7,line=.25)
    axis(side=3,at=mid_clusts*length(cell_ord),labels=F,tck=-.075)
    
  }
  
  if(iter==3){
    mtext(expression("Log"[2]*"(1+#ADT)"),side=4,line=1.5,cex=.7,xpd=NA)
  }
}


dev.off()


adt_mat <- adt_mat[,ldm$dataset$cell_to_sample[colnames(adt_mat)]%in%c("Naive","KP")]

lin_col <- rgb(t(col2rgb(c(7,4,5))*.7),max=255)
names(lin_col) <- annots$node[1:3]

png(file.path(figure_dir,"fig1g.png"),height=2.5,width=3,units="in",res=300,pointsize=8)
par(oma=c(0,0,0,0),mar=c(5,2.5,5,1),cex.axis=.5,mgp=c(1,.25,0),tcl=-.15,cex=.5)

layout(matrix(1:2,nrow=1))

mask <- ldm$dataset$cell_to_cluster[colnames(adt_mat)]%in%clusts &
  ldm$dataset$cell_to_sample[colnames(adt_mat)]=="Naive"
plot(log2(1+adt_mat["CD11c",mask]),(log2(1+adt_mat["I-A/I-E",mask])),bty="n",col="white",xlab="CD11c",ylab="I-A/I-E",xlim=c(0,12))
mtext("Naive")
mask <- ldm$dataset$cell_to_cluster[colnames(adt_mat)]==clusts[3] & ldm$dataset$cell_to_sample[colnames(adt_mat)]=="Naive"
points(log2(1+adt_mat["CD11c",mask]),(log2(1+adt_mat["I-A/I-E",mask])),pch=20,col=paste(lin_col[cell2cluster[colnames(adt_mat)[mask]]],"2D",sep=""),cex=.5)
mask <- ldm$dataset$cell_to_cluster[colnames(adt_mat)]==clusts[2] & ldm$dataset$cell_to_sample[colnames(adt_mat)]=="Naive"
points(log2(1+adt_mat["CD11c",mask]),(log2(1+adt_mat["I-A/I-E",mask])),pch=20,col=paste(lin_col[cell2cluster[colnames(adt_mat)[mask]]],"2D",sep=""),cex=.5)
mask <- ldm$dataset$cell_to_cluster[colnames(adt_mat)]==clusts[1] & ldm$dataset$cell_to_sample[colnames(adt_mat)]=="Naive"
points(log2(1+adt_mat["CD11c",mask]),(log2(1+adt_mat["I-A/I-E",mask])),pch=20,col=paste(lin_col[cell2cluster[colnames(adt_mat)[mask]]],"2D",sep=""),cex=.5)

mask <- ldm$dataset$cell_to_cluster[colnames(adt_mat)]%in%clusts &
  ldm$dataset$cell_to_sample[colnames(adt_mat)]=="KP"
plot(log2(1+adt_mat["CD11c",mask]),(log2(1+adt_mat["I-A/I-E",mask])),bty="n",col="white",xlab="CD11c",ylab="I-A/I-E",xlim=c(0,12))
mtext("Tumor")
mask <- ldm$dataset$cell_to_cluster[colnames(adt_mat)]==clusts[3] & ldm$dataset$cell_to_sample[colnames(adt_mat)]=="KP"
points(log2(1+adt_mat["CD11c",mask]),(log2(1+adt_mat["I-A/I-E",mask])),pch=20,col=paste(lin_col[cell2cluster[colnames(adt_mat)[mask]]],"2D",sep=""),cex=.5)
mask <- ldm$dataset$cell_to_cluster[colnames(adt_mat)]==clusts[2] & ldm$dataset$cell_to_sample[colnames(adt_mat)]=="KP"
points(log2(1+adt_mat["CD11c",mask]),(log2(1+adt_mat["I-A/I-E",mask])),pch=20,col=paste(lin_col[cell2cluster[colnames(adt_mat)[mask]]],"2D",sep=""),cex=.5)
mask <- ldm$dataset$cell_to_cluster[colnames(adt_mat)]==clusts[1] & ldm$dataset$cell_to_sample[colnames(adt_mat)]=="KP"
points(log2(1+adt_mat["CD11c",mask]),(log2(1+adt_mat["I-A/I-E",mask])),pch=20,col=paste(lin_col[cell2cluster[colnames(adt_mat)[mask]]],"2D",sep=""),cex=.5)

dev.off()
}
