library(scales)

plot_4e_s4a = function(){
  

# Load the CITEseq matrices

cells <- colnames(human_dc$filtered_umitab)

adtmat <- adt_list_to_matrix(human_dc$adt_matrix_by_sample)

# Figures 4e

clust2col <- rgb(t(col2rgb(c(7,4,5))*.7/255))
names(clust2col) <- c("mregDC","DC1","DC2")
clust2col[2:3] <- alpha(clust2col[2:3],.2)

png(file.path(figure_dir,"fig4e.png"),height=1.61,width=4.5,units="in",res=300,pointsize=6)
layout(matrix(1:3,nrow=1))
par(cex.lab=1.5,pin=c(1,1),xpd=NA,mgp=c(1.9,0.75,0))
plot(log2(1+adtmat["CD1c",]),log2(1+adtmat["CD141",]),
     pch=16,
     col=clust2col[human_dc$cell_to_annot[colnames(adtmat)]],
     bty="n",xlab="CD1c",ylab="CD141")

plot(log2(1+adtmat["CD11c",]),log2(1+adtmat["HLADR",]),
     pch=16,
     col=clust2col[human_dc$cell_to_annot[colnames(adtmat)]],
     bty="n",xlab="CD11c",ylab="HLADR")

plot(log2(1+adtmat["CD274PDL1",]),log2(1+adtmat["CD40",]),
     pch=16,
     col=clust2col[human_dc$cell_to_annot[colnames(adtmat)]],
     bty="n",xlab="PD-L1",ylab="CD40")
dev.off()




# Figure S4A
s <- split(colnames(adtmat),list(
           factor(human_dc$cell_to_annot[colnames(adtmat)],c("mregDC","DC1","DC2")),
           factor(human_dc$cell_to_tissue[colnames(adtmat)],c("Tumor","Normal")),
           factor(human_dc$cell_to_patient[colnames(adtmat)])))
marker_means <- lapply(s,function(x){rowMeans(adtmat[,x],na.rm=T)})
marker_means <- do.call(cbind,marker_means)
marker_means[is.nan(marker_means)] <- NA

x <- array(NA,ncol(marker_means))
x[grep("mregDC.N",colnames(marker_means))] <- 2
x[grep("mregDC.T",colnames(marker_means))] <- 1
x[grep("DC1.N",colnames(marker_means))] <- 4
x[grep("DC1.T",colnames(marker_means))] <- 3
x[grep("DC2.N",colnames(marker_means))] <- 6
x[grep("DC2.T",colnames(marker_means))] <- 5

pch <- array(1,ncol(marker_means))
pch[grep("Tumor",colnames(marker_means))] <- 16

col <- array(NA,ncol(marker_means))
col[x%in%c(1,2)] <- rgb(t(col2rgb(7)*.7/255))
col[x%in%c(3,4)] <- rgb(t(col2rgb(4)*.7/255))
col[x%in%c(5,6)] <- rgb(t(col2rgb(5)*.7/255))

group <- substr(colnames(marker_means),1,nchar(colnames(marker_means))-4)
s_group <- split(colnames(marker_means),factor(group,c("mregDC.Tumor","mregDC.Normal","DC1.Tumor","DC1.Normal","DC2.Tumor","DC2.Normal")))
group_means <- lapply(s_group,function(x){rowMeans(log2(1+marker_means[,x]),na.rm=T)})
group_means <- do.call(cbind,group_means)

x0 <- 1:6-0.25; x1=1:6+0.25

markers <- c("CD40","CD141","CD103","CD1c","HLADR","XCR1","CD273PDL2","CD274PDL1","CD86")
names(markers) <- markers
names(markers)[markers=="CD273PDL2"] <- "PD-L2"
names(markers)[markers=="CD274PDL1"] <- "PD-L1"

png(file.path(figure_dir,"fig_s4a.png"),height=3.65,width=2.73,units="in",res=300,pointsize=6)
layout(matrix(1:9,nrow=3,byrow=T))
par(mar=c(5.1,3,4.1,2.1),oma=c(2,0,0,0))
for(marker_iter in markers){
  plot(x,log2(1+marker_means[marker_iter,]),pch=pch,col=col,
       xaxt="n",bty="L",xlim=c(0.5,6.5),ylab="",xlab="")
  segments(x0=x0,x1=x1,y0=group_means[marker_iter,],y1=group_means[marker_iter,])
  mtext(at=1:6,
        c(expression("mregDC"["Tumor"]),
          expression("mregDC"["nLung"]),
          expression("DC1"["Tumor"]),
          expression("DC1"["nLung"]),
          expression("DC2"["Tumor"]),
          expression("DC2"["nLung"])),side=1,las=2,cex=0.75,line=0.25)
  mtext(names(markers)[markers==marker_iter],cex=1.5)
}
dev.off()
}