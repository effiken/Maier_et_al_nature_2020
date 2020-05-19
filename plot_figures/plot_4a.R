
plot_4a <- function(){


clust2annot <- c("mregDC","DC1","DC2")
names(clust2annot) <- c("41","45","14")

cells <- colnames(human_dc$filtered_umitab)[colnames(human_dc$filtered_umitab)%in%colnames(human_dc$filtered_ds)]
cells <- split(cells,human_dc$cell_to_annot[cells])
cells <- unlist(lapply(cells,function(x){sample(x,min(unlist(lapply(cells,length))))}))
cell_to_tissue <- human_dc$cell_to_tissue[cells]
names(cell_to_tissue) <- cells
cells <- cells[order(factor(human_dc$cell_to_annot[cells],c("mregDC","DC1","DC2")),factor(cell_to_tissue[cells],c("Normal","Tumor")))]


lin_col <- rgb(t(col2rgb(c(7,4,5)))*.7,max=255)

breaks <- which(!diff(as.numeric(factor(human_dc$cell_to_annot[cells])))==0)/length(cells)


genes <- rev(strsplit("MARCKSL1,CD80,TRAF1,RELB,CCL22,BIRC2,IL12B,ENO3,LAMP3,MARCKS,FSCN1,CD274,PDCD1LG2,CD200,FAS,SOCS2,CCL17,IL4I1,CCL19,BIRC3,CCR7,IDO1,NAAA,C1orf54,DNASE1L3,CLEC9A,XCR1,IRF8,IRF4,SIRPA,CD1C,CD1E,FCER1A,FCGR2B,CLEC10A",",")[[1]])

mat <- log2(1+human_dc$filtered_ds[genes,cells])

mat[mat > 4] <- 4

s <- seq(1,0,-.02)
cell2annot <- human_dc$cell_to_annot[cells]


pdf(file.path(figure_dir,"fig_4a.pdf"),height=3.74,width=2.55,pointsize = 12,compress=FALSE)
#png("fig_4a.png",height=3.74,width=2.55,units="in",res=1000)
layout(matrix(1:3,nrow=3),heights=c(.5,10,.5))
par(mar=c(0,0,0,0),oma=c(2,5,2,2))
image(as.matrix(as.numeric(factor(cell2annot,c("mregDC","DC1","DC2")))),xaxt="n",yaxt="n",col=lin_col)
box()
abline(v=breaks,col="red")

s <- seq(1,0,-1/49)
image(t(as.matrix(mat)),col=rgb(s,s,s),xaxt="n",yaxt="n")
box()
abline(v=breaks,col="red")
mtext(genes,side=2,at=seq(0,1,1/(length(genes)-1)),las=2,line=.25,font=3,cex=.6)

image(as.matrix(as.numeric(factor(cell_to_tissue[cells],c("Normal","Tumor")))),xaxt="n",yaxt="n",col=c(rgb(255,192,203,max=255),rgb(165,75,42,max=255)))
abline(v=breaks,col="red")
box()

dev.off()
}
