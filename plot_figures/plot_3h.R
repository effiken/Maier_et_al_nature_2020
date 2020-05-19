
plot_3h <- function(){
  


pathways <- list()
pathways$'Th2 response' <- c("Bcl2l1","Stat6","Tnfrsf4","Ccl22","Ccl17","Il4i1","Il4ra")


s <- split(colnames(ldm$dataset$umitab)[ldm$dataset$cell_to_sample%in%c("KP","Naive")],ldm$dataset$cell_to_cluster[ldm$dataset$cell_to_sample%in%c("KP","Naive")])

mat_avg <- matrix(NA,nrow=nrow(ldm$dataset$umitab),ncol=ncol(ldm$model$models),dimnames=list(rownames(ldm$dataset$umitab),colnames(ldm$model$models)))
for(iter in names(s)){
  rs <- rowSums(ldm$dataset$umitab[,s[[iter]]])
  mat_avg[,iter] <- rs/sum(rs)
}
colnames(mat_avg) <- annots[match(colnames(mat_avg),annots$node),"parent"]

mat_avg <- mat_avg[,c("mregDC","DC1","DC2")]

pdf(file.path(figure_dir,"fig3h.pdf"),width=1.36,height=1.01,pointsize=9)
#png("fig3h.png",res=300,width=1.36,height=1.01,units="in",pointsize=9)
par(oma=c(0,0,0,0),mar=c(2,2,1.5,2))
for(iter in 1:length(pathways)){
  genes <- pathways[[iter]]
  mat <- mat_avg
  mat <- (1e-6+mat)/(1e-6+rowMeans(mat))
  mat <- log2(mat)
  mat[mat > 2] <- 2
  mat[mat < -2] <- -2
  mat <- round((mat+2)/4*50)
  col <- bluered(50)
  col <- col[min(mat):max(mat)]
  image(t(mat[genes,]),col=col,xaxt="n",yaxt="n")
  box()
  mtext(genes,at=seq(0,1,1/(length(genes)-1)),side=2,las=2,line=.25,cex=.4,font=3)
  mtext(c("mregDC","DC1","DC2"),side=1,at=seq(0,1,1/2),las=1,line=.1,cex=.6)
  mtext(names(pathways)[iter],font=2,cex=.7)
}

dev.off()
}