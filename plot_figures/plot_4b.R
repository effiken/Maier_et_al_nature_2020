
plot_4b <- function(){
s <- split(colnames(human_dc$filtered_umitab),human_dc$cell_to_annot)
exprs <- lapply(s,function(x){
  rs <- rowSums(human_dc$filtered_umitab[,x]);
  return(rs/sum(rs))
})
exprs <- do.call(cbind,exprs)
exprs <- exprs[,c("mregDC","DC1","DC2")]



pathways <- list()
pathways$'Maturation' <- rev(c("Cd40","Cd80","Cd86","Relb","Cd83"))
pathways$'Regulatory' <- rev(c("Cd274","Pdcd1lg2","Cd200","Fas","Aldh1a2","Socs1","Socs2"))
pathways$'TLRs & adaptors' <- strsplit("Tlr1,Tlr2,Tlr3,Tlr4,Tlr5,Tlr6,Tlr7,Tlr8,Tlr9,Mavs,Myd88",",")[[1]]
pathways$'Migration' <- rev(c("Ccr7","Myo1g","Cxcl16","Adam8","Icam1","Fscn1","Marcks","Marcksl1"))
pathways$'Th2 response' <- strsplit("Bcl2l1,Stat6,Tnfrsf4,Ccl22,Ccl17,Il4i1,Il4r",",")[[1]]

pathways <- lapply(pathways,toupper)

pdf(file.path(figure_dir,"fig4b.pdf"),width=7.07,height=1.35)
#png("fig4b.png",res=1000,width=7.07,height=1.35,units="in")
layout(matrix(c(1:length(pathways),array(length(pathways)+1,length(pathways))),byrow=T,ncol=length(pathways)),heights=c(4,1))
par(oma=c(0,0,0,0),mar=c(2,2,1.5,2))
for(iter in 1:length(pathways)){
  genes <- pathways[[iter]]
  mat <- exprs
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
  mtext(c("mregDC","DC1","DC2"),side=1,at=seq(0,1,1/2),las=1,line=.1,cex=.5)
  mtext(names(pathways)[iter],font=2,cex=.7)
}

par(mar=c(1,20,.2,20))
image(as.matrix(1:50),col=bluered(50),xaxt="n",yaxt="n")
mtext(c(-2,2),at=c(0,1),side=1,cex=.5)
mtext("Log2(Expression / Row Mean)",cex=.7)
dev.off()
}