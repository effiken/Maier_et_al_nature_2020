
library(scales)
library(viridis)
library(Matrix)

plot_4d <- function(){

  s <- split(names(human_dc$cell_to_annot),human_dc$cell_to_annot)
  exprs <- lapply(s,function(x){rs <- rowSums(human_dc$filtered_umitab[,x]); return(rs/sum(rs))})
  
  reg <- 1e-6
  dc2_vs_dc1 <- log2((reg+exprs$DC2)/(reg+exprs$DC1))
  mreg_vs_dc1 <- log2((reg+exprs$mregDC)/(reg+exprs$DC1))
  mreg_vs_dc2 <- log2((reg+exprs$mregDC)/(reg+exprs$DC2))
  
  gene_sets <- list()
  gene_sets$dc2 <- names(dc2_vs_dc1)[dc2_vs_dc1 > 1 & exprs$DC2>1e-5]
  gene_sets$dc1 <- names(dc2_vs_dc1)[dc2_vs_dc1 < -1 & exprs$DC1>1e-5]
  gene_sets$mreg <- names(dc2_vs_dc1)[mreg_vs_dc1 > 1.5 & mreg_vs_dc2 > 1.5 & exprs$mregDC>1e-5]
  
  numi <- colSums(human_dc$filtered_umitab)
  
  gene_scores <- lapply(gene_sets,function(x){log10(1e-4+colSums(human_dc$filtered_umitab[x,])/numi)})
  
  annot2col <- rgb(t(col2rgb(c(7,4,5))*.7/255)); names(annot2col) <- c("mregDC","DC1","DC2")
  
  colvec <- log10(1e-4+human_dc$filtered_umitab["IL12B",]/numi)
  colvec <- round((colvec-min(colvec))/diff(range(colvec))*49+1)
  colvec <- rev(magma(50))[colvec]
  
  xlim=c(-0.9,1.4); ylim=c(-2.8,-0.9)
  
  png(file.path(figure_dir,"fig4d.png"),height=1.4,width=4.49,
      units="in",res=1000,pointsize=5)
  layout(matrix(1:2,nrow=1))
  par(bg=NA,bty="n",pin=c(1.1,0.9))
  plot(gene_scores$dc2-gene_scores$dc1,gene_scores$mreg,pch=16,col=alpha(annot2col[human_dc$cell_to_annot],0.6),
       xlim=xlim,ylim=ylim,xlab="",ylab="")
  
  plot(gene_scores$dc2-gene_scores$dc1,gene_scores$mreg,pch=16,col=alpha(colvec,0.6),
       xlim=xlim,ylim=ylim,xlab="",ylab="")
  dev.off()
  
}




