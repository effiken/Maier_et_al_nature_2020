

plot_s3f=function(){
  
if("figS2b_DE_WT_10k.csv"%in%list.files(file.path(wd,"DE_results/"))){
  DE <- read.csv(file.path(wd,"DE_results/figS2b_DE_WT_10k.csv"),r=1,h=1,stringsAsFactors = F)
}else{
  error("Need to run fig S2b first")
}

genes <- c("Tnfrsf4","Il4i1","Ccl22","Ccl17","Il4ra","Bcl2l1","Stat6")

sig_genes_up <- rownames(DE)[DE$adj.p.value<0.01 & DE$log2_FC>0]
sig_genes_down <- rownames(DE)[DE$adj.p.value<0.01 & DE$log2_FC<0]
rs <- rowSums(ldm$dataset$umitab[,ldm$dataset$cell_to_sample%in%c("Naive","KP") & ldm$dataset$cell_to_cluster%in%annots$node[annots$parent%in%c("mregDC","DC1","DC2")]])


mask_fg <- ldm$dataset$cell_to_cluster%in%c("5") & ldm$dataset$cell_to_sample%in%c("Naive","KP")
mask_bg <- ldm$dataset$cell_to_cluster%in%c("3","4") & ldm$dataset$cell_to_sample%in%c("Naive","KP")

rs_mdc <- rowSums(ldm$dataset$umitab[,mask_fg])
rs_rdc <- rowSums(ldm$dataset$umitab[,mask_bg])
rs_mdc <- rs_mdc/sum(rs_mdc)
rs_rdc <- rs_rdc/sum(rs_rdc)
l2fc <- log2((1e-6+rs_mdc)/(1e-6+rs_rdc))

#m <- rs/sum(rs)
m <- (rs_mdc+rs_rdc)/2
m <- log10(1e-6+m)

l2fc <- l2fc[names(m)]

sig <- c(sig_genes_up,sig_genes_down)
not_sig <- names(m)[!names(m)%in%sig]

png(file.path(figure_dir,"figS3f.png"),width=2.645,height=4.6,units="in",res=1000,pointsize=8)
par(mfrow=c(2,1))

for(iter in 1:2){
  plot(m[not_sig],l2fc[not_sig],
       bty="n",
       pch=".",col=rgb(.6,.6,.6,.2),bty="n",
       xlab="",
       ylab="",tcl=.15,mgp=c(1,.1,0),
       ylim=c(-8,8),xlim=c(-6,-1))
title(xlab=expression("Log"[10]*"(Expression)"), line=1.5, family="sans")
title(ylab=expression("Log"[2]*"(mregDC/resting DC)"),line=1.5, family="sans")

title(main="mregDC vs. resting DC",line=1, family="sans")
  
  points(m[sig],l2fc[sig],pch=".",col=rgb(0,.7,0,alpha=0.2))
  
  abline(h=0,lty=2,col="grey")
  
  points(m[genes],l2fc[genes],pch=16,cex=.8,col=rgb(1,.5,0,alpha=1))
}

text(x=m[genes],y=l2fc[genes], adj=c(1,0),labels=genes,cex=.3,font=3,ps=5)

dev.off()
}