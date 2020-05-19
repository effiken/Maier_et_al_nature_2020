
plot_s2i <- function(){


if(!"DE_mreg_KPvsNaive_10k.csv"%in%list.files(file.path(wd,"DE_results/"))){
  mask_fg <- colnames(ldm$dataset$umitab)[ldm$dataset$cell_to_sample=="KP" & ldm$dataset$cell_to_cluster%in%c("5")]
  mask_bg <- colnames(ldm$dataset$umitab)[ldm$dataset$cell_to_sample=="Naive" & ldm$dataset$cell_to_cluster%in%c("5")]
  DE <- DE_between_two_sets(ldm,mask_bg,mask_fg,1,3,1e-6,10,1000)
  write.csv(DE,file.path(wd,"DE_results/DE_mreg_KPvsNaive_10k.csv"))
}else{
  DE <- read.csv(file.path(wd,"DE_results/DE_mreg_KPvsNaive_10k.csv"),r=1,h=1,stringsAsFactors = F)
}

genes <- strsplit("Ciita,Pdcd1lg2,Il4ra,Cd274,Il12b",",")[[1]]

sig_genes_up <- rownames(DE)[DE$adj.p.value<0.15 & DE$log2_FC>0]
sig_genes_down <- rownames(DE)[DE$adj.p.value<0.15 & DE$log2_FC<0]
rs <- rowSums(ldm$dataset$umitab[,ldm$dataset$cell_to_sample%in%c("Naive","KP") & ldm$dataset$cell_to_cluster==annots$node[annots$parent%in%c("mregDC")]])
m <- rs/sum(rs)

m <- log10(1e-6+m)

rs_mdc <- rowSums(ldm$dataset$umitab[,ldm$dataset$cell_to_cluster=="5" & ldm$dataset$cell_to_sample%in%c("KP")])
rs_rdc <- rowSums(ldm$dataset$umitab[,ldm$dataset$cell_to_cluster%in%c("5") & ldm$dataset$cell_to_sample%in%c("Naive")])
rs_mdc <- rs_mdc/sum(rs_mdc)
rs_rdc <- rs_rdc/sum(rs_rdc)
l2fc <- log2((1e-6+rs_mdc)/(1e-6+rs_rdc))
l2fc <- l2fc[names(m)]

png(file.path(figure_dir,"fig_S2i.png"),width=2.645,height=4.2,units="in",res=300,pointsize=8)
par(mfrow=c(2,1))

for(iter in 1:2){
  plot(m,l2fc,
       bty="n",
       pch=20,cex=.5,col=rgb(.6,.6,.6,.2),bty="n",
       xlab="",
       ylab="",tcl=.15,mgp=c(1,.1,0))
title(xlab=expression("Log"[10]*"(Expression)"), line=1.5, family="sans")
#title(ylab=expression("Log"[2]*"(mregDC"["Tumor"]*"/mregDC"["Naive"]*")"),line=1.5, family="sans")
mtext(side=2,expression("Log"[2]*"(mregDC"["Tumor"]*"/"),line=2.5)
mtext(side=2,expression("mregDC"["Naive"]*")"),line=1.5)
title(main="Tumor vs. Naive mregDC",line=1, family="sans")
  
  points(m[c(sig_genes_up,sig_genes_down)],l2fc[c(sig_genes_up,sig_genes_down)],pch=20,cex=.5,col=rgb(0,.7,0,alpha=0.2))
  
  abline(h=0,lty=2,col="grey")
  
  points(m[genes],DE[genes,]$log2_FC,pch=20,cex=1,col=rgb(1,.5,0,alpha=1))
}

text(x=m[genes],y=DE[genes,]$log2_FC, adj=c(1,0),labels=genes,cex=.3,font=3,ps=5)

dev.off()
}