

plot_1h = function(){

#mreg_gated <- read.csv("mDC_gating_naive_190424.csv",r=1)

if(!"mregDC1v2_DE_40k.csv"%in%list.files(file.path(wd,"DE_results/"))){
    
#mreg_gated generated from plotting fig 1f
    if(!exists("mreg_gated")){
      error("need to run fig1f")
    }
mask_bg <- names(mreg_gated)[mreg_gated=="DC2"]
mask_fg <- names(mreg_gated)[mreg_gated=="DC1"]

# mask_bg <- sample(mask_bg,pmin(length(mask_bg),length(mask_fg)))
# mask_fg <- sample(mask_fg,pmin(length(mask_bg),length(mask_fg)))
  DE <- DE_between_two_sets(ldm,mask_bg = mask_bg,
                            mask_fg=mask_fg,nmin_umi_thresh = 0,nmin_cells_with_min_umi = 5,nchunks = 20,n_per_chunk = 2000,reg=1e-7)
  write.csv(DE,file.path(wd,"DE_results/mregDC1v2_DE_40k.csv"))
}else{
  DE <- read.csv(file.path(wd,"DE_results/mregDC1v2_DE_40k.csv"),r=1,h=1,stringsAsFactors = F)
}

mask_bg <- names(mreg_gated)[mreg_gated=="DC2"]
mask_fg <- names(mreg_gated)[mreg_gated=="DC1"]
  
genes <- strsplit("Il12b,Ccl17,Irf8,Fcer1g,Cst3,Xcr1,Sirpa,Cadm1,Ccl22,Tap1,Il4ra",",")[[1]]
sig_genes_up <- rownames(DE)[DE$adj.p.value<0.15 & DE$log2_FC>0]
sig_genes_down <- rownames(DE)[DE$adj.p.value<0.15 & DE$log2_FC<0]
#m <- log10(1e-6+rowMeans(cbind(DE$freq_mDC2,DE$freq_mDC1)))
rs <- rowSums(ldm$dataset$umitab[,ldm$dataset$cell_to_sample=="Naive" & ldm$dataset$cell_to_cluster==annots$node[annots$parent=="mregDC"]])
m <- log10(5e-6+rs/sum(rs))

reg <- 5e-6
rsDC1 <- rowSums(ldm$dataset$umitab[,mask_fg])
rsDC2 <- rowSums(ldm$dataset$umitab[,mask_bg])
m_DC1 <- rsDC1/sum(rsDC1)
m_DC2 <- rsDC2/sum(rsDC2)
l2fc <- log2((m_DC1+reg)/(m_DC2+reg))

png(file.path(figure_dir,"fig1h.png"),width=2.645,height=4.2,units="in",res=300)
par(mfrow=c(2,1),oma=c(0,0,0,0),mar=c(2,2,2,.5))

for(iter in 1:2){
  plot(m,l2fc,
       bty="n",
       pch=20,cex=.5,col=rgb(.6,.6,.6,.2),bty="n",
       xlab="",
       ylab="",cex.axis=.5,tcl=.15,cex.lab=.5,mgp=c(1,.1,0))
title(xlab=expression("Log"[10]*"(Expression)"), line=.75, cex.lab=.7, family="sans")
title(ylab=expression("Log"[2]*"(mregDC1/mregDC2)"),line=.75, cex.lab=.7, family="sans")
title(main="mregDC",line=1, cex.main=.7, family="sans")
  
  points(m[c(sig_genes_up,sig_genes_down)],l2fc[c(sig_genes_up,sig_genes_down)],pch=20,cex=.5,col=rgb(0,.7,0,alpha=0.2))
  
  abline(h=0,lty=2,col="grey")
  
  points(m[genes],l2fc[genes],pch=20,cex=1,col=rgb(1,.5,0,alpha=1))
}

text(x=m[genes],y=l2fc[genes], adj=c(1,0),labels=genes,cex=.3,font=3,ps=5)

dev.off()
}