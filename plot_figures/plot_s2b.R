
plot_s2b <- function(){

if(!"figS2b_DE_WT_10k.csv"%in%list.files(file.path(wd,"DE_results/"))){
  mask_fg <- colnames(ldm$dataset$umitab)[ldm$dataset$cell_to_sample=="Naive" & ldm$dataset$cell_to_cluster%in%c("5")]
  mask_bg <- colnames(ldm$dataset$umitab)[ldm$dataset$cell_to_sample=="Naive" & ldm$dataset$cell_to_cluster%in%c("3","4")]
  DE_WT <- DE_between_two_sets(ldm,mask_bg,mask_fg,1,3,1e-6,10,1000)
  write.csv(DE_WT,file.path(wd,"DE_results/figS2b_DE_WT_10k.csv"))
}else{
  DE_WT <- read.csv(file.path(wd,"DE_results/figS2b_DE_WT_10k.csv"),r=1,h=1,stringsAsFactors=F)
}

if(!"figS2b_DE_KO_10k.csv"%in%list.files(file.path(wd,"DE_results/"))){
  mask_fg <- colnames(ldm$dataset$umitab)[ldm$dataset$cell_to_sample=="Naive_CCR7KO" & ldm$dataset$cell_to_cluster%in%c("5")]
  mask_bg <- colnames(ldm$dataset$umitab)[ldm$dataset$cell_to_sample=="Naive_CCR7KO" & ldm$dataset$cell_to_cluster%in%c("3","4")]
  DE_KO <- DE_between_two_sets(ldm,mask_bg,mask_fg,1,3,1e-6,10,1000)
  write.csv(DE_KO,file.path(wd,"DE_results/figS2b_DE_KO_10k.csv"))
}else{
  DE_KO <- read.csv(file.path(wd,"DE_results/figS2b_DE_KO_10k.csv"),r=1,h=1,stringsAsFactors=F)
}


rest_DC_WT <- rowSums(ldm$dataset$umitab[,ldm$dataset$cell_to_sample=="Naive" & ldm$dataset$cell_to_cluster%in%c("3","4")])
rest_DC_WT <- rest_DC_WT/sum(rest_DC_WT)
rest_DC_KO <- rowSums(ldm$dataset$umitab[,ldm$dataset$cell_to_sample=="Naive_CCR7KO" & ldm$dataset$cell_to_cluster%in%c("3","4")])
rest_DC_KO <- rest_DC_KO/sum(rest_DC_KO)
mgDC_WT <- rowSums(ldm$dataset$umitab[,ldm$dataset$cell_to_sample=="Naive" & ldm$dataset$cell_to_cluster%in%c("5")])
mgDC_WT <- mgDC_WT/sum(mgDC_WT)
mgDC_KO <- rowSums(ldm$dataset$umitab[,ldm$dataset$cell_to_sample=="Naive_CCR7KO" & ldm$dataset$cell_to_cluster%in%c("5")])
mgDC_KO <- mgDC_KO/sum(mgDC_KO)

include <- rest_DC_WT > 1e-6 & rest_DC_KO > 1e-6

reg <- 1e-6

l2fc_WT <- log2((reg+mgDC_WT)/(reg+rest_DC_WT))
l2fc_KO <- log2((reg+mgDC_KO)/(reg+rest_DC_KO))

l2fc_WT <- l2fc_WT[include]
l2fc_KO <- l2fc_KO[include]

sig_wt <- rownames(DE_WT)[DE_WT$adj.p.value<0.01]
sig_ko <- rownames(DE_KO)[DE_KO$adj.p.value<0.01]
sig_both <- intersect(sig_wt,sig_ko)
#sig_both <- sig_both[DE_WT[sig_both,]$log2_FC*DE_KO[sig_both,]$log2_FC>0]

png(file.path(figure_dir,"figS2b.png"),height=2,width=2,units="in",res=1000,pointsize=6)
par(mar=c(5.1,4.1,4.1,6))
plot(l2fc_WT,l2fc_KO,bty="n",xlab="",ylab="",pch=16,col=rgb(0,0,0,.4),cex=.2)
points(l2fc_WT[sig_wt],l2fc_KO[sig_wt],col="pink",pch=16,cex=.2)
points(l2fc_WT[sig_ko],l2fc_KO[sig_ko],col="green",pch=16,cex=.2)
points(l2fc_WT[sig_both],l2fc_KO[sig_both],col="purple",pch=16,cex=.2)
title(xlab=expression("Log"[2]*"(mregDC/resting DC)"["WT"]),family="sans")
title(ylab=expression("Log"[2]*"(mregDC/resting DC)"[italic("Ccr7"^"-/-")]),family="sans")
title(main=expression("mregDC vs. resting DC in WT & "*italic("Ccr7"^"-/-")))
abline(h=0,col="grey",lwd=.5)
abline(v=0,col="grey",lwd=.5)
abline(0,1,col="grey",lty=2,lwd=.5)
legend(x=1,y=-2,legend=c(
  expression("Adj. P"["WT"]*" < 0.01"),
  expression("Adj. P"["Ccr7"^"-/-"]*" < 0.01"),
  expression("Adj. P"["Both"]*" < 0.01")),col=c("pink","green","purple"),bty="n",xpd="",pch=16)
dev.off()
}