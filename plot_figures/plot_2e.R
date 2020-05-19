# Fig 2e:
# Compare sorted bulk lung GFP+ vs. GFP- DC1 with mature DC & DC1 in single cell data

library(scales)

plot_2e <- function(){


sort_expt1 <- read.csv(file.path(wd,"additional_input_data/voom_exprs.csv"),r=1,h=1,stringsAsFactors = F)

if(!"DE_mregDC_vs_DC1_20k.csv"%in%list.files(file.path(wd,"DE_results/"))){
  mask_fg <- colnames(ldm$dataset$umitab)[ldm$dataset$cell_to_sample%in%c("Naive","KP") & ldm$dataset$cell_to_cluster%in%c("5")]
  mask_bg <- colnames(ldm$dataset$umitab)[ldm$dataset$cell_to_sample%in%c("Naive","KP") & ldm$dataset$cell_to_cluster%in%c("4")]
  DE <- DE_between_two_sets(ldm,mask_bg,mask_fg,1,3,1e-6,20,1000)
  write.csv(DE,file.path(wd,"DE_results/DE_mregDC_vs_DC1_20k.csv"))
}else{
  DE <- read.csv(file.path(wd,"DE_results/DE_mregDC_vs_DC1_20k.csv"),r=1,h=1,stringsAsFactors = F)
}

mdc_avg <- rowSums(ldm$dataset$umitab[,ldm$dataset$cell_to_sample%in%c("KP","naive") & ldm$dataset$cell_to_cluster=="5"])
mdc_avg <- mdc_avg/sum(mdc_avg)
dc1_avg <- rowSums(ldm$dataset$umitab[,ldm$dataset$cell_to_sample%in%c("KP","naive") & ldm$dataset$cell_to_cluster=="4"])
dc1_avg <- dc1_avg/sum(dc1_avg)

l2fc_mdc_dc1 <- log2((1e-6+mdc_avg)/(1e-6+dc1_avg))

l2fc_mdc_dc1 <- l2fc_mdc_dc1[dc1_avg>1e-6]


sort_expt1 <- 2^sort_expt1; sort_expt1 <- t(t(sort_expt1)/rowSums(t(sort_expt1)))

l2fc_gfp <- log2((1e-6+rowMeans(sort_expt1[,1:3]))/(1e-6+rowMeans(sort_expt1[,4:5])))

overlap <- intersect(names(l2fc_mdc_dc1),names(l2fc_gfp))

l2fc_mdc_dc1 <- l2fc_mdc_dc1[overlap]
l2fc_gfp <- l2fc_gfp[overlap]

#genes <- read.csv("/users/andrew leader/Dropbox/merad_lab_analyses/AL/AL_C126_mDC_mouse_DEGs/DE_mDC_naive.csv",r=1,h=1,stringsAsFactors = F)
genes <- rownames(DE)[DE$log2_FC>2 & DE$adj.p.value<0.01]

genes <- intersect(genes,overlap)

#genes_select <- c("Fscn1","Nudt17","Il12b","Ccr7","Socs2","Tnfrsf4","Cd274","Itgb8","Vsig10","Cd40")
genes_select <- c("Fscn1","Il12b","Ccr7","Cd274","Cd40","Ccl22","Fas","Il4ra","Ccl17","Pdcd1lg2")


png(file.path(figure_dir,"fig2e.png"),height=5,width=2.2,units="in",pointsize = 7,res=1000)
par(oma=c(1,0,0,0),mfrow=c(2,1))

for(iter in 1:2){
  plot(l2fc_mdc_dc1,l2fc_gfp,pch=16,col=rgb(0,0,0,.05),bty="n",
       xlab="",
       ylab="",cex=.5,ylim=c(-3,4))
  title(xlab=expression("Log"[2]*"(mregDC/DC1)"),line=2)
  title(ylab=expression("Log"[2]*"(GFP"^"+"*" DC1/GFP"^"-"*" DC1)"),line=2)
  title(main=expression("GFP"^"+"*" vs. mregDC signatures"))
  abline(h=0,col=alpha("grey",.4))
  abline(v=0,col=alpha("grey",.4))
  points(l2fc_mdc_dc1[genes],l2fc_gfp[genes],pch=16,col=alpha(rgb(t(col2rgb("yellow"))*.7,max=255),.3),cex=.5)
  points(l2fc_mdc_dc1[genes_select],l2fc_gfp[genes_select],pch=16,col="red",cex=.5)
#  points(l2fc_mdc_dc1[c("Apoe","C1qb","C1qc")],l2fc_gfp[c("Apoe","C1qb","C1qc")],pch=16,col=alpha("black",.2),cex=.5)
  
  legend(x=-4,y=-1.5,legend="FDR < 0.01 & Log2FC\n> 2 in mregDC vs. DC1",pch=16,col=rgb(t(col2rgb("yellow"))*.7,max=255),bty="n",xpd="")
}

text(labels=genes_select,x=l2fc_mdc_dc1[genes_select],y=l2fc_gfp[genes_select],cex=.5)
#text(labels=c("Apoe","C1qb","C1qc"),x=l2fc_mdc_dc1[c("Apoe","C1qb","C1qc")],y=l2fc_gfp[c("Apoe","C1qb","C1qc")],cex=.5)

dev.off()

}
