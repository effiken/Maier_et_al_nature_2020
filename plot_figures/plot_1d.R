
plot_1d <- function(){
  

rs <- rowSums(ldm$dataset$umitab[,ldm$dataset$cell_to_cluster==annots$node[annots$parent=="mregDC"]])
mdc_avg <- rs/sum(rs)

rs <- rowSums(ldm$dataset$umitab[,ldm$dataset$cell_to_cluster==annots$node[annots$parent%in%c("DC1","DC2")]])
rest_avg <- rs/sum(rs)

rs <- rowSums(ldm$dataset$umitab[,ldm$dataset$cell_to_cluster==annots$node[annots$parent%in%c("mregDC","DC1","DC2")]])
total_avg <- rs/sum(rs)

reg <- 1e-6
l2fc <- log2((mdc_avg+reg)/(rest_avg+reg))

miller <- read.csv(file.path(wd,"additional_input_data/miller.csv"))
miller_up <- miller$ï..miller_up
miller_down <- miller$miller_down
miller_up <- intersect(miller_up,names(l2fc))
miller_down <- intersect(miller_down,names(l2fc)) 

png(file.path(figure_dir,"fig1d.png"),height=1.96,width=2.11,units="in",res=1000,pointsize=6)
par(oma=c(0,0,0,1))
plot(log10(total_avg+reg),l2fc,
     bty="n",
     pch=20,
     cex=.5,
     col=rgb(.6,.6,.6,.2),
     xlab="",ylab="",
     ylim=c(-8,8),
     tcl=.15,mgp=c(1,.25,0))
title(xlab=expression("Log"[10]*"(Expression)"), line=1.5, cex.lab=1.2, family="sans")
title(ylab=expression("Log"[2]*"(mregDC/resting DC)"),line=1.5, cex.lab=1.2, family="sans")
title(main="Lung DC",line=1, cex.lab=1.4, family="sans")

abline(h=0,lty=2,col="grey")

points(log10(total_avg[miller_up]+reg),l2fc[miller_up],pch=20,cex=.5,col="red")
#points(log10(total_avg[miller_down]+reg),l2fc[miller_down],pch=20,cex=.5,col="blue")
#legend(x=-3,y=8,xpd="",legend=c("migDC genes","tissue DC genes"),pch=20,col=c("red","blue"),bty="n")
legend(x=-3,y=8,xpd="",legend="migDC genes",pch=20,col="red",bty="n")

dev.off()
}
