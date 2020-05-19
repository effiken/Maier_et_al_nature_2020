



plot_2a <- function(){
  

s <- split(colnames(ldm$dataset$umitab),
           list(ldm$dataset$cell_to_cluster,ldm$dataset$cell_to_sample))
s <- s[c("5.Naive","5.Naive_CCR7KO")]
exprs <- lapply(s,function(x){rs <- rowSums(ldm$dataset$umitab[,x]); return(rs/sum(rs))})
exprs <- do.call(cbind,exprs)

reg <- 1e-6
png(file.path(figure_dir,"Fig2a.png"),height=1.91,width=1.91,units="in",res=600,pointsize = 6)
#par(mar=c(2,2,2,2))
plot(log10(reg+exprs[,"5.Naive"]),log10(reg+exprs[,"5.Naive_CCR7KO"]),bty="L",pch=".",xlab="",ylab="",col=alpha("black",.3))
title(xlab=expression("Log"[10]*"(mregDC"["WT"]*" Expression)"))
title(ylab=expression("Log"[10]*"(mregDC"[italic("Ccr7"^"-/-")]*" Expression)"))
title(main=expression("mregDC: WT vs. "*italic("Ccr7"^"-/-")))
abline(0,1,col="grey",lwd=.5)
abline(-log10(2),1,col="red",lwd=.5)
abline(log10(2),1,col="red",lwd=.5)
legend(x=-6,y=-1,legend=expression("Log"[2]*"FC = \u00B1 1"),col="red",lwd=.5,bty="n")
dev.off()

}
