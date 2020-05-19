
plot_s3e <- function(){
  

DC_expr_mat <- matrix(NA,nrow=nrow(ldm$dataset$umitab),ncol=ncol(ldm$model$models),
                      dimnames=list(rownames(ldm$dataset$umitab),colnames(ldm$model$models)))
s <- split(colnames(ldm$dataset$umitab),ldm$dataset$cell_to_cluster)
for(iter in names(s)){
  rs <- rowSums(ldm$dataset$umitab[,s[[iter]]])
  DC_expr_mat[,iter] <- rs/sum(rs)
}

mat <- DC_expr_mat[,annots$node[match(c("mregDC","DC1","DC2"),annots$parent)]]
colnames(mat) <- c("mregDC","DC1","DC2")

mat <- log10(1e-6+mat[rev(c("Axl","Tyro3","Mertk")),])

s <- seq(0,1,1/50)
png(file.path(figure_dir,"figS3e.png"),height=1.88,width = 2.35,units="in",res=1000,pointsize=8)
par(oma=c(0,0,0,1))
layout(matrix(c(1,2),nrow=1),widths=c(5,1))
#par(mar=c(5,4,10,1))
image(t(mat),col=rev(rgb(s,s,s)),xaxt="n",yaxt="n")
box()
mtext(rownames(mat),at=seq(0,1,1/(nrow(mat)-1)),side=2,las=2,font=3,line=.5)
mtext(colnames(mat),at=seq(0,1,1/(ncol(mat)-1)),side=1,line=.5)
mtext("TAM receptors",font=2,line=1,cex=1.2)

par(mar=c(5.1,.5,4.1,1.5))

image(as.matrix(t(rev(s))),col=rgb(s,s,s),xaxt="n",yaxt="n"); box();
mtext(c(min(mat),round(max(mat),1)),at=c(0,1),side=4,las=2,line=.5)
mtext("Log10(1e-6+\nExpression)",line=.5)
dev.off()
}