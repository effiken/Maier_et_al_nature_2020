#Needs to have Fig 1b generated first

plot_1f <- function(){
  


library(sp)
  
clusts <- c("5","4","3")
lin_col <- rgb(t(col2rgb(c(4,5,7))*0.7/255))

cells <- names(ldm$dataset$cell_to_cluster)[
  ldm$dataset$cell_to_cluster%in%clusts &
    ldm$dataset$cell_to_sample%in%c("Naive","KP")
]


adt_mat <- adt_mat[,cells]


png(file.path(figure_dir,"fig1f.png"),height=2.16,width=4.52,units="in",res=300)
layout(matrix(1:3,nrow=1))

par(oma=c(0,0,0,0),mar=c(5,2.5,5,1),cex.axis=.5,mgp=c(1,.25,0),tcl=-.15,cex=.5)

plot(log2(1+adt_mat["CD103",]),(log2(1+adt_mat["CD11b",])),bty="n",col="white",xlab="CD103",ylab="CD11b")
mtext("Naive")
plot_mask <- ldm$dataset$cell_to_cluster[cells]==clusts[3] & ldm$dataset$cell_to_sample[cells]=="Naive"
points(log2(1+adt_mat["CD103",plot_mask]),(log2(1+adt_mat["CD11b",plot_mask])),pch=20,col=paste(lin_col[2],"2D",sep=""))
plot_mask <- ldm$dataset$cell_to_cluster[cells]==clusts[2] & ldm$dataset$cell_to_sample[cells]=="Naive"
points(log2(1+adt_mat["CD103",plot_mask]),(log2(1+adt_mat["CD11b",plot_mask])),pch=20,col=paste(lin_col[1],"2D",sep=""))
plot_mask <- ldm$dataset$cell_to_cluster[cells]==clusts[1] & ldm$dataset$cell_to_sample[cells]=="Naive"
plot_mask_naive <- plot_mask
points(log2(1+adt_mat["CD103",plot_mask]),(log2(1+adt_mat["CD11b",plot_mask])),pch=20,col=paste(lin_col[3],sep=""))

DC2_gate <- list()
DC2_gate$x <- as.numeric(strsplit("-0.327141491733754,-0.448006574886619,0.0354537577248405,1.24410458925349,4.16903960155281,9.34206516049542,10.3331588423489,10.1156016926738,9.8980445429986,8.23010639548907,5.90949679895407,3.9998284851388,2.04181413806239",",")[[1]])
DC2_gate$x <- DC2_gate$x-1
DC2_gate$y <- as.numeric(strsplit("3.72164383990655,5.84985659293229,7.97806934595803,10.7065472344526,11.5523753798859,9.72429519459453,8.30548669257737,7.32323465271934,6.36826739174625,2.93038525224313,2.65753746339368,2.63025268450873,2.73939180004851",",")[[1]])

DC1_gate <- list()
DC1_gate$x <- as.numeric(strsplit("10.2606397924572,9.70466040995402,8.23010639548907,7.5774349464636,6.9489365140687,6.87641746417698,12.8954986051896,14.2733605531323,14.7568208857438,13.5965160874763",",")[[1]])
DC1_gate$y <- as.numeric(strsplit("7.59608244156879,5.27687623634843,2.82124613670335,1.81170931796037,0.911311614757175,-0.371072992835257,-0.316503435065366,2.41197445342917,6.36826739174625,7.59608244156879",",")[[1]])
#polygon(DC1_gate)
#polygon(DC2_gate$x,DC2_gate$y)

naive_gated_type <- array("ungated",sum(plot_mask))
naive_gated_type[which(as.logical(point.in.polygon(log2(1+adt_mat["CD103",which(plot_mask)]),
                        log2(1+adt_mat["CD11b",which(plot_mask)]),
                        DC1_gate$x, DC1_gate$y)))] <- "DC1"
naive_gated_type[which(as.logical(point.in.polygon(log2(1+adt_mat["CD103",which(plot_mask)]),
                               log2(1+adt_mat["CD11b",which(plot_mask)]),
                               DC2_gate$x, DC2_gate$y)))] <- "DC2"
names(naive_gated_type) <- colnames(adt_mat)[plot_mask]
mreg_gated <<- naive_gated_type
#write.csv(naive_gated_type,"/users/andrew leader/Dropbox/Andrew&Barbara/figure_R_scripts/fig1/mDC_gating_naive_190424.csv")

plot(log2(1+adt_mat["CD103",]),(log2(1+adt_mat["CD11b",])),bty="n",col="white",xlab="CD103",ylab="CD11b")
mtext("KP")
plot_mask <- ldm$dataset$cell_to_cluster[cells]==clusts[3] & ldm$dataset$cell_to_sample[cells]=="KP"
points(log2(1+adt_mat["CD103",plot_mask]),(log2(1+adt_mat["CD11b",plot_mask])),pch=20,col=paste(lin_col[2],"2D",sep=""))
plot_mask <- ldm$dataset$cell_to_cluster[cells]==clusts[2] & ldm$dataset$cell_to_sample[cells]=="KP"
points(log2(1+adt_mat["CD103",plot_mask]),(log2(1+adt_mat["CD11b",plot_mask])),pch=20,col=paste(lin_col[1],"2D",sep=""))
plot_mask <- ldm$dataset$cell_to_cluster[cells]==clusts[1] & ldm$dataset$cell_to_sample[cells]=="KP"
plot_mask_tumor <- plot_mask
points(log2(1+adt_mat["CD103",plot_mask]),(log2(1+adt_mat["CD11b",plot_mask])),pch=20,col=paste(lin_col[3],sep=""))


#polygon(DC1_gate)
#polygon(DC2_gate$x,DC2_gate$y)



tumor_gated_type <- array("ungated",sum(plot_mask))
tumor_gated_type[which(as.logical(point.in.polygon(log2(1+adt_mat["CD103",which(plot_mask)]),
                                                   log2(1+adt_mat["CD11b",which(plot_mask)]),
                                                   DC1_gate$x, DC1_gate$y)))] <- "DC1"
tumor_gated_type[which(as.logical(point.in.polygon(log2(1+adt_mat["CD103",which(plot_mask)]),
                                                   log2(1+adt_mat["CD11b",which(plot_mask)]),
                                                   DC2_gate$x, DC2_gate$y)))] <- "DC2"
names(tumor_gated_type) <- colnames(adt_mat)[plot_mask]
write.csv(tumor_gated_type,"DE_results/mDC_gating_KP_190424.csv")

gated_type <- cbind(table(naive_gated_type),table(tumor_gated_type))
gated_type <- t(gated_type)/rowSums(t(gated_type))

par(mar=c(5,2.5,5,5))
barplot(t(gated_type),col=c(lin_col[2],lin_col[3],"grey"),names.arg=c("Naive","Tumor"),las=2,ylab="% of mDC",cex.names=1)

dev.off()
}
