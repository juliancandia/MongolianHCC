rm(list=ls())

PROJECT_DIR = "~/ACTIVE/LCS/Mongolia_HCC_Pipeline" # replace this line with your local path

setwd(file.path(PROJECT_DIR,"RESULTS","cluster"))

maxK=8
mad_T_thr = c(1,1.25,1.5,1.75,2,2.25,2.5,2.75,3)
n_thres = length(mad_T_thr)

ratio_mean = matrix(rep(NA,n_thres*maxK),ncol=maxK)
ratio_wmean = matrix(rep(NA,n_thres*maxK),ncol=maxK)
for (i_thres in 1:n_thres) {
    thres = mad_T_thr[i_thres]
    infile = file.path(paste0("T_mad",thres),"patient_sample_metadata_w_clustering.txt")
    data = as.matrix(read.table(infile,header=T,check.names=F,stringsAsFactors=F,sep="\t"))
    data = data[!is.na(data[,"RNASeq_ID"]),] # only patients with RNASeq data
    for (i_clust in 2:maxK) {
        ratio = NULL
        size = NULL
        for (k in 1:i_clust) {
            patient_sel = as.numeric(data[,paste0("K",i_clust,"_class")])==k
            ratio = c(ratio,mean(as.numeric(data[patient_sel,paste0("K",i_clust,"_ratio")])))
            size = c(size,sum(patient_sel))
        }
        ratio_mean[i_thres,i_clust] = mean(ratio)
        ratio_wmean[i_thres,i_clust] = sum(ratio*size)/sum(size)
    }
}

ratio_mean[,1] = mad_T_thr
output = rbind(c("mad_T_thr",paste0("K",2:maxK)),ratio_mean)
dir.create("FINAL")
outfile = "FINAL/InOutRatio.txt"
write(t(output),ncol=ncol(output),file=outfile,sep="\t")

outfile = "FINAL/InOutRatio.pdf"
pdf(outfile,width=7,height=5)
x = mad_T_thr
xrange = range(x)
x_legend = xrange[1] + (xrange[2]-xrange[1])*1.035
xrange[2] = xrange[1] + (xrange[2]-xrange[1])*1.15
yrange = range(as.numeric(ratio_mean[,-1]),na.rm=T)
plot(xrange,yrange,type="n",xlab="mad_T_thr",ylab="in/(in+out) ratio",cex.lab=1)
title(main="comparison of clustering solutions",cex.main=0.95)
color = c("black","red","green","blue","cyan","magenta","orange") #"brown","darkgreen"
symbol = c(8,6,2,20,0,1,5) # other useful symbols: 17,15,3,4
for (i_clust in 2:maxK) {
    points(x,as.numeric(ratio_mean[,i_clust]),type="b",lty=2,lwd=2,pch=symbol[i_clust-1],col=color[i_clust-1])
}
y_legend = yrange[1] + (yrange[2]-yrange[1])*1
legend_txt = paste0("Cl ",2:maxK)
legend_pch = symbol
legend_lwd = 2
legend_lty = 2
legend_col = color
legend_text_size = 0.8
legend_symbol_size = 0.8
legend(x_legend,y_legend,legend_txt,pch=legend_pch,lwd=legend_lwd,col=legend_col,lty=legend_lty,
cex=legend_text_size,pt.cex=legend_symbol_size)
dev.off()
