library(survival)

rm(list=ls())

PROJECT_DIR = "~/ACTIVE/LCS/Mongolia_HCC_Pipeline" # replace this line with your local path

setwd(file.path(PROJECT_DIR,"RESULTS","cluster"))

mad_T_thr = c(1,1.25,1.5,1.75,2,2.25,2.5,2.75,3)
n_thres = length(mad_T_thr)
maxK=8
risk_color = c("blue","red","orange","brown","magenta","cyan","green","black")

for (i_thres in 1:n_thres) {
    thres = mad_T_thr[i_thres]
    infile = file.path(paste0("T_mad",thres),"patient_sample_metadata_w_clustering.txt")
    data = as.matrix(read.table(infile,header=T,check.names=F,stringsAsFactors=F,sep="\t"))
    select = (!is.na(data[,"RNASeq_ID"]))&(!is.na(data[,"survival.status"]))&(!is.na(data[,"survival.time"]))
    data = data[select,]
    # Kaplan-Meier plots and log-rank test
    outfile = file.path(paste0("T_mad",thres),"KaplanMeier.pdf")
    pdf(outfile,width=7,height=5)
    for (i_clust in 2:maxK) {
        class_surv = data.frame(as.numeric(data[,"survival.time"]),
        as.numeric(data[,"survival.status"]),data[,paste0("K",i_clust,"_class")])
        colnames(class_surv) = c("surv_time","surv_st","class")
        fit <- survfit(Surv(surv_time,surv_st)~class, data=class_surv)
        logrank_pval = 1-pchisq(survdiff(Surv(surv_time,surv_st)~class,data=class_surv)$chisq,i_clust-1)
        class_surv_table = table(class_surv[,"class"],class_surv[,"surv_st"])
        risk_label = paste0("Cl",1:i_clust," (n",apply(class_surv_table,1,sum),"/d",class_surv_table[,2],"/a",class_surv_table[,1],")")
        plot(fit,col=risk_color[1:i_clust],xlab="survival time [days from Dx]",ylab="survival probability",conf.int=F,lty=1,lwd=2,mark.time=T)
        title = paste0("mad",mad_T_thr[i_thres]," (n",nrow(class_surv),"/d",table(class_surv[,"surv_st"])[2],"/a",table(class_surv[,"surv_st"])[1],
        ") ; logrank p-val=",signif(logrank_pval,digits=2))
        title(title,cex.main=1,line=1)
        legend(35,0.5,risk_label,col=risk_color,lty=1,lwd=2,cex=0.8)
    }
    dev.off()
}
