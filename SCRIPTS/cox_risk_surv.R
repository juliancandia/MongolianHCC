library(survival)

rm(list=ls())

PROJECT_DIR = "~/ACTIVE/LCS/Mongolia_HCC_Pipeline" # replace this line with your local path

setwd(file.path(PROJECT_DIR,"RESULTS","cox"))

infile = file.path(PROJECT_DIR,"RESULTS","cluster","FINAL","patient_sample_metadata_w_clustering_risk.txt")
data = as.matrix(read.table(infile,header=T,check.names=F,stringsAsFactors=F,sep="\t"))
select = !is.na(data[,"risk"])
data = data[select,]
class_color = c("#4DAF4A","#E41A1C")

for (method in c("risk_bin","risk")) {
    # Kaplan-Meier plots and log-rank test
    outfile = paste0("KaplanMeier_",method,".pdf")
    pdf(outfile,width=7,height=5)
    class = trimws(data[,method])
    if (method=="risk") {
        tmp = as.numeric(class)
        class = rep("1",length(class))
        class[tmp<=median(tmp)] = "0"
    }
    class_surv = data.frame(as.numeric(data[,"survival.time"]),
    as.numeric(data[,"survival.status"]),class)
    colnames(class_surv) = c("surv_time","surv_st","class")
    class_surv[,"surv_time"] = class_surv[,"surv_time"]/(365.25/12) # convert time to months
    fit <- survfit(Surv(surv_time,surv_st)~class, data=class_surv)
    logrank_pval = 1-pchisq(survdiff(Surv(surv_time,surv_st)~class,data=class_surv)$chisq,1)
    class_surv_table = table(class_surv[,"class"],class_surv[,"surv_st"])
    #risk_label = paste0(levels(class_surv[,"class"]),": n=",apply(class_surv_table,1,sum)," (dead=",class_surv_table[,2],", alive=",class_surv_table[,1],")")
    #risk_label = paste0(levels(class_surv[,"class"])," (n=",apply(class_surv_table,1,sum),")")
    risk_label = paste0(paste(c("Low","High"),"Risk")," (n=",apply(class_surv_table,1,sum),")")

    plot(fit,xlim=range(0,56),col=class_color,xlab="Months from Diagnosis",ylab="Overall Survival (%)",conf.int=F,lty=1,lwd=3,mark.time=T,xaxt="n",yaxt="n",cex.lab=1.25)

    title = paste0("(n",nrow(class_surv),"/d",table(class_surv[,"surv_st"])[2],"/a",table(class_surv[,"surv_st"])[1],
    ") ; logrank p-val=",signif(logrank_pval,digits=2))
    title(title,cex.main=1,line=1)
    axis(side=1,at=c(0,12,24,36,48),labels=c(0,12,24,36,48),las=1,cex.axis=1)
    axis(side=1,at=c(6,18,30,42),labels=rep("",4),las=1,cex.axis=1,tck=-0.01)
    axis(side=2,at=seq(0,1,by=0.2),labels=seq(0,100,by=20),las=2,cex.axis=1)
    legend(3,0.28,risk_label,col=class_color,lty=1,lwd=3,cex=1)

    # we plot again without the title
    plot(fit,xlim=range(0,56),col=class_color,xlab="Months from Diagnosis",ylab="Overall Survival (%)",conf.int=F,lty=1,lwd=3,mark.time=T,xaxt="n",yaxt="n",cex.lab=1.25)
    axis(side=1,at=c(0,12,24,36,48),labels=c(0,12,24,36,48),las=1,cex.axis=1)
    axis(side=1,at=c(6,18,30,42),labels=rep("",4),las=1,cex.axis=1,tck=-0.01)
    axis(side=2,at=seq(0,1,by=0.2),labels=seq(0,100,by=20),las=2,cex.axis=1)
    legend(3,0.28,risk_label,col=class_color,lty=1,lwd=3,cex=1)
    dev.off()
}
