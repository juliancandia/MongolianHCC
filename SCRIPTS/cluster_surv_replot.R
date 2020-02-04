library(survival)

rm(list=ls())

PROJECT_DIR = "~/ACTIVE/LCS/Mongolia_HCC_Pipeline" # replace this line with your local path

setwd(file.path(PROJECT_DIR,"RESULTS","cluster","FINAL"))

class_color = c("#4DAF4A","#377EB8","orange","#E41A1C")
best_K = 4 # number of classes

infile = "patient_sample_metadata_w_clustering.txt"
data = as.matrix(read.table(infile,header=T,check.names=F,stringsAsFactors=F,sep="\t"))
select = (!is.na(data[,"RNASeq_ID"]))&(!is.na(data[,"survival.status"]))&(!is.na(data[,"survival.time"]))
data = data[select,]

# Kaplan-Meier plots and log-rank test
outfile = "KaplanMeier.pdf"
pdf(outfile,width=7,height=5)
class_surv = data.frame(as.numeric(data[,"survival.time"]),
as.numeric(data[,"survival.status"]),data[,"class"])
colnames(class_surv) = c("surv_time","surv_st","class")
class_surv[,"surv_time"] = class_surv[,"surv_time"]/(365.25/12) # convert time to months
fit <- survfit(Surv(surv_time,surv_st)~class, data=class_surv)
logrank_pval = 1-pchisq(survdiff(Surv(surv_time,surv_st)~class,data=class_surv)$chisq,best_K-1)
class_surv_table = table(class_surv[,"class"],class_surv[,"surv_st"])
#risk_label = paste0(levels(class_surv[,"class"]),": n=",apply(class_surv_table,1,sum)," (dead=",class_surv_table[,2],", alive=",class_surv_table[,1],")")
risk_label = paste0(levels(class_surv[,"class"])," (n=",apply(class_surv_table,1,sum),")")

plot(fit,xlim=range(0,50),col=class_color[1:best_K],xlab="Months from Diagnosis",ylab="Overall Survival (%)",conf.int=F,lty=1,lwd=3,mark.time=T,xaxt="n",yaxt="n",cex.lab=1.25)

#title = paste0("(n",nrow(class_surv),"/d",table(class_surv[,"surv_st"])[2],"/a",table(class_surv[,"surv_st"])[1],
#") ; logrank p-val=",signif(logrank_pval,digits=2))
axis(side=1,at=c(0,12,24,36,48),labels=c(0,12,24,36,48),las=1,cex.axis=1)
axis(side=1,at=c(6,18,30,42),labels=rep("",4),las=1,cex.axis=1,tck=-0.01)
axis(side=2,at=seq(0,1,by=0.2),labels=seq(0,100,by=20),las=2,cex.axis=1)
legend(3,0.38,risk_label,col=class_color,lty=1,lwd=3,cex=1)
dev.off()

class = data[,"class"]
class_label = names(table(class))
n_class = length(class_label)
logrank_pval_pair = matrix(rep(NA,n_class**2),ncol=n_class)
rownames(logrank_pval_pair) = class_label
colnames(logrank_pval_pair) = class_label
for (i_class in 1:(n_class-1)) {
    for (j_class in (i_class+1):n_class) {
        class_surv_pair = class_surv[class_surv[,"class"]%in%class_label[c(i_class,j_class)],]
        logrank_pval_pair[i_class,j_class] = 1-pchisq(survdiff(Surv(surv_time,surv_st)~class,data=class_surv_pair)$chisq,1)
        logrank_pval_pair[j_class,i_class] = logrank_pval_pair[i_class,j_class]
    }
}
outfile = "KaplanMeier_logrank_pairwise.txt"
output = rbind(class_label,logrank_pval_pair)
write(t(output),ncol=ncol(output),file=outfile,sep="\t")

# high vs low survival
class_H_L = rep("H",length(class))
class_H_L[class%in%c("MO3","MO4")] = "L"
class_surv = data.frame(class_surv,class_H_L)
pval = 1-pchisq(survdiff(Surv(surv_time,surv_st)~class_H_L,data=class_surv)$chisq,1)
outfile = "KaplanMeier_MO1MO2_vs_MO3MO4.txt"
write(pval,ncol=1,file=outfile,sep="\t")
