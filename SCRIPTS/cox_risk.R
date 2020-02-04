library(survival)

rm(list=ls())

PROJECT_DIR = "~/ACTIVE/LCS/Mongolia_HCC_Pipeline" # replace this line with your local path

setwd(file.path(PROJECT_DIR,"RESULTS","cox"))

load("eNet.Robj") # loads eNetXplorer object
fit = eNet
predictor = as.matrix(fit$predictor)
rownames(predictor) = fit$instance
response = fit$response
feat_pval = as.matrix(fit$feature_freq_model_vs_null_pval)
pval_thres = 0.05
feat_index = order(feat_pval[,"a1"])[1:sum(feat_pval[,"a1"]<pval_thres)] # we choose the top features (based on p-value for the lasso solution).

x = predictor[,feat_index]
myData = data.frame(response,x)
myFormula = as.formula(paste0("Surv(time, status) ~ ",paste(colnames(myData)[-c(1:2)],collapse="+")))
cox.fit = coxph(myFormula,data=myData)
coeff = cox.fit$coefficients
#score = predict(cox.fit,type="lp")
score = predict(cox.fit,type="risk") # risk = exp(lp)

outfile = "risk_scores_distrib.pdf"
pdf(outfile,width=7,height=5)
plot(density(score),xlab="Risk Score",ylab="Probability Density",main="",cex.lab=1.25)
score_thres = 2.77 # based on density plot
abline(v=score_thres,lty=3)
dev.off()

infile = file.path(PROJECT_DIR,"RESULTS","cluster","FINAL","patient_sample_metadata_w_clustering.txt")
metadata = as.matrix(read.table(infile,header=T,check.names=F,stringsAsFactors=F,sep="\t"))
metadata2 = metadata[!is.na(metadata[,"RNASeq_ID"]),] # only patients with RNASeq data
# we check that subject IDs match
check_passed = (nrow(x)==nrow(metadata2))&
(sum(trimws(metadata2[,"Patient"])==rownames(x))==nrow(x))
cat("Consistency check passed: ",check_passed,"\n")
if (!check_passed) {stop("Rows don't match")}

risk = rep(NA,nrow(metadata))
risk[!is.na(metadata[,"RNASeq_ID"])] = score

score_bin = rep(NA,length(score))
score_bin[score<score_thres] = 0
score_bin[score>=score_thres] = 1
risk_bin = rep(NA,nrow(metadata))
risk_bin[!is.na(metadata[,"RNASeq_ID"])] = score_bin

output = rbind(c(colnames(metadata),"risk","risk_bin"),cbind(metadata,risk,risk_bin))
outfile = file.path(PROJECT_DIR,"RESULTS","cluster","FINAL","patient_sample_metadata_w_clustering_risk.txt")
write(t(output),ncol=ncol(output),file=outfile,sep="\t")

# we also place a copy in the DATA/PROCESSED folder
outfile = file.path(PROJECT_DIR,"DATA","PROCESSED","patient_sample_metadata_w_clustering_risk.txt")
write(t(output),ncol=ncol(output),file=outfile,sep="\t")
