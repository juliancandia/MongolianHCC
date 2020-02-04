library(calibrate)
library(survival)

rm(list=ls())

PROJECT_DIR = "~/ACTIVE/LCS/Mongolia_HCC_Pipeline" # replace this line with your local path

setwd(file.path(PROJECT_DIR,"RESULTS","cox"))

class_color = c("#4DAF4A","#377EB8","orange","#E41A1C")

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

infile = file.path(PROJECT_DIR,"RESULTS","cluster","FINAL","patient_sample_metadata_w_clustering.txt")
metadata = as.matrix(read.table(infile,header=T,check.names=F,stringsAsFactors=F,sep="\t"))
metadata = metadata[!is.na(metadata[,"RNASeq_ID"]),] # only patients with RNASeq data
# we check that subject IDs match
check_passed = (nrow(x)==nrow(metadata))&
(sum(trimws(metadata[,"Patient"])==rownames(x))==nrow(x))
cat("Consistency check passed: ",check_passed,"\n")
if (!check_passed) {stop("Rows don't match")}

# boxplots
outfile = "risk_scores_boxplot.pdf"
pdf(outfile,width=5,height=5)
class = metadata[,"class"]
class_label = names(table(class))
n_class = length(class_label)
class_size = table(class)
BoxplotData = data.frame(class,score)
sample_label = rownames(x)
#boxplot(as.formula("score ~ class"),data=BoxplotData,main="",xlab="class", ylab="risk scores",outline=F,boxwex=0.8,range=0,las=2)
boxplot(as.formula("score ~ class"),data=BoxplotData,main="",xlab="", ylab="Risk Scores (Cox Regression)",outline=F,boxwex=0.8,range=0,las=2,xaxt="n",log="y")
#mtext(side=2, text="Risk Scores (Cox Regression)", line=2.3)
axis(side=3,at=1:n_class,labels=rep("",n_class),las=2,cex.axis=1) # to add tick labels
mtext(side=3, text=class_label , at=1:n_class, las=1, line=1, col=class_color, cex=1) # line=distance between axis and tick labels

#library(plotteR)
# Example:
#ddf = data.frame(NUMS = rnorm(500), GRP = sample(LETTERS[1:5],500,replace=T))
#boxplot(NUMS ~ GRP, data = ddf, lwd = 2, ylab = 'NUMS')
#spreadPointsMultiple(data=ddf, responseColumn="NUMS", categoriesColumn="GRP",
#col="blue", plotOutliers=TRUE)
#spreadPointsMultiple(data=BoxplotData, responseColumn="score", categoriesColumn="class",
#col=class_color, plotOutliers=TRUE)
transparency = 80
scatter_span = 0.35
wilcox_pval = matrix(rep(NA,n_class**2),ncol=n_class)
ttest_pval = matrix(rep(NA,n_class**2),ncol=n_class)
for (i_class in 1:n_class) {
    scatter = runif(class_size[i_class],-scatter_span,scatter_span)
    x = i_class+scatter-mean(scatter)
    y = BoxplotData[BoxplotData$class==class_label[i_class],"score"]
    col_par = as.numeric(col2rgb(class_color[i_class]))
    #points(x,y,cex=1.5,col=rgb(col_par[1],col_par[2],col_par[3],transparency,maxColorValue=255),pch=16)
    points(x,y,col=class_color[i_class],cex=1.5,pch=16)
    #textxy(x,y,sample_label[BoxplotData$class==class_label[i_class]],cex=0.45,offset=0.6)
    for (j_class in 1:n_class) {
        y2 = BoxplotData[BoxplotData$class==class_label[j_class],"score"]
        wilcox_pval[i_class,j_class] = wilcox.test(y,y2)$p.value
        ttest_pval[i_class,j_class] = t.test(y,y2)$p.value
    }
}
dev.off()

output = rbind(class_label,wilcox_pval)
outfile = "risk_scores_wilcox.txt"
write(t(output),ncol=ncol(output),file=outfile,sep="\t")

output = rbind(class_label,ttest_pval)
outfile = "risk_scores_ttest.txt"
write(t(output),ncol=ncol(output),file=outfile,sep="\t")

output = rbind(c("feature","coeff"),cbind(names(coeff),coeff))
outfile = "risk_coeff.txt"
write(t(output),ncol=ncol(output),file=outfile,sep="\t")
