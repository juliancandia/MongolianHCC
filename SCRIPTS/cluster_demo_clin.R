library("gplots")

rm(list=ls())

PROJECT_DIR = "~/ACTIVE/LCS/Mongolia_HCC_Pipeline" # replace this line with your local path

setwd(file.path(PROJECT_DIR,"RESULTS","cluster","FINAL"))

infile = "patient_sample_metadata_w_clustering.txt"
data = as.matrix(read.table(infile,header=T,check.names=F,stringsAsFactors=F,sep="\t"))
data = data[!is.na(data[,"RNASeq_ID"]),]
n_subj = nrow(data)
class = data[,"class"]
class_label = names(table(class))
n_class = length(class_label)

var = c("age_bin","sex","hcv","hbv","hdv","stage","tumor_size","cirrhosis","obesity","alcohol","smoker","FHx_LC","alb","bil","alt","afp","multinodular")
var_title = c("Age","Sex","HCV","HBV","HDV","Stage","Tumor Size","Cirrhosis","Obesity","Alcohol","Smoker","Fam Hx Liv C","Albumin","Bilirubin","ALT","AFP","Multinodular")
n_var = length(var)
var_label = vector("list",n_var)
var_label[[1]] = c("below median","above median")
var_label[[2]] = c("female","male")
var_label[[3]] = c("no","yes")
var_label[[4]] = c("no","yes")
var_label[[5]] = c("no","yes")
var_label[[6]] = c("early","advanced")
var_label[[7]] = c("small","large")
var_label[[8]] = c("no","yes")
var_label[[9]] = c("no","yes")
var_label[[10]] = c("no","yes")
var_label[[11]] = c("no","yes")
var_label[[12]] = c("no","yes")
var_label[[13]] = c("normal","abnormal")
var_label[[14]] = c("normal","abnormal")
var_label[[15]] = c("normal","abnormal")
var_label[[16]] = c("normal","abnormal")
var_label[[17]] = c("no","yes")

# Fisher's exact tests (All Classes)
pdf("DemoClin_All_Classes.pdf",width=5,height=5)
fisher.pval = NULL
for (i_var in 1:n_var) {
    var_subj = as.numeric(data[,var[i_var]])
    var_subj[var_subj==0] = var_label[[i_var]][1]
    var_subj[var_subj==1] = var_label[[i_var]][2]
    select = !is.na(var_subj)
    dt = table(class[select],var_subj[select])
    pval = fisher.test(dt)$p.value
    fisher.pval = c(fisher.pval,pval)
    balloonplot(t(dt),main="",xlab="",ylab="",label=T,show.margins=F)
    title(var_title[i_var],line=2)
    title(paste0("Fisher's exact test p-value = ",signif(pval,digits=2)),line=0.5,cex.main=0.95)
}
dev.off()

o = order(fisher.pval)
res = cbind(var_title,fisher.pval)[o,]
output = rbind(c("Demographic/Clinical Variable","Fisher pval"),res)
write(t(output),ncol=ncol(output),file="DemoClin_All_Classes_pval.txt",sep="\t")

##############
# Fisher's exact tests
pdf("DemoClin_One_Class.pdf",width=7,height=5)
fisher.pval = matrix(rep(NA,n_var*n_class),ncol=n_class)
fisher.odds = matrix(rep(NA,n_var*n_class),ncol=n_class)
for (i_var in 1:n_var) {
    if (var[i_var]=="age") {
        var_subj = as.numeric(data[,var[i_var]])
        thres = median(var_subj,na.rm=T)
        var_subj[which(var_subj<thres)] = 0
        var_subj[which(var_subj>=thres)] = 1
    } else {
        var_subj = as.numeric(data[,var[i_var]])
    }
    select = !is.na(var_subj)
    for (i_class in 1:n_class) {
        var_class = rep(0,n_subj)
        var_class[class==class_label[i_class]] = 1
        dt = table(var_class[select],var_subj[select])
        rownames(dt) = paste0(c("not ",""),class_label[i_class])
        colnames(dt) = var_label[[i_var]]
        fisher.pval[i_var,i_class] = fisher.test(dt)$p.value
        fisher.odds[i_var,i_class] = fisher.test(dt)$estimate
        balloonplot(t(dt),main="",xlab="",ylab="",label=T,show.margins=F)
        title(paste0(var_title[i_var]," / ",class_label[i_class]),line=2)
        title(paste0("Fisher's exact test p-value = ",signif(fisher.pval[i_var,i_class],digits=2),
        " ; odds ratio = ",signif(fisher.odds[i_var,i_class],digits=2)),line=0.5,cex.main=0.95)
    }
}
dev.off()

res = cbind(var_title,fisher.pval)[o,]
output = rbind(c("Demographic/Clinical Variable",class_label),res)
write(t(output),ncol=ncol(output),file="DemoClin_One_Class_pval.txt",sep="\t")

res = cbind(var_title,fisher.odds)[o,]
output = rbind(c("Demographic/Clinical Variable",class_label),res)
write(t(output),ncol=ncol(output),file="DemoClin_One_Class_odds.txt",sep="\t")

# binarized classes
class_bin = rep(NA,n_subj)
class_bin[(class=="MO1")|(class=="MO2")] = "MO1_MO2"
class_bin[(class=="MO3")|(class=="MO4")] = "MO3_MO4"
class_bin_label = names(table(class_bin))

# Fisher's exact tests
pdf("DemoClin_MO1MO2_MO3MO4.pdf",width=7,height=5)
fisher_res = matrix(rep(NA,n_var*2),ncol=2)
for (i_var in 1:n_var) {
    if (var[i_var]=="age") {
        var_subj = as.numeric(data[,var[i_var]])
        thres = median(var_subj,na.rm=T)
        var_subj[which(var_subj<thres)] = 0
        var_subj[which(var_subj>=thres)] = 1
    } else {
        var_subj = as.numeric(data[,var[i_var]])
    }
    select = !is.na(var_subj)
    dt = table(class_bin[select],var_subj[select])
    colnames(dt) = var_label[[i_var]]
    fisher_res[i_var,1] = fisher.test(dt)$p.value
    fisher_res[i_var,2] = fisher.test(dt)$estimate
    balloonplot(t(dt),main="",xlab="",ylab="",label=T,show.margins=F)
    title(var_title[i_var],line=2)
    title(paste0("Fisher's exact test p-value = ",signif(fisher_res[i_var,1],digits=2),
    " ; odds ratio = ",signif(fisher_res[i_var,2],digits=2)),line=0.5,cex.main=0.95)
}
dev.off()

o = order(fisher_res[,1])
res = cbind(var_title,fisher_res)[o,]
output = rbind(c("Demographic/Clinical Variable","p-value","odds"),res)
write(t(output),ncol=ncol(output),file="DemoClin_MO1MO2_MO3MO4.txt",sep="\t")
