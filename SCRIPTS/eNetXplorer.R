library(eNetXplorer)

rm(list=ls())

PROJECT_DIR = "~/ACTIVE/LCS/Mongolia_HCC_Pipeline" # replace this line with your local path

setwd(file.path(PROJECT_DIR,"RESULTS","cox"))

data = as.matrix(read.table("data.txt",header=F,check.names=F,stringsAsFactors=F,sep="\t"))
instances = as.matrix(read.table("instances.txt",header=F,check.names=F,stringsAsFactors=F,sep="\t"))[,1]
predictors = as.matrix(read.table("predictors.txt",header=F,check.names=F,stringsAsFactors=F,sep="\t"))[,1]
response = as.matrix(read.table("response.txt",header=T,check.names=F,stringsAsFactors=F,sep="\t"))

x = data
rownames(x) = instances
colnames(x) = predictors
y = response
fit = eNetXplorer(x=x,y=y,family="cox",alpha=seq(0,1,by=0.1),n_run=1000,n_perm_null=250,seed=123,save_obj=T,dest_obj="eNet.Robj")
summaryPDF(fit, dest_file="eNet.pdf")
export(fit, dest_dir="eNet")
