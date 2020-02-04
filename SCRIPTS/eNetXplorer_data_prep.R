library(missForest)

rm(list=ls())

PROJECT_DIR = "~/ACTIVE/LCS/Mongolia_HCC_Pipeline" # replace this line with your local path

infile = file.path(PROJECT_DIR,"DATA","PROCESSED","gene_expr_voom_T.txt")
expr = as.matrix(read.table(infile,header=T,check.names=F,stringsAsFactors=F,sep="\t"))
gene = expr[,1]
patient = colnames(expr[,-c(1,2)])
n_pat = length(patient)
expr = t(matrix(as.numeric(expr[,-c(1,2)]),ncol=ncol(expr)-2))

infile = file.path(PROJECT_DIR,"DATA","ORIGINAL","patient_sample_metadata.txt")
metadata = read.table(infile,header=T,check.names=F,stringsAsFactors=F,sep="\t")
index = rep(NA,n_pat)
for (i_pat in 1:n_pat) {
    index[i_pat] = which(metadata[,"Patient"]==patient[i_pat])
}
metadata = metadata[index,c("age","sex","hcv","hbv","hdv","survival.status","survival.time")]
# NOTE: other clinical variables (such as  "stage", "tumor_size" and "cirrhosis") have too many missing values to be useful for this analysis.

select = (!is.na(metadata[,"survival.status"]))&(!is.na(metadata[,"survival.time"]))
metadata = metadata[select,]
expr = expr[select,]
patient = patient[select]

# we read pathway info
infile = file.path(PROJECT_DIR,"DATA","Reactome","ReactomePathways.gmt")
con  = file(infile,open="r")
pathway_name = NULL
pathway_id = NULL
pathway_genes = list()
dataList = list()
while (length(oneLine <- readLines(con, n=1, warn=FALSE)) > 0) {
    myVector <- strsplit(oneLine, "\t")
    pathway_name = c(pathway_name,myVector[[1]][1])
    pathway_id = c(pathway_id,myVector[[1]][2])
    pathway_genes <- c(pathway_genes,list(myVector[[1]][-c(1,2)]))
}
close(con)
n_pathway = length(pathway_id)

infile = file.path(PROJECT_DIR,"DATA","Reactome","reactome_metadata.txt")
pathway_metadata = read.table(infile,header=T,sep="\t")
pathway_class = rep("NA",n_pathway)
for (i_pathway in 1:n_pathway) {
    index = which(as.character(pathway_metadata[,"reactomeID"])==pathway_id[i_pathway])
    if (length(index)==1) {
        pathway_class[i_pathway] = as.character(pathway_metadata[index,"categoryNameAbbr"])
    }
}
pathway_label = paste0(substring(pathway_id,7)," [",pathway_class,"]")

# PC1 expression on pathways
pathway_expr = matrix(rep(NA,nrow(expr)*n_pathway),ncol=n_pathway)
pathway_stats_header = c("pathway_size","genes_found_n","genes_found_frac","genes_found","genes_not_found")
pathway_stats = matrix(rep(NA,n_pathway*length(pathway_stats_header)),ncol=length(pathway_stats_header))
colnames(pathway_stats) = pathway_stats_header
for (i_pathway in 1:n_pathway) {
    index = which(gene%in%pathway_genes[[i_pathway]])
    pathway_stats[i_pathway,"pathway_size"] = length(pathway_genes[[i_pathway]])
    pathway_stats[i_pathway,"genes_found_n"] = length(unique(gene[index]))
    pathway_stats[i_pathway,"genes_found_frac"] = as.numeric(pathway_stats[i_pathway,"genes_found_n"])/as.numeric(pathway_stats[i_pathway,"pathway_size"])
    pathway_stats[i_pathway,"genes_found"] =
    paste(pathway_genes[[i_pathway]][pathway_genes[[i_pathway]]%in%unique(gene[index])],collapse="|")
    pathway_stats[i_pathway,"genes_not_found"] =
    paste(pathway_genes[[i_pathway]][!pathway_genes[[i_pathway]]%in%unique(gene[index])],collapse="|")
    if (length(index)>1) {
        PC1 = (prcomp(expr[,index],scale=F))$x[,1]
        correl = NULL
        for (i in 1:length(index)) {
            correl = c(correl,cor(expr[,index[i]],PC1))
        }
        pathway_expr[,i_pathway] = sign(median(correl))*PC1
    } else if (length(index)==1) {
        pathway_expr[,i_pathway] = expr[,index]
    }
}
pathway_select = pathway_stats[,"genes_found_n"]>0
pathway_expr = pathway_expr[,pathway_select]
pathway_stats = pathway_stats[pathway_select,]
pathway_label = pathway_label[pathway_select]
pathway_id = pathway_id[pathway_select]
pathway_name = pathway_name[pathway_select]
pathway_class = pathway_class[pathway_select]

# we prep the input data to feed into eNetXplorer
x = cbind(metadata[,c("age","sex","hcv","hbv","hdv")],pathway_expr)
set.seed(123) # imputation algorithm is stochastic
x = missForest(x)$ximp # to impute missing values in clinical/demographics data
y = metadata[,c("survival.time","survival.status")]

destination = file.path(PROJECT_DIR,"RESULTS","cox")
dir.create(destination)
setwd(destination)

write(t(x),ncol=ncol(x),file="data.txt",sep="\t")
write(patient,ncol=1,file="instances.txt",sep="\t")
predictors = c("age","sex","hcv","hbv","hdv",pathway_label)
write(predictors,ncol=1,file="predictors.txt",sep="\t")
output = rbind(c("label","id","name","class",pathway_stats_header),cbind(pathway_label,pathway_id,pathway_name,pathway_class,pathway_stats))
write(t(output),ncol=ncol(output),file="pathway_metadata.txt",sep="\t")
output = rbind(c("time","status"),y)
write(t(output),ncol=ncol(output),file="response.txt",sep="\t")
