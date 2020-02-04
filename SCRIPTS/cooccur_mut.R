library("gplots")

rm(list=ls())

PROJECT_DIR = "~/ACTIVE/LCS/Mongolia_HCC_Pipeline" # replace this line with your local path

infile = file.path(PROJECT_DIR,"DATA","PROCESSED","maf_only_genes_in_oncoplot.maf")
maf_data <- read.table(infile,sep="\t",quote="",stringsAsFactors=F,header=T) # 146 alterations
mutated_samples = unique(maf_data[,"Tumor_Sample_Barcode"]) # 57/71 samples with mutations

infile = file.path(PROJECT_DIR,"DATA","PROCESSED","patient_sample_metadata_w_clustering_risk.txt")
metadata = as.matrix(read.table(infile,header=T,check.names=F,stringsAsFactors=F,sep="\t"))
metadata = metadata[!is.na(metadata[,"WES_T"]),]
n_sample = nrow(metadata) # 71 samples

# we add class-specific binary variables
class = levels(factor(metadata[,"class"]))
n_class = length(class)
for (i_class in 1:n_class) {
    class_bin = rep("Undef",n_sample)
    class_bin[which(metadata[,"class"]==class[i_class])] = 1
    class_bin[which(metadata[,"class"]%in%class[-i_class])] = 0
    metadata = cbind(metadata,rep(NA,nrow(metadata)))
    colnames(metadata)[ncol(metadata)] = paste0("class_",class[i_class])
    metadata[,paste0("class_",class[i_class])] = class_bin
}

class_bin = rep("Undef",n_sample)
class_bin[c(which(metadata[,"class_MO1"]=="1"),which(metadata[,"class_MO2"]=="1"))] = 1
class_bin[c(which(metadata[,"class_MO3"]=="1"),which(metadata[,"class_MO4"]=="1"))] = 0
metadata = cbind(metadata,rep(NA,nrow(metadata)))
colnames(metadata)[ncol(metadata)] = "class_MO1_MO2"
metadata[,"class_MO1_MO2"] = class_bin

class_bin = rep("Undef",n_sample)
class_bin[c(which(metadata[,"class_MO1"]=="1"),which(metadata[,"class_MO2"]=="1"))] = 0
class_bin[c(which(metadata[,"class_MO3"]=="1"),which(metadata[,"class_MO4"]=="1"))] = 1
metadata = cbind(metadata,rep(NA,nrow(metadata)))
colnames(metadata)[ncol(metadata)] = "class_MO3_MO4"
metadata[,"class_MO3_MO4"] = class_bin

metadata[is.na(metadata[,"risk_bin"]),"risk_bin"] = "Undef"
metadata[,"risk_bin"] = trimws(metadata[,"risk_bin"])

group = c("all","low_risk","high_risk","MO1","MO2","MO3","MO4","MO1_MO2","MO3_MO4")
for (i_group in 1:length(group)) {
    if (group[i_group]=="all") {
        select = rep(T,nrow(metadata))
    } else if (group[i_group]=="low_risk") {
        select = metadata[,"risk_bin"]=="0"
    } else if (group[i_group]=="high_risk") {
        select = metadata[,"risk_bin"]=="1"
    } else if (group[i_group]=="MO1") {
        select = metadata[,"class_MO1"]=="1"
    } else if (group[i_group]=="MO2") {
        select = metadata[,"class_MO2"]=="1"
    } else if (group[i_group]=="MO3") {
        select = metadata[,"class_MO3"]=="1"
    } else if (group[i_group]=="MO4") {
        select = metadata[,"class_MO4"]=="1"
    } else if (group[i_group]=="MO1_MO2") {
        select = metadata[,"class_MO1_MO2"]=="1"
    } else if (group[i_group]=="MO3_MO4") {
        select = metadata[,"class_MO3_MO4"]=="1"
    }
    title = c("gene","gene","p-value","odds")
    gene = unique(maf_data[maf_data[,"Tumor_Sample_Barcode"]%in%mutated_samples[mutated_samples%in%metadata[select,"WES_T"]],"Hugo_Symbol"])
    n_gene = length(gene)
    res = matrix(rep(NA,n_gene*(n_gene-1)/2*length(title)),ncol=length(title))
    index = 0
    for (i_gene in 1:(n_gene-1)) {
        select1 = metadata[select,"WES_T"]%in% unique(as.character(maf_data[maf_data[,"Hugo_Symbol"]==gene[i_gene],"Tumor_Sample_Barcode"]))
        for (j_gene in (i_gene+1):n_gene) {
            select2 = metadata[select,"WES_T"]%in% unique(as.character(maf_data[maf_data[,"Hugo_Symbol"]==gene[j_gene],"Tumor_Sample_Barcode"]))
            dt = table(select1,select2)
            fisher.pval = fisher.test(dt)$p.value
            odds = fisher.test(dt)$estimate
            index = index+1
            res[index,] = c(gene[i_gene],gene[j_gene],fisher.pval,odds)
        }
    }
    outfile = file.path(PROJECT_DIR,"RESULTS","mut_driver",paste0("cooccur_mut_",group[i_group],".txt"))
    output = rbind(title,res)
    write(t(output),ncol=ncol(output),file=outfile,sep="\t")
}
