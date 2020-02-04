library("gplots")

rm(list=ls())

PROJECT_DIR = "~/ACTIVE/LCS/Mongolia_HCC_Pipeline" # replace this line with your local path

infile = file.path(PROJECT_DIR,"DATA","PROCESSED","maf_only_genes_in_oncoplot.maf")
maf_data <- read.table(infile, sep="\t", quote="", header=T) # 146 alterations
#n_sample = length(unique(maf_data[,"Tumor_Sample_Barcode"])) # 57/71 samples with mutations
gene = as.character(unique(maf_data[,"Hugo_Symbol"]))
n_gene = length(gene) # 19 genes

infile = file.path(PROJECT_DIR,"DATA","PROCESSED","patient_sample_metadata_w_clustering_risk.txt")
metadata = as.matrix(read.table(infile,header=T,check.names=F,stringsAsFactors=F,sep="\t"))
metadata = metadata[!is.na(metadata[,"WES_T"]),]
n_sample = nrow(metadata) # 71 samples

# we add class-specific binary variables
class = levels(factor(metadata[,"class"]))
n_class = length(class)
for (i_class in 1:n_class) {
    class_bin = rep(NA,n_sample)
    class_bin[which(metadata[,"class"]==class[i_class])] = 1
    class_bin[which(metadata[,"class"]%in%class[-i_class])] = 0
    metadata = cbind(metadata,rep(NA,nrow(metadata)))
    colnames(metadata)[ncol(metadata)] = paste0("class_",class[i_class])
    metadata[,paste0("class_",class[i_class])] = class_bin
}

var = c("age_bin","sex","hcv","hbv","hdv","stage","tumor_size","cirrhosis","obesity","alcohol","smoker","FHx_LC","alb","bil","alt","afp","multinodular","risk_bin",paste0("class_",class))
var_title = c("Age","Sex","HCV","HBV","HDV","Stage","Tumor Size","Cirrhosis","Obesity","Alcohol","Smoker","Fam Hx Liv C","Albumin","Bilirubin","ALT","AFP","Multinodular","Risk Group",paste("Class",class))
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
var_label[[18]] = c("low risk","high risk")
var_label[[19]] = paste0(c("not ",""),class[1])
var_label[[20]] = paste0(c("not ",""),class[2])
var_label[[21]] = paste0(c("not ",""),class[3])
var_label[[22]] = paste0(c("not ",""),class[4])

# Fisher's exact tests
outfile = file.path(PROJECT_DIR,"RESULTS","mut_driver","assoc_mut_annot.pdf")
pdf(outfile,width=7,height=5)
title = c("gene","annotation","p-value","odds")
res = matrix(rep(NA,n_gene*n_var*length(title)),ncol=length(title))
for (i_gene in 1:n_gene) {
    sel_gene = metadata[,"WES_T"]%in%unique(as.character(maf_data[maf_data[,"Hugo_Symbol"]==gene[i_gene],"Tumor_Sample_Barcode"]))
    gene_mut = rep(0,n_sample)
    gene_mut[sel_gene] = 1
    for (i_var in 1:n_var) {
        select = !is.na(metadata[,var[i_var]])
        if ((sum(select&sel_gene)==0)|(sum(select&(!sel_gene))==0)) {
            res[(i_gene-1)*n_var+i_var,] = c(gene[i_gene],var_title[i_var],NA,NA)
        } else {
            dt = table(gene_mut[select],metadata[select,var[i_var]])
            rownames(dt) = paste0(c("not ",""),"mutated")
            colnames(dt) = var_label[[i_var]]
            fisher.pval = fisher.test(dt)$p.value
            odds = fisher.test(dt)$estimate
            res[(i_gene-1)*n_var+i_var,] = c(gene[i_gene],var_title[i_var],fisher.pval,odds)
            balloonplot(t(dt),main="",xlab="",ylab="",label=T,show.margins=F)
            title(paste(gene[i_gene],"vs",var_title[i_var]),line=2)
            par(cex=0.85)
            title(paste0("Fisher's exact test p-value = ",signif(fisher.pval,digits=2),
            " ; odds ratio = ",signif(odds,digits=2)),line=0.5,cex.main=0.95)
        }
    }
}
dev.off()
outfile = file.path(PROJECT_DIR,"RESULTS","mut_driver","assoc_mut_annot.txt")
#output = rbind(title,res)
output = rbind(c("plot",title),cbind(1:nrow(res),res)) # we add plot (page) numbers
write(t(output),ncol=ncol(output),file=outfile,sep="\t")
