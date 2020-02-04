rm(list=ls())

PROJECT_DIR = "~/ACTIVE/LCS/Mongolia_HCC_Pipeline" # replace this line with your local path

infile = file.path(PROJECT_DIR,"DATA","PROCESSED","maf_only_genes_in_oncoplot.maf")
mo <- read.table(infile, sep="\t", quote="", header=T)

n_MO = 71 # total number of samples analyzed

# analyze only candidate "Mongolia genes"
target = c("TP53","CTNNB1","ALB","GTF2IRD2B","RB1","COL11A1","PNRC2","AK2","VPS13A","SPTA1",            "PCLO","CSMD2","APOB","SMC6","CSMD3","DYNC2H1","FKBP9","LRP1B","PCDH7") # same order as in the oncoplot.
n_target = length(target)

MO_freq = rep(NA,n_target)
for (i_target in 1:n_target) {
    
    MO_freq[i_target] = length(unique(mo[mo[,"Hugo_Symbol"]==target[i_target],"Tumor_Sample_Barcode"]))/n_MO
}

infile =  file.path(PROJECT_DIR,"RESULTS","mut_driver","mut_freq_TCGA_byGene.txt")
data2 = as.matrix(read.table(infile,header=T,stringsAsFactors=F,sep="\t"))
ref_label = colnames(data2)[-1]
n_ref = length(ref_label)
ref = matrix(as.numeric(data2[,-1]),ncol=ncol(data2)-1)
gene = data2[,"gene"]
n_gene = length(gene)

infile =  file.path(PROJECT_DIR,"RESULTS","mut_driver","size_TCGA.txt")
data3 = as.matrix(read.table(infile,header=T,stringsAsFactors=F,sep="\t"))
ref_size = as.numeric(data3[,"size"])

hyper_pval = matrix(rep(NA,n_gene*n_ref),ncol=n_ref)
for (i_ref in 1:n_ref) {
    if (ref_size[i_ref]>n_MO) {
        for (i_gene in 1:n_gene) {
            # phyper(x, m, n, k) = probability of getting x or fewer white balls in a sample of
            # size k from an urn containing m white balls and n black balls
            m = round(ref_size[i_ref]*ref[i_gene,i_ref])
            n = ref_size[i_ref]-m
            x = round(n_MO*MO_freq[i_gene])
            k = n_MO
            hyper_pval[i_gene,i_ref] = 1-phyper(x, m, n, k)
        }
    } else if (ref_size[i_ref]<n_MO) {
        for (i_gene in 1:n_gene) {
            # phyper(x, m, n, k) = probability of getting x or fewer white balls in a sample of
            # size k from an urn containing m white balls and n black balls
            m = round(n_MO*MO_freq[i_gene])
            n = n_MO-m
            x = round(ref_size[i_ref]*ref[i_gene,i_ref])
            k = ref_size[i_ref]
            hyper_pval[i_gene,i_ref] = phyper(x, m, n, k)
        }
    }
}

output = rbind(c("gene",ref_label),cbind(gene,hyper_pval))
outfile =  file.path(PROJECT_DIR,"RESULTS","mut_driver","mut_pval_TCGA_byGene.txt")
write(t(output),ncol=ncol(output),file=outfile,sep="\t")

