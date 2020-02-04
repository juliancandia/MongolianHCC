rm(list=ls())

PROJECT_DIR = "~/ACTIVE/LCS/Mongolia_HCC_Pipeline" # replace this line with your local path

infile = file.path(PROJECT_DIR,"DATA","PROCESSED","Hotspot_Loci_MANUAL.txt")
target = read.table(infile, sep="\t", stringsAsFactors=F, header=T)
n_target = nrow(target)

infile = file.path(PROJECT_DIR,"DATA","PROCESSED","hotspot.maf")
mo <- read.table(infile, sep="\t", quote="", header=T)
n_MO = 71 # total number of samples analyzed
MO_freq = rep(NA,n_target)
for (i_target in 1:n_target) {
    hit = grep(target[i_target,"Locus"],mo[mo[,"Hugo_Symbol"]==target[i_target,"Hugo_Symbol"],"HGVSp_Short"],value=T)
    MO_freq[i_target] = length(hit)/n_MO
}

infile =  file.path(PROJECT_DIR,"RESULTS","mut_driver","mut_freq_TCGA_byLocus.txt")
data2 = as.matrix(read.table(infile,header=T,stringsAsFactors=F,sep="\t"))
ref_label = colnames(data2)[-(1:2)]
n_ref = length(ref_label)
ref = matrix(as.numeric(data2[,-(1:2)]),ncol=ncol(data2)-2)
gene = data2[,"gene"]
n_gene = length(gene)
locus = data2[,"locus"]

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

output = rbind(c("gene","locus",ref_label),cbind(gene,locus,hyper_pval))
outfile =  file.path(PROJECT_DIR,"RESULTS","mut_driver","mut_pval_TCGA_byLocus.txt")
write(t(output),ncol=ncol(output),file=outfile,sep="\t")

