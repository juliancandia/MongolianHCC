rm(list=ls())

PROJECT_DIR = "~/ACTIVE/LCS/Mongolia_HCC_Pipeline" # replace this line with your local path

setwd(file.path(PROJECT_DIR,"DATA","PROCESSED"))

# we read in gene data (tumor and non-tumor in a single file)
data = as.matrix(read.table("gene_expr_voom.txt",header=T,check.names=F,stringsAsFactors=F,sep="\t"))
gene = cbind(data[,1],matrix(unlist(strsplit(data[,1],'|',fixed=T)),byrow=T,ncol=2)) # NOTE: transcript IDs are unique, gene IDs are not unique.
n_gene = nrow(gene)
sample = colnames(data)[-1]
n_sample = length(sample)
expr = matrix(as.numeric(data[,-1]),ncol=n_sample)

### we remove identical records
tmp = table(gene[,3])
multipl = names(tmp[tmp>=2])
for (i in 1:length(multipl)) {
    index = which(gene[,3]==multipl[i])
    keep = T
    for (j in 2:length(index)) {
        if (sum(expr[index[1],]==expr[index[j],])!=ncol(expr)) {
            keep = F
        }
    }
    if (keep) {
        expr = expr[-index[2:length(index)],]
        gene = gene[-index[2:length(index)],]
    }
}
n_gene = nrow(gene)

infile = file.path(PROJECT_DIR,"DATA","ORIGINAL","patient_sample_metadata.txt")
metadata = read.table(infile,header=T,check.names=F,stringsAsFactors=F,sep="\t")
metadata = metadata[!is.na(metadata[,"RNASeq_ID"]),] # only patients with RNASeq data

n_pat = nrow(metadata)
index_T = rep(NA,n_pat)
index_NT = rep(NA,n_pat)
for (i_pat in 1:n_pat) {
    index_T[i_pat] = which(sample==metadata[i_pat,"RNASeq_T"])
    index_NT[i_pat] = which(sample==metadata[i_pat,"RNASeq_NT"])
}

# T vs NT ratio
T_NT_ratio = matrix(rep(NA,n_gene*length(index_T)),nrow=n_gene)
for (i_gene in 1:n_gene) {
    T_NT_ratio[i_gene,] = expr[i_gene,index_T]-expr[i_gene,index_NT]
}
output = rbind(c("gene",metadata[,"LHC_ID"]),cbind(gene[,3],T_NT_ratio))
outfile = "gene_expr_voom_T_vs_NT_ratio.txt"
write(t(output),ncol=ncol(output),file=outfile,sep="\t")

# T vs NT analysis
ttest_T_NT_pval = rep(NA,n_gene)
FC_T_NT = rep(NA,n_gene)
for (i_gene in 1:n_gene) {
    ttest_T_NT_pval[i_gene] = t.test(expr[i_gene,index_T],expr[i_gene,index_NT],paired=T)$p.value
    FC_T_NT[i_gene] = median(2**(expr[i_gene,index_T]-expr[i_gene,index_NT]))
}
ttest_T_NT_pval_adj = p.adjust(ttest_T_NT_pval,method="fdr")
output = rbind(c("Ensembl|Hugo","Ensembl","Hugo","FC_T_NT","ttest_T_NT_pval","ttest_T_NT_pval_adj"),cbind(gene,FC_T_NT,ttest_T_NT_pval,ttest_T_NT_pval_adj))
outfile = "gene_expr_voom_T_vs_NT.txt"
write(t(output),ncol=ncol(output),file=outfile,sep="\t")

mad_T = apply(expr[,index_T],1,mad)

mad_T_thr = c(1,1.25,1.5,1.75,2,2.25,2.5,2.75,3)
gene_tot = NULL
for (thres in mad_T_thr) {
    select = (mad_T>thres)&(ttest_T_NT_pval_adj<0.05) # we apply filters
    gene_tot = c(gene_tot,sum(select))
    
    outfile = paste0("T_mad",thres,"_genes.txt")
    write(t(gene[select,]),ncol=ncol(gene),file=outfile,sep="\t")
    
    outfile = paste0("T_mad",thres,"_expr.txt")
    write(t(expr[select,index_T]),ncol=ncol(expr[,index_T]),file=outfile,sep="\t")
}
outfile = "T_mad_stats.txt"
output = rbind(c("mad_T_thr","gene_tot"),cbind(mad_T_thr,gene_tot))
write(t(output),ncol=ncol(output),file=outfile,sep="\t")

outfile = "T_mad_patients.txt"
write(metadata[,"Patient"],ncol=1,file=outfile)

# we write the full gene expression matrix for tumor samples (plain text)
output = rbind(c("NAME","Description",metadata[,"Patient"]),cbind(gene[,3],gene[,3],expr[,index_T]))
outfile = file.path(PROJECT_DIR,"DATA","PROCESSED","gene_expr_voom_T.txt")
write(t(output),ncol=ncol(output),file=outfile,sep="\t")

# also in gct format (after removing multiple entries associated with the same gene name)
# for downstream analysis in GenePattern
### we remove records with repeated genes
tmp = table(gene[,3])
multipl = names(tmp[tmp>=2])
for (i in 1:length(multipl)) {
    index = which(gene[,3]==multipl[i])
    expr = expr[-index,]
    gene = gene[-index,]
}
output = rbind(c("NAME","Description",metadata[,"Patient"]),cbind(gene[,3],gene[,3],expr[,index_T]))
outfile = "gene_expr_voom_T.gct"
write("#1.2",file=outfile,sep="\t")
write(c(nrow(expr[,index_T]),ncol(expr[,index_T])),ncol=2,file=outfile,sep="\t",append=T)
write(t(output),ncol=ncol(output),file=outfile,sep="\t",append=T)
