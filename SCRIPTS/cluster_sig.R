rm(list=ls())

PROJECT_DIR = "~/ACTIVE/LCS/Mongolia_HCC_Pipeline" # replace this line with your local path

setwd(file.path(PROJECT_DIR,"RESULTS","cluster","FINAL"))

best_thres = 2
best_K = 4

infile = file.path(PROJECT_DIR,"DATA","PROCESSED",paste0("T_mad",best_thres,"_expr.txt"))
expr = as.matrix(read.table(infile,header=F,check.names=F,stringsAsFactors=F,sep="\t"))
expr_Z = expr
expr_mean = apply(expr,1,mean)
expr_sd = apply(expr,1,sd)
for (i in 1:nrow(expr)) {
    expr_Z[i,] = (expr[i,]-expr_mean[i])/expr_sd[i] # if using Euclidean distance, gene z-score transformations matter
}

infile = file.path(PROJECT_DIR,"DATA","PROCESSED",paste0("T_mad",best_thres,"_genes.txt"))
gene = as.matrix(read.table(infile,header=F,check.names=F,stringsAsFactors=F,sep="\t"))
n_gene = nrow(gene)

infile = file.path(PROJECT_DIR,"DATA","PROCESSED","T_mad_patients.txt")
patient = as.matrix(read.table(infile,header=F,check.names=F,stringsAsFactors=F,sep="\t"))[,1]
n_pat = length(patient)

infile = "patient_sample_metadata_w_clustering.txt"
metadata = as.matrix(read.table(infile,header=T,check.names=F,stringsAsFactors=F,sep="\t"))
metadata = metadata[!is.na(metadata[,"RNASeq_ID"]),] # only patients with RNASeq data
check = (nrow(metadata)==n_pat)&(sum(trimws(metadata[,"Patient"])==patient)==n_pat)
if (!check) {
    stop("Patient mismatch!")
}

pval_thres = 0.001
class = names(table(metadata[,"class"]))
n_class = length(class)
gene_sig = rep(NA,n_gene)
gene_stats = matrix(rep(NA,n_gene*n_class*2),nrow=n_gene)
gene_stats_label = NULL
for (i_class in 1:n_class) {
    select = metadata[,"class"]==class[i_class]
    gene_stats[,2*(i_class-1)+1] = apply(expr_Z[,select],1,mean)
    gene_stats[,2*(i_class-1)+2] = apply(expr_Z[,select],1,sd)
    gene_stats_label = c(gene_stats_label,paste0(class[i_class],"_",c("mean","sd")))
    pval = rep(NA,n_gene)
    diff = rep(NA,n_gene)
    for (i_gene in 1:n_gene) {
        x = expr_Z[i_gene,select]
        y = expr_Z[i_gene,!select]
        pval[i_gene] = t.test(x,y)$p.value
        diff[i_gene] = median(x)-median(y)
    }
    select1 = (is.na(gene_sig))&(pval<pval_thres)
    select2 = (!is.na(gene_sig))&(pval<pval_thres)
    gene_sig[select1] = i_class*sign(diff[select1])
    gene_sig[select2] = "MC" # multi-cluster flag (to be later removed)
}
gene_sel = !is.na(as.numeric(gene_sig)) # we remove multi-cluster genes
gene_sig = abs(as.numeric(gene_sig[gene_sel])) # we remove the up/down-regulated distinction
gene_class = rep(NA,length(gene_sig))
for (i_class in 1:n_class) {
    gene_class[gene_sig==i_class] = class[i_class] # we restore class labels
}
expr_Z = expr_Z[gene_sel,]
gene = gene[gene_sel,]
output = rbind(c("gene","gene_cluster",paste0("PID_",patient)),cbind(gene[,3],gene_class,expr_Z))
outfile = "gene_expr_sig.txt"
write(t(output),ncol=ncol(output),file=outfile,sep="\t")

gene_stats = gene_stats[gene_sel,]
output = rbind(c("gene","gene_cluster",gene_stats_label),cbind(gene[,3],gene_class,gene_stats))
outfile = "gene_expr_sig_stats.txt"
write(t(output),ncol=ncol(output),file=outfile,sep="\t")
