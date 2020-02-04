library(maftools)
library(BSgenome)       ### Genomic data package for extracting sequence info
library(deconstructSigs)

rm(list=ls())

PROJECT_DIR = "~/ACTIVE/LCS/Mongolia_HCC_Pipeline" # replace this line with your local path

tri.counts.method = "default"
#tri.counts.method = "exome"
#tri.counts.method = "exome2genome"

maf_file = file.path(PROJECT_DIR,"DATA","PROCESSED","maf_after_all_filters_HDV_pos.maf")
mo_hdv_pos = read.maf(maf_file)
maf_file = file.path(PROJECT_DIR,"DATA","PROCESSED","maf_after_all_filters_HDV_neg.maf")
mo_hdv_neg = read.maf(maf_file)

maf_label = c("HDV_pos","HDV_neg")
n_maf = length(maf_label)
maf = vector("list",n_maf)
maf[[1]] = mo_hdv_pos
maf[[2]] = mo_hdv_neg

#### tri-nucleotide matrix for calculating mutational signatures
# NOTE: reference genomes must be installed individually from Bioconductor
trinucl = vector("list",n_maf)
for (i_maf in 1:2) {
    trinucl[[i_maf]] = as.data.frame(trinucleotideMatrix(maf[[i_maf]],
    ref_genome="BSgenome.Hsapiens.UCSC.hg38")$nmf_matrix)
}

signatures = c("Reference","Signature_Index","Signature_Label","p-value")
prevalence = "Prevalence"
sample_weights = NULL
for (sig_ref in c("COSMIC_v2","COSMIC_v3","EnvirAg")) {
    infile = file.path(PROJECT_DIR,"DATA","Mutational_Signatures",paste0(sig_ref,"_signature.txt"))
    data = as.matrix(read.table(infile,sep="\t",header=T))
    substitution_labels = data[,3]
    data = data[,-(1:3)]
    sig_labels = colnames(data)
    n_sig = length(sig_labels)
    sig_ref_data = t(matrix(as.numeric(data),ncol=ncol(data)))
    rownames(sig_ref_data) = sig_labels
    colnames(sig_ref_data) = substitution_labels
    sig_ref_data = as.data.frame(sig_ref_data)

    weights.mat = vector("list",n_maf)
    for (i_maf in 1:n_maf) {
        sample = rownames(trinucl[[i_maf]])
        n_sample = length(sample)
        weights.mat[[i_maf]] = matrix(rep(NA,n_sample*n_sig),ncol=n_sig)
        rownames(weights.mat[[i_maf]]) = sample
        colnames(weights.mat[[i_maf]]) = rownames(sig_ref_data)
        for (i_sample in 1:n_sample) {
            weights.mat[[i_maf]][i_sample,] = as.numeric(whichSignatures(tumor.ref = trinucl[[i_maf]],
            signatures.ref = sig_ref_data,
            sample.id = sample[i_sample], contexts.needed = T, tri.counts.method = tri.counts.method)$weights)
        }
    }

    wilcox = rep(NA,n_sig)
    for (i_sig in 1:n_sig) {
        wilcox[i_sig] = wilcox.test(weights.mat[[1]][,i_sig],weights.mat[[2]][,i_sig])$p.value
    }
    sig_index = which(wilcox<0.05)
    
    if (length(sig_index)>1) {
        signatures = rbind(signatures,cbind(rep(sig_ref,length(sig_index)),sig_index,rep("",length(sig_index)),wilcox[sig_index]))
        prev = rep(maf_label[2],length(sig_index))
        prev[apply(weights.mat[[1]][,sig_index],2,mean)>apply(weights.mat[[2]][,sig_index],2,mean)] = maf_label[1]
        prevalence = c(prevalence,prev)
        sample_weights = rbind(sample_weights,cbind(t(weights.mat[[1]][,sig_index]),t(weights.mat[[2]][,sig_index])))
    } else if (length(sig_index)==1) {
        signatures = rbind(signatures,c(sig_ref,sig_index,"",wilcox[sig_index]))
        if (mean(weights.mat[[1]][,sig_index])>mean(weights.mat[[2]][,sig_index])) {
            prevalence = c(prevalence,maf_label[1])
        } else {
            prevalence = c(prevalence,maf_label[2])
        }
        sample_weights = rbind(sample_weights,c(weights.mat[[1]][,sig_index],weights.mat[[2]][,sig_index]))
    }
}
sample_weights = rbind(c(rownames(weights.mat[[1]]),rownames(weights.mat[[2]])),sample_weights)
output = cbind(signatures,prevalence,sample_weights)
outfile = paste0(file.path(PROJECT_DIR,"RESULTS","mut_sig"),"/weights_HDV_",tri.counts.method,".txt")
write(t(output),ncol=ncol(output),file=outfile,sep="\t")
