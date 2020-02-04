library(gplots)
require(RColorBrewer)

rm(list=ls())

PROJECT_DIR = "~/ACTIVE/LCS/Mongolia_HCC_Pipeline" # replace this line with your local path

setwd(file.path(PROJECT_DIR,"RESULTS","mut_driver"))

freq = as.matrix(read.table("mut_freq_TCGA_byLocus.txt",header=T,stringsAsFactors=F,sep="\t"))
pval = as.matrix(read.table("mut_pval_TCGA_byLocus.txt",header=T,stringsAsFactors=F,sep="\t"))
ref_label = colnames(freq)[-(1:2)]
n_ref = length(ref_label)
gene = freq[,"gene"]
n_gene = length(gene)
locus = freq[,"locus"]
freq = matrix(as.numeric(freq[,-(1:2)]),ncol=ncol(freq)-2)
pval = matrix(as.numeric(pval[,-(1:2)]),ncol=ncol(pval)-2)

# We add Mongolia; remove "Unknown"
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
MO_pval = rep(1,n_gene) # by definition
freq[,which(ref_label=="Unknown")] = MO_freq
pval[,which(ref_label=="Unknown")] = MO_pval
ref_label[which(ref_label=="Unknown")] = "Mongolia"

outfile = "mut_heatmap_TCGA_byLocus.pdf"
pdf(outfile,width=7,height=5)

annot_mat = matrix(rep("",length(pval)),ncol=ncol(pval))
pval_thres = c(0.05)
pval_annot = c("*")
for (i in 1:length(pval_thres)) {
    index = which(pval<pval_thres[i],arr.ind=T)
    if (length(index)>0) {
        annot_mat[index] = pval_annot[i]
    }
}

heatmap_mat = 100*freq
heatmap_mat[which(heatmap_mat==0,arr.ind=T)] = NA # to enforce a different color to strict zero frequencies.
rownames(heatmap_mat) = paste0(gene," [",locus,"]")
colnames(heatmap_mat) = ref_label
myBreaks = 100*c(0,0.02,0.04,0.06,0.08,0.10,0.12)
myCol = colorRampPalette(brewer.pal(length(myBreaks)-1,"Blues"))
hm <- heatmap.2(heatmap_mat, scale="none", Rowv=F, Colv=F, na.rm=T, na.color="white",
col = myCol, breaks = myBreaks, dendrogram = "none", margins=c(7,7), cexRow=0.82, srtRow=0,srtCol=45, cexCol=0.67, key=T, trace="none", cellnote = annot_mat, notecol = "black", notecex = 1.3,
    density.info="none",key.title="",key.xlab="% Mutated samples",keysize=1.31)
#subtitle = "P-value significance codes:  <0.001 (***), <0.01 (**), <0.05 (*), <0.1 (.)"
#title(main="Mutation Frequency: Mongolia vs all TCGA Cancers", line=-1, cex.main=1)
#title(main=subtitle, line=-2, cex.main=0.55)

o = as.matrix(read.table("mut_heatmap_TCGA_byGene_order.txt",header=F,stringsAsFactors=F,sep="\t"))[,1]
hm <- heatmap.2(heatmap_mat[,o], scale="none", Rowv=F, Colv=F, na.rm=T, na.color="white",
col = myCol, breaks = myBreaks, dendrogram = "none", margins=c(7,7), cexRow=0.82, srtRow=0,srtCol=45, cexCol=0.67, key=F, trace="none", cellnote = annot_mat[,o], notecol = "black", notecex = 1.3)

dev.off()
