library(gplots)
require(RColorBrewer)

rm(list=ls())

PROJECT_DIR = "~/ACTIVE/LCS/Mongolia_HCC_Pipeline" # replace this line with your local path

setwd(file.path(PROJECT_DIR,"RESULTS","mut_driver"))

freq = as.matrix(read.table("mut_freq_TCGA_byGene.txt",header=T,stringsAsFactors=F,sep="\t"))
pval = as.matrix(read.table("mut_pval_TCGA_byGene.txt",header=T,stringsAsFactors=F,sep="\t"))
ref_label = colnames(freq)[-1]
n_ref = length(ref_label)
gene = freq[,"gene"]
n_gene = length(gene)
freq = matrix(as.numeric(freq[,-1]),ncol=ncol(freq)-1)
pval = matrix(as.numeric(pval[,-1]),ncol=ncol(pval)-1)

# We add Mongolia; remove "Unknown"
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
MO_pval = rep(1,n_gene) # by definition
unk_index = which(ref_label=="Unknown")
freq[,unk_index] = MO_freq
pval[,unk_index] = MO_pval
ref_label[unk_index] = "Mongolia"
ref_color = rep("black",n_ref)
ref_color[which(ref_label=="Mongolia")] = "red"

outfile = "mut_heatmap_TCGA_byGene.pdf"
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
rownames(heatmap_mat) = gene
colnames(heatmap_mat) = ref_label
myBreaks = 100*c(0,0.05,0.10,0.15,0.20,0.25,0.30,0.40)
myCol = colorRampPalette(brewer.pal(length(myBreaks)-1,"Blues"))
hm <- heatmap.2(heatmap_mat, scale="none", Rowv=F, Colv=T, na.rm=T, na.color="white",
col = myCol, breaks = myBreaks, dendrogram = "none", margins=c(7,7), cexRow=0.82, srtRow=0,srtCol=45, cexCol=0.67, key=T, trace="none", cellnote = annot_mat, notecol = "black", notecex = 1.3,
    density.info="none",key.title="",key.xlab="% Mutated samples",keysize=1.31)
#subtitle = "P-value significance codes:  <0.001 (***), <0.01 (**), <0.05 (*), <0.1 (.)"
#title(main="Mutation Frequency: Mongolia vs all TCGA Cancers", line=-1, cex.main=1)
#title(main=subtitle, line=-2, cex.main=0.55)

hm <- heatmap.2(heatmap_mat, scale="none", Rowv=F, Colv=T, na.rm=T, na.color="white",
col = myCol, breaks = myBreaks, dendrogram = "none", margins=c(7,7), cexRow=0.82, srtRow=0,srtCol=45, cexCol=0.67, key=F, trace="none", cellnote = annot_mat, notecol = "black", notecex = 1.3)

#heatmap_mat = 100*freq
#heatmap_mat[which(heatmap_mat==0,arr.ind=T)] = NA # to enforce a different color to strict zero frequencies.
#rownames(heatmap_mat) = gene
#hc = as.hclust(hm$colDendrogram)
#hm <- heatmap.2(heatmap_mat, scale="none", Rowv=F, Colv=T, na.rm=T, na.color="white",
#col = myCol, breaks = myBreaks, dendrogram = "none", margins=c(7,7), cexRow=0.82, srtRow=0,srtCol=45, cexCol=0.67, key=F, trace="none", cellnote = annot_mat, notecol = "black", notecex = 1.3,
#add.expr=mtext(side=1, text=ref_label[hc$order], at=1:n_ref, las=2, line = 0.5,col = 1:ncol(module2), cex = 0.8))
#)

#text(par("usr")[2] + 0.25, 5.5, srt=-45, adj = 0, labels = "dependent",
#xpd = TRUE)

# hm$colDendrogram -> extract order of columns to import into Locus heatmap
hc = as.hclust(hm$colDendrogram) # extract order of columns to import into Locus heatmap
write(hc$order, ncol=1, file="mut_heatmap_TCGA_byGene_order.txt")


dev.off()
