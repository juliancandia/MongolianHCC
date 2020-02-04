library(dendextend)

rm(list=ls())

PROJECT_DIR = "~/ACTIVE/LCS/Mongolia_HCC_Pipeline" # replace this line with your local path

setwd(file.path(PROJECT_DIR,"RESULTS","cox"))

outfile = "risk_scores_dendro.pdf"
pdf(outfile,width=7,height=3.5)
data = as.matrix(read.table("risk_scores_wilcox.txt",header=T,check.names=F,stringsAsFactors=F,sep="\t"))
data = -log10(data)
#rownames(data) = colnames(data)
rownames(data) = rep("",ncol(data)) # to avoid labels
colnames(data) = rep("",ncol(data)) # to avoid labels
dend <- as.dist(data) %>% hclust("ave") %>% as.dendrogram
labels_colors(dend) = c("#4DAF4A","#377EB8","#ffe102","#E41A1C")
plot(dend)
dev.off()
