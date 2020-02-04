library(maftools)
library(BSgenome)       ### Genomic data package for extracting sequence info

rm(list=ls())

PROJECT_DIR = "~/ACTIVE/LCS/Mongolia_HCC_Pipeline" # replace this line with your local path

destination = file.path(PROJECT_DIR,"RESULTS","mut_sig")
dir.create(destination)
setwd(destination)

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
    trinucl[[i_maf]] = trinucleotideMatrix(maf[[i_maf]],
    ref_genome="BSgenome.Hsapiens.UCSC.hg38")
}

tmp = apply(rbind(trinucl[[1]]$nmf_matrix,trinucl[[2]]$nmf_matrix),2,sum)
trinucl_freq_merged = tmp/sum(tmp)

trinucl_freq = vector("list",n_maf)
for (i_maf in 1:2) {
    tmp = apply(trinucl[[i_maf]]$nmf_matrix,2,sum)
    trinucl_freq[[i_maf]] = tmp/sum(tmp)
}
trinucl_freq_diff = (trinucl_freq[[1]]-trinucl_freq[[2]])/(trinucl_freq[[1]]+trinucl_freq[[2]])

pdf(file="mut_sig_trinucl_freq.pdf",width=7,height=5)
par(mfrow = c(2,1), oma = c(5,4,1,0) + 0.1, mar = c(0,0,2.5,0) + 0.1, las=1, tcl=-.25, font.main=4, xpd = NA)
color = c('coral4', 'lightcyan4', 'deeppink3', 'lightsalmon1', 'forestgreen', 'cornflowerblue')
colors = rep(color, each=16)

barplot(trinucl_freq_merged, xaxt="n", col=colors, beside=T, ylim = c(0,0.05),
cex.main = 1, border = NA, font.axis = 2, font.lab = 2, adj = 0.25, cex.axis=0.8, ylab="frequency (f)")

rect(xleft = seq(0, 6*19.2, 19.2), ybottom = 0.052, xright = 6*19.2, ytop = 0.054, col = color, border = 'gray70')
y = 0.065
text(labels = c("C>A","C>G","C>T","T>A","T>C","T>G"),x=19.2*(1:6)-19.2/2,y=rep(y,6),cex=1,font = 2, font.lab = 2, pos = 1.2, col="red")
text(labels = "pyr.",x=0,y=y,cex=0.7,font=2, font.lab=2, pos=1.2, col="red")
y = 0.072
text(labels = c("G>T","G>C","G>A","A>T","A>G","A>C"),x=19.2*(1:6)-19.2/2,y=rep(y,6),cex=1,font = 2, font.lab = 2, pos = 1.2, col="blue")
text(labels = "pur.",x=0,y=y,cex=0.7,font=2, font.lab=2, pos=1.2, col="blue")

barplot(trinucl_freq_diff, xaxt="n", col=colors, beside=T, ylim=c(-0.6,0.6),
cex.main = 1, border = NA, font.axis = 2, font.lab = 2, adj = 0.25, cex.axis=0.8, ylab="(fHDV+-fHDV-)/(fHDV++fHDV-)")

rect(xleft = seq(0, 6*19.2, 19.2), ybottom = -0.70, xright = 6*19.2, ytop = -0.65, col = color, border = 'gray70')

offset = -0.1
step = 19.2/16
epsilon = 0.1*step
x = (1:96)*step-step/2+epsilon
y = -0.675+offset
label=rep(c("T","G","C","A"),16)
text(labels=label, x=x, y=rep(y,96), cex=0.4, font=2, font.lab=2, pos=1.2, col="blue")
text(labels="preceded by 5'", x=-5.5, y=y, cex=0.4, font=2, font.lab=2, pos=1.2, col="blue")
y = -0.745+offset
label = rep(c(rep("T",4),rep("G",4),rep("C",4),rep("A",4)),6)
text(labels=label, x=x, y=rep(y,96), cex=0.4, font=2, font.lab=2, pos=1.2, col="blue")
text(labels="followed by 3'", x=-5.5, y=y, cex=0.4, font=2, font.lab=2, pos=1.2, col="blue")
y = -0.705+offset
text(labels = "pur.", x=-15, y=y, cex=0.7, font=2, font.lab=2, pos=1.2, col="blue")

y = -0.855+offset
label = rep(c(rep("A",4),rep("C",4),rep("G",4),rep("T",4)),6)
text(labels=label, x=x, y=rep(y,96), cex=0.4, font=2, font.lab=2, pos=1.2, col="red")
text(labels="preceded by 5'", x=-5.5, y=y, cex=0.4, font=2, font.lab=2, pos=1.2, col="red")
y = -0.925+offset
label=rep(c("A","C","G","T"),16)
text(labels=label, x=x, y=rep(y,96), cex=0.4, font=2, font.lab=2, pos=1.2, col="red")
text(labels="followed by 3'", x=-5.5, y=y, cex=0.4, font=2, font.lab=2, pos=1.2, col="red")
y = -0.885+offset
text(labels = "pyr.", x=-15, y=y, cex=0.7, font=2, font.lab=2, pos=1.2, col="red")

dev.off()
