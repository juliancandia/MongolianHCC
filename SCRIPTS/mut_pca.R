library(maftools)
library(TCGAmutations)
library(BSgenome)
library(gplots)
library(calibrate)
library(FactoMineR)

rm(list=ls())

PROJECT_DIR = "~/ACTIVE/LCS/Mongolia_HCC_Pipeline" # replace this line with your local path

maf_file = file.path(PROJECT_DIR,"DATA","PROCESSED","maf_after_all_filters.maf")
mo = read.maf(maf_file)

maf_file = file.path(PROJECT_DIR,"DATA","PROCESSED","maf_after_all_filters_TCGA_LIHC_Asian.maf")
tcga.asian = read.maf(maf_file)

maf_file = file.path(PROJECT_DIR,"DATA","PROCESSED","maf_after_all_filters_TCGA_LIHC_Cauc.maf")
tcga.cauc = read.maf(maf_file)

maf_label = c("TCGA-Cauc","TCGA-Asian","Mongolia")
maf_color = c("green","blue","red")
maf_symbol = c(17,15,16)
n_maf = length(maf_label)
maf = vector("list",n_maf)
maf[[1]] = tcga.cauc
maf[[2]] = tcga.asian
maf[[3]] = mo

#### tri-nucleotide matrix for calculating mutational signatures
# NOTE: reference genomes must be installed individually from Bioconductor
trinucl = vector("list",n_maf)
for (i_maf in 1:2) {
    trinucl[[i_maf]] = trinucleotideMatrix(maf[[i_maf]], ref_genome="BSgenome.Hsapiens.UCSC.hg19", prefix="chr", add=T)$nmf_matrix
}
i_maf = 3
trinucl[[i_maf]] = trinucleotideMatrix(maf[[i_maf]], ref_genome="BSgenome.Hsapiens.UCSC.hg38")$nmf_matrix

mut_all = NULL
mut_sample_class = NULL
for (i_maf in 1:3) {
    mut_all = rbind(mut_all,trinucl[[i_maf]])
    mut_sample_class = c(mut_sample_class,rep(maf_label[i_maf],nrow(trinucl[[i_maf]])))
}

destination = file.path(PROJECT_DIR,"RESULTS","mut_pca")
dir.create(destination)
setwd(destination)

pdf(file.path(destination,"mut_pca.pdf"),width=7,height=7)

maf_legend = c("TCGA (Cauc)","TCGA (Asian)","Mongolia")
maf_color = c("#4DAF4A","blue","red")
maf_symbol = c(17,15,16)
maf_size = c(1.35,1.1,1)

scaled=T
data = mut_all
if (scaled) {
    for (i in 1:ncol(data)) {
        data[,i] = (data[,i]-mean(data[,i]))/sd(data[,i])
    }
}
# PCA plots
PCA_RES = prcomp(data)
pr.var = PCA_RES$sdev^2
pve=pr.var/sum(pr.var)
cum_pve=cumsum(pve)
x = PCA_RES$x[,1]
y = PCA_RES$x[,2]

# outliers
mad_thres = 10
sel_outlier = ((abs(x-median(x))/mad(x))>mad_thres)|((abs(y-median(y))/mad(y))>mad_thres)

# plot with outliers removed
x = x[!sel_outlier]
y = y[!sel_outlier]
mut_sample_class = mut_sample_class[!sel_outlier]
xrange = range(x)
yrange = range(y)
yrange[1] = yrange[1] - (yrange[2]-yrange[1])*0.05
plot(xrange,yrange,xaxt="n",yaxt="n",type="n",xlab="",ylab="")
mtext("PC1", 1, line=1, cex=1.2)
mtext("PC2", 2, line=1, cex=1.2)
maf_size = c(1.35,1.15,1.35)
###
aa = cbind.data.frame(factor(mut_sample_class),x,y)
bb = coord.ellipse(aa,bary=F,level.conf = 0.95) # level.conf = 0.95, npoint = 100, bary = FALSE
###
x_mean = rep(NA,n_maf)
y_mean = rep(NA,n_maf)
x_sd = rep(NA,n_maf)
y_sd = rep(NA,n_maf)
x_median = rep(NA,n_maf)
y_median = rep(NA,n_maf)
x_mad = rep(NA,n_maf)
y_mad = rep(NA,n_maf)
for (i_maf in 1:n_maf) {
    select = mut_sample_class==maf_label[i_maf]
    #col_par = as.numeric(col2rgb(maf_color[i_maf]))
    points(x[select],y[select],cex=maf_size[i_maf],col=maf_color[i_maf],pch=maf_symbol[i_maf])
    x_mean[i_maf] = mean(x[select])
    y_mean[i_maf] = mean(y[select])
    x_sd[i_maf] = sd(x[select])
    y_sd[i_maf] = sd(y[select])
    x_median[i_maf] = median(x[select])
    y_median[i_maf] = median(y[select])
    x_mad[i_maf] = mad(x[select])
    y_mad[i_maf] = mad(y[select])
}
for (i_maf in 1:n_maf) {
    sel = bb[[1]][,1]==maf_label[i_maf]
    #polygon(bb[[1]][sel,2:3],col=adjustcolor(maf_color[i_maf],alpha.f=0.25),border=maf_color[i_maf])
    polygon(bb[[1]][sel,2:3],col=NA,border=maf_color[i_maf],lwd=2)
    #points(x_median[i_maf],y_median[i_maf],cex=2,col=maf_color[i_maf],pch=4)
}
#x_legend = xrange[1] + (xrange[2]-xrange[1])*0.0
x_legend = 0.1
y_legend = yrange[1] + (yrange[2]-yrange[1])*1
legend(x_legend,y_legend,maf_legend,pch=maf_symbol,col=maf_color,cex=0.8,pt.cex=0.8)
dev.off()

# plot just barycenters and ellipses
pdf(file.path(destination,"mut_pca_inset.pdf"),width=7,height=7)
plot(range(bb[[1]][,2]),range(bb[[1]][,3]),xaxt="n",yaxt="n",type="n",xlab="",ylab="")
mtext("PC1", 1, line=1, cex=2.2)
mtext("PC2", 2, line=1, cex=2.2)
for (i_maf in 1:n_maf) {
    sel = bb[[1]][,1]==maf_label[i_maf]
    polygon(bb[[1]][sel,2:3],col=adjustcolor(maf_color[i_maf],alpha.f=0.05),border=maf_color[i_maf])
    polygon(bb[[1]][sel,2:3],col=NA,border=maf_color[i_maf],lwd=5)
    points(x_median[i_maf],y_median[i_maf],cex=4,col=maf_color[i_maf],pch=4,lwd=5)
}
dev.off()

# save PC1, PC2, PC3
# each COL is a PCA mode, each ROW is the contribution to that PCA mode from one analyte.
loading_vectors = PCA_RES$rotation
PC1_order=order(-abs(loading_vectors[,1]))
PC2_order=order(-abs(loading_vectors[,2]))
PC3_order=order(-abs(loading_vectors[,3]))
header = c("substitution_PC1","PC1","substitution_PC2","PC2","substitution_PC3","PC3")
substitution = colnames(data)
write(t(rbind(header,cbind(substitution[PC1_order],loading_vectors[PC1_order,1],substitution[PC2_order],loading_vectors[PC2_order,2],
substitution[PC3_order],loading_vectors[PC3_order,3]))),ncol=6,file=file.path(destination,"mut_pca.txt"),sep="\t")
