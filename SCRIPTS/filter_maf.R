
rm(list=ls())

PROJECT_DIR = "~/ACTIVE/LCS/Mongolia_HCC_Pipeline" # replace this line with your local path

infile = file.path(PROJECT_DIR,"DATA","ORIGINAL","mutect2_merged.maf")
maf_data <- read.table(infile, sep="\t", quote="", header=T) # 194682 alterations
n_sample = length(unique(maf_data[,"Tumor_Sample_Barcode"])) # 71 samples

# Read list of blacklisted genes
flags <- c("TTN","MUC16","OBSCN","AHNAK2","SYNE1","FLG","MUC5B","DNAH17","PLEC","DST","SYNE2","NEB","HSPG2","LAMA5","AHNAK","HMCN1","USH2A","DNAH11","MACF1","MUC17","DNAH5","GPR98","FAT1","PKD1","MDN1","RNF213","RYR1","DNAH2","DNAH3","DNAH8","DNAH1","DNAH9","ABCA13","SRRM2","CUBN","SPTBN5","PKHD1","LRP2","FBN3","CDH23","DNAH10","FAT4","RYR3","PKHD1L1","FAT2","CSMD1","PCNT","COL6A3","FRAS1","FCGBP","RYR2","HYDIN","XIRP2","LAMA1")
sel_not_flagged = !maf_data[,"Hugo_Symbol"]%in%flags # 97.4%

# Filters based on normal allele frequencies
tmp = as.numeric(as.character(maf_data[,"ExAC_AF"]))
sel_ExAC_AF = rep(T,length(tmp))
sel_ExAC_AF[which(tmp>0.001)] = F # 99.9%
tmp = as.numeric(as.character(maf_data[,"gnomAD_AF"]))
sel_gnomAD_AF = rep(T,length(tmp))
sel_gnomAD_AF[which(tmp>0.001)] = F # 99.8%
tmp = as.numeric(as.character(maf_data[,"AF"]))
sel_AF01 = rep(T,length(tmp))
sel_AF01[which(tmp>0.01)] = F # 99.8%
tmp = as.numeric(as.character(maf_data[,"AF"]))
sel_AF001 = rep(T,length(tmp))
sel_AF001[which(tmp>0.001)] = F # 99.4%

# Filter based on variant type (MAFTOOLS hard-coded list of nonsilent mutations)
sel_var = maf_data[,"Variant_Classification"]%in%c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site",
"Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del",
"In_Frame_Ins", "Missense_Mutation") # 53.6%

# Filters based on read counts in normal and tumor
sel_tdepth20 = maf_data[,"t_depth"]>20 # 99.9%
frac_nalt = maf_data[,"n_alt_count"]/maf_data[,"n_depth"]
sel_frac_nalt_1 = frac_nalt<0.01 # 85.2%
sel_frac_nalt_2 = frac_nalt<0.02 # 95.9%
frac_talt = maf_data[,"t_alt_count"]/maf_data[,"t_depth"] # This is exactly "tumor_freq" (last column)
sel_frac_talt_5 = frac_talt>0.05 # 10.4%

# All filters combined (NOTE: we don't apply filters on normal frequencies).
select = sel_not_flagged&sel_ExAC_AF&sel_gnomAD_AF&sel_AF001&sel_var&sel_tdepth20&sel_frac_talt_5 # 4.19% = 8149 variants

# We filter the maf and extract gene lists based on the fraction of mutated samples in the cohort
maf_data_filt = maf_data[select,]
target = unique(maf_data_filt[,"Hugo_Symbol"])
n_target = length(target) # 5797 genes
mut_frac = rep(NA,n_target)
names(mut_frac) = target
for (i_target in 1:n_target) {
    mut_frac[i_target] = length(unique(maf_data_filt[maf_data_filt[,"Hugo_Symbol"]==target[i_target],
    "Tumor_Sample_Barcode"]))/n_sample
}
# Genes above 5% and 10% thresholds:
mut_genes_5perc = names(mut_frac[mut_frac>(floor(n_sample*0.05)/n_sample)]) # freq>5%
mut_genes_10perc = names(mut_frac[mut_frac>(floor(n_sample*0.10)/n_sample)]) # freq>10%

# mutsigcv generated from the full unfiltered maf
infile = file.path(PROJECT_DIR,"DATA","ORIGINAL","mut_full.sig_genes.txt")
mutsigcv = read.table(infile, sep="\t", quote="", header=T) # 18862 genes
mutsigcv_filt = mutsigcv[mutsigcv[,"gene"]%in%maf_data_filt[,"Hugo_Symbol"],] # 5197 genes
#mutsigcv_filt_p_adj = p.adjust(mutsigcv_filt[,"p"],method="fdr") # adjust p-values based on this list.
mutsigcv_filt_p_adj = mutsigcv_filt[,"q"] # NOTE: here we don't adjust the p-values, we still use the q values reported by mutsigcv.
select = (mutsigcv_filt[,"gene"]%in%mut_genes_5perc)&(mutsigcv_filt_p_adj<0.1)
mutsigcv_filt_q10 = mutsigcv_filt[select,] # 10 genes

mut_genes_10perc_add = mut_genes_10perc[!mut_genes_10perc%in%mutsigcv_filt_q10[,"gene"]] # 9 additional genes

# genes from top panel
top = mut_frac[names(mut_frac)%in%mutsigcv_filt_q10[,"gene"]]
top = top[order(-top)]
# genes from bottom panel
bot = mut_frac[names(mut_frac)%in%mut_genes_10perc_add]
bot = bot[order(-bot)]

# we save results
destination = file.path(PROJECT_DIR,"RESULTS","mut_driver")
dir.create(destination)
setwd(destination)
output = rbind(c("gene","freq","panel"),cbind(names(top),top,rep("top",length(top)))
,cbind(names(bot),bot,rep("bottom",length(bot))))
outfile = file.path(destination,"genes_for_oncoplot.txt")
write(t(output),ncol=ncol(output),file=outfile,sep="\t")

# we save the filtered maf
output = maf_data_filt
outfile = file.path(PROJECT_DIR,"DATA","PROCESSED","maf_after_all_filters.maf")
write(colnames(output),ncol=ncol(output),file=outfile,sep="\t")
write(t(output),ncol=ncol(output),file=outfile,sep="\t",append=T)

# maf with just the final genes in oncoplot
output = maf_data_filt[maf_data_filt[,"Hugo_Symbol"]%in%c(names(top),names(bot)),]
outfile = file.path(PROJECT_DIR,"DATA","PROCESSED","maf_only_genes_in_oncoplot.maf")
write(colnames(output),ncol=ncol(output),file=outfile,sep="\t")
write(t(output),ncol=ncol(output),file=outfile,sep="\t",append=T)

# we generate separate mafs based on HDV status
infile = file.path(PROJECT_DIR,"DATA","ORIGINAL","patient_sample_metadata.txt")
sample_metadata = as.matrix(read.table(infile,sep="\t",stringsAsFactors=F,header=T))
sample_metadata = sample_metadata[!is.na(sample_metadata[,"WES_T"]),]
hdv_pos = sample_metadata[which(sample_metadata[,"hdv"]==1),"WES_T"]
hdv_neg = sample_metadata[which(sample_metadata[,"hdv"]==0),"WES_T"]
maf_pos = maf_data_filt[maf_data_filt[,"Tumor_Sample_Barcode"]%in%hdv_pos,]
maf_neg = maf_data_filt[maf_data_filt[,"Tumor_Sample_Barcode"]%in%hdv_neg,]

# we save the hdv split mafs
output = maf_pos
outfile = file.path(PROJECT_DIR,"DATA","PROCESSED","maf_after_all_filters_HDV_pos.maf")
write(colnames(output),ncol=ncol(output),file=outfile,sep="\t")
write(t(output),ncol=ncol(output),file=outfile,sep="\t",append=T)

output = maf_neg
outfile = file.path(PROJECT_DIR,"DATA","PROCESSED","maf_after_all_filters_HDV_neg.maf")
write(colnames(output),ncol=ncol(output),file=outfile,sep="\t")
write(t(output),ncol=ncol(output),file=outfile,sep="\t",append=T)

