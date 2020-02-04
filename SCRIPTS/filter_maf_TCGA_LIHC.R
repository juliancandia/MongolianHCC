library(maftools)       ### For read/writing/plotting/manipulating maf file
library(TCGAmutations)  ### Basically an interface to a Github repo with pre-packaged maf objects for TCGA data (https://github.com/PoisonAlien/TCGAmutations)

rm(list=ls())

PROJECT_DIR = "~/ACTIVE/LCS/Mongolia_HCC_Pipeline" # replace this line with your local path

tcga_load(study = "LIHC") # loads maf object named 'tcga_lihc_mc3'
maf_data = as.matrix(tcga_lihc_mc3@data)
n_sample = length(unique(maf_data[,"Tumor_Sample_Barcode"]))

# Read list of blacklisted genes
flags <- c("TTN","MUC16","OBSCN","AHNAK2","SYNE1","FLG","MUC5B","DNAH17","PLEC","DST","SYNE2","NEB","HSPG2","LAMA5","AHNAK","HMCN1","USH2A","DNAH11","MACF1","MUC17","DNAH5","GPR98","FAT1","PKD1","MDN1","RNF213","RYR1","DNAH2","DNAH3","DNAH8","DNAH1","DNAH9","ABCA13","SRRM2","CUBN","SPTBN5","PKHD1","LRP2","FBN3","CDH23","DNAH10","FAT4","RYR3","PKHD1L1","FAT2","CSMD1","PCNT","COL6A3","FRAS1","FCGBP","RYR2","HYDIN","XIRP2","LAMA1")
sel_not_flagged = !maf_data[,"Hugo_Symbol"]%in%flags

# Filters based on normal allele frequencies
tmp = as.numeric(as.character(maf_data[,"ExAC_AF"]))
sel_ExAC_AF = rep(T,length(tmp))
sel_ExAC_AF[which(tmp>0.001)] = F

# Filter based on variant type (MAFTOOLS hard-coded list of nonsilent mutations)
sel_var = maf_data[,"Variant_Classification"]%in%c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site",
"Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del",
"In_Frame_Ins", "Missense_Mutation")

# Filters based on read counts in tumor
sel_tdepth20 = as.numeric(maf_data[,"t_depth"])>20
frac_talt = as.numeric(maf_data[,"t_alt_count"])/as.numeric(maf_data[,"t_depth"])
sel_frac_talt_5 = frac_talt>0.05

# All filters combined (NOTE: we don't apply filters on normal frequencies).
select = sel_not_flagged&sel_ExAC_AF&sel_var&sel_tdepth20&sel_frac_talt_5

maf = maf_data[select,]
# we save the filtered maf
output = maf
outfile = file.path(PROJECT_DIR,"DATA","PROCESSED","maf_after_all_filters_TCGA_LIHC.maf")
write(colnames(output),ncol=ncol(output),file=outfile,sep="\t")
write(t(output),ncol=ncol(output),file=outfile,sep="\t",append=T)

# we save ethnicity-based subcohorts.
tcga_clindata <- getClinicalData(tcga_lihc_mc3)
asian <- levels(droplevels(tcga_clindata$Tumor_Sample_Barcode[tcga_clindata$race_list=="ASIAN"]))
output = maf[maf[,"Tumor_Sample_Barcode"]%in%asian,]
outfile = file.path(PROJECT_DIR,"DATA","PROCESSED","maf_after_all_filters_TCGA_LIHC_Asian.maf")
write(colnames(output),ncol=ncol(output),file=outfile,sep="\t")
write(t(output),ncol=ncol(output),file=outfile,sep="\t",append=T)

cauc <- levels(droplevels(tcga_clindata$Tumor_Sample_Barcode[tcga_clindata$race_list=="WHITE"]))
output = maf[maf[,"Tumor_Sample_Barcode"]%in%cauc,]
outfile = file.path(PROJECT_DIR,"DATA","PROCESSED","maf_after_all_filters_TCGA_LIHC_Cauc.maf")
write(colnames(output),ncol=ncol(output),file=outfile,sep="\t")
write(t(output),ncol=ncol(output),file=outfile,sep="\t",append=T)
