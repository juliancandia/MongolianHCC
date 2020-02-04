library(maftools)       ### For read/writing/plotting/manipulating maf file
library(TCGAmutations)  ### Basically an interface to a Github repo with pre-packaged maf objects for TCGA data (https://github.com/PoisonAlien/TCGAmutations)

rm(list=ls())

PROJECT_DIR = "~/ACTIVE/LCS/Mongolia_HCC_Pipeline" # replace this line with your local path

infile = file.path(PROJECT_DIR,"DATA","PROCESSED","maf_only_genes_in_oncoplot.maf")
mo <- read.table(infile, sep="\t", quote="", header=T)

# analyze only candidate "Mongolia genes"
target = c("TP53","CTNNB1","ALB","GTF2IRD2B","RB1","COL11A1","PNRC2","AK2","VPS13A","SPTA1",            "PCLO","CSMD2","APOB","SMC6","CSMD3","DYNC2H1","FKBP9","LRP1B","PCDH7") # same order as in the oncoplot.
n_target = length(target)

study = tcga_available()[[1]]
#study = study[-length(study)] # we remove "unknown"
n_study = length(study)
for (i_study in 1:n_study) {
    tcga_load(study=study[i_study])
}

maf = vector("list",n_study)
maf[[1]] = tcga_acc_mc3
maf[[2]] = tcga_blca_mc3
maf[[3]] = tcga_brca_mc3
maf[[4]] = tcga_cesc_mc3
maf[[5]] = tcga_chol_mc3
maf[[6]] = tcga_coad_mc3
maf[[7]] = tcga_dlbc_mc3
maf[[8]] = tcga_esca_mc3
maf[[9]] = tcga_gbm_mc3
maf[[10]] = tcga_hnsc_mc3
maf[[11]] = tcga_kich_mc3
maf[[12]] = tcga_kirc_mc3
maf[[13]] = tcga_kirp_mc3
maf[[14]] = tcga_laml_mc3
maf[[15]] = tcga_lgg_mc3
maf[[16]] = tcga_lihc_mc3
maf[[17]] = tcga_luad_mc3
maf[[18]] = tcga_lusc_mc3
maf[[19]] = tcga_meso_mc3
maf[[20]] = tcga_ov_mc3
maf[[21]] = tcga_paad_mc3
maf[[22]] = tcga_pcpg_mc3
maf[[23]] = tcga_prad_mc3
maf[[24]] = tcga_read_mc3
maf[[25]] = tcga_sarc_mc3
maf[[26]] = tcga_skcm_mc3
maf[[27]] = tcga_stad_mc3
maf[[28]] = tcga_tgct_mc3
maf[[29]] = tcga_thca_mc3
maf[[30]] = tcga_thym_mc3
maf[[31]] = tcga_ucec_mc3
maf[[32]] = tcga_ucs_mc3
maf[[33]] = tcga_uvm_mc3
maf[[34]] = tcga_unknown

n_sample = rep(NA,n_study)
for (i_study in 1:n_study) {
    maf_data = as.matrix(maf[[i_study]]@data)
    n_sample[i_study] = length(unique(maf_data[,"Tumor_Sample_Barcode"]))
    
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
    
    # Genes in target list
    sel_gene_in_list = maf_data[,"Hugo_Symbol"]%in%target
    
    # All filters combined (NOTE: we don't apply filters on normal frequencies).
    select = sel_not_flagged&sel_ExAC_AF&sel_var&sel_tdepth20&sel_frac_talt_5&sel_gene_in_list
    
    maf[[i_study]] = maf_data[select,]
}

mut_freq = matrix(rep(0,n_target*n_study),ncol=n_study)
for (i_study in 1:n_study) {
    tmp = table(maf[[i_study]][,"Hugo_Symbol"],maf[[i_study]][,"Tumor_Sample_Barcode"]) # integer matrix genes (rows) x subjects (columns). Only mutated subjects appear listed.
    freq = apply(tmp>0,1,sum)/n_sample[i_study]
    gene = names(freq)
    for (i_target in 1:n_target) {
        index = which(gene==target[i_target])
        if (length(index)>0) {
            mut_freq[i_target,i_study] = freq[index]
        }
    }
}

output = rbind(c("gene",study),cbind(target,mut_freq))
outfile =  file.path(PROJECT_DIR,"RESULTS","mut_driver","mut_freq_TCGA_byGene.txt")
write(t(output),ncol=ncol(output),file=outfile,sep="\t")

output = rbind(c("study","size"),cbind(study,n_sample))
outfile =  file.path(PROJECT_DIR,"RESULTS","mut_driver","size_TCGA.txt")
write(t(output),ncol=ncol(output),file=outfile,sep="\t")
