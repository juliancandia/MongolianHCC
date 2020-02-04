rm(list=ls())

PROJECT_DIR = "~/ACTIVE/LCS/Mongolia_HCC_Pipeline" # replace this line with your local path

# To explore the existence of hotspots
infile = file.path(PROJECT_DIR,"DATA","PROCESSED","maf_only_genes_in_oncoplot.maf")
maf_data <- read.table(infile, sep="\t", quote="", header=T) # 146 alterations
n_sample = length(unique(maf_data[,"Tumor_Sample_Barcode"])) # 61/71 samples with mutations
var = paste0(maf_data[,"Hugo_Symbol"],"|",maf_data[,"Variant_Classification"],"|",maf_data[,"Start_Position"])
var_t = table(var)
var_t = var_t[order(-var_t)]
var_t = var_t[var_t>1] # threshold on variant hotspots

# we save the aggregate results
var_desc = c("gene","type","location")
var_split = matrix(unlist(strsplit(names(var_t),"|",fixed=T)),byrow=T,ncol=length(var_desc))
colnames(var_split) = var_desc
output = rbind(c(var_desc,"frequency"),cbind(var_split,var_t))
outfile = file.path(PROJECT_DIR,"DATA","PROCESSED","hotspot_stats.txt")
write(t(output),ncol=ncol(output),file=outfile,sep="\t")

index = NULL
for (i in 1:length(var_t)) {
    index = c(index, which((maf_data[,"Hugo_Symbol"]==var_split[i,"gene"])&(maf_data[,"Variant_Classification"]==var_split[i,"type"])&(maf_data[,"Start_Position"]==var_split[i,"location"])))
}
# we save the maf subset
output = maf_data[index,]
outfile = file.path(PROJECT_DIR,"DATA","PROCESSED","hotspot.maf")
write(colnames(output),ncol=ncol(output),file=outfile,sep="\t")
write(t(output),ncol=ncol(output),file=outfile,sep="\t",append=T)
