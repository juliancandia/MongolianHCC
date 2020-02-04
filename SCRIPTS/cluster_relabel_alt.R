rm(list=ls())

PROJECT_DIR = "~/ACTIVE/LCS/Mongolia_HCC_Pipeline" # replace this line with your local path

# original sample metadata file
infile = file.path(PROJECT_DIR,"DATA","ORIGINAL","patient_sample_metadata.txt")
metadata = read.table(infile,header=T,check.names=F,stringsAsFactors=F,sep="\t")
col_range = 1:ncol(metadata)

setwd(file.path(PROJECT_DIR,"RESULTS","cluster"))

# Here: we select an "alternative best solution" based on clustering stats
#best_thres = 1.25
best_thres = 1
best_K = 2

infile = file.path(paste0("T_mad",best_thres),"patient_sample_metadata_w_clustering.txt")
data = as.matrix(read.table(infile,header=T,check.names=F,stringsAsFactors=F,sep="\t"))

cl_ID = paste0("MO",c(4,3,1,2))
tmp = data[,paste0("K",best_K,"_class")]
for (i_cl in 1:length(cl_ID)) {
    tmp[as.numeric(tmp)==i_cl] = cl_ID[i_cl]
}
output = rbind(c(colnames(data)[1:ncol(metadata)],"class","class_score"),cbind(data[,1:ncol(metadata)],tmp,data[,paste0("K",best_K,"_ratio")]))
outfile = "FINAL/patient_sample_metadata_w_clustering_alt.txt"
write(t(output),ncol=ncol(output),file=outfile,sep="\t")
