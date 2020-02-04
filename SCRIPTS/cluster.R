library(ConsensusClusterPlus)

rm(list=ls())

PROJECT_DIR = "~/ACTIVE/LCS/Mongolia_HCC_Pipeline" # replace this line with your local path

destination = file.path(PROJECT_DIR,"RESULTS","cluster")
dir.create(destination)
setwd(destination)

mad_T_thr = c(1,1.25,1.5,1.75,2,2.25,2.5,2.75,3)
n_thres = length(mad_T_thr)

seed=111
maxK=8
reps=1000
pItem=0.8
pFeature=0.8
clusterAlg="km"
distance="euclidean"
innerLinkage="average"
finalLinkage="average"

infile = file.path(PROJECT_DIR,"DATA","PROCESSED","T_mad_patients.txt")
patient = as.matrix(read.table(infile,header=F,check.names=F,stringsAsFactors=F,sep="\t"))[,1]
n_pat = length(patient)

infile = file.path(PROJECT_DIR,"DATA","ORIGINAL","patient_sample_metadata.txt")
metadata = read.table(infile,header=T,check.names=F,stringsAsFactors=F,sep="\t")
index = rep(NA,n_pat)
for (i_pat in 1:n_pat) {
    index[i_pat] = which(metadata[,"Patient"]==patient[i_pat])
}

for (i_thres in 1:n_thres) {
    thres = mad_T_thr[i_thres]
    infile = file.path(PROJECT_DIR,"DATA","PROCESSED",paste0("T_mad",thres,"_expr.txt"))
    expr = as.matrix(read.table(infile,header=F,check.names=F,stringsAsFactors=F,sep="\t"))
    expr_Z = expr
    expr_mean = apply(expr,1,mean)
    expr_sd = apply(expr,1,sd)
    for (i in 1:nrow(expr)) {
        expr_Z[i,] = (expr[i,]-expr_mean[i])/expr_sd[i] # if using Euclidean distance, gene z-score transformations matter
    }
    colnames(expr_Z) = patient
    
    dir.create(paste0("T_mad",thres))
    cClust_Z = ConsensusClusterPlus(d=expr_Z,maxK=maxK,reps=reps,pItem=pItem,
    pFeature=pFeature, clusterAlg=clusterAlg,distance=distance, plot="pdf", title=paste0("T_mad",thres),
    seed=seed, writeTable=T, innerLinkage=innerLinkage, finalLinkage=finalLinkage)
    
    patient_summary = NULL
    patient_summary_desc = NULL
    for (i_clust in 2:maxK) {
        class = cClust_Z[[i_clust]]$consensusClass # cluster labels (integers)
        cMat = cClust_Z[[i_clust]]$consensusMatrix
        cIn_mean = rep(NA,n_pat)
        cOut_mean = rep(NA,n_pat)
        for (i_pat in 1:n_pat) {
            class_in = class==class[i_pat]
            class_out = !class_in
            class_in[i_pat] = F # to avoid self
            cIn_this_pat = cMat[class_in,i_pat]
            cOut_this_pat = cMat[class_out,i_pat]
            cIn_mean[i_pat] = mean(cIn_this_pat)
            cOut_mean[i_pat] = mean(cOut_this_pat)
        }
        ratio = cIn_mean/(cIn_mean+cOut_mean)
        patient_summary = cbind(patient_summary,class,ratio)
        patient_summary_desc = c(patient_summary_desc,paste0("K",i_clust,"_",c("class","ratio")))
    }
    
    metadata2 = matrix(rep(NA,nrow(metadata)*ncol(patient_summary)),ncol=ncol(patient_summary))
    metadata2[index,] = patient_summary
    colnames(metadata2) = patient_summary_desc
    metadata2 = cbind(metadata,metadata2)
    # Patient summary info
    outfile = file.path(paste0("T_mad",thres),"patient_sample_metadata_w_clustering.txt")
    output = rbind(colnames(metadata2),metadata2)
    write(t(output),ncol=ncol(output),file=outfile,sep="\t")
}
