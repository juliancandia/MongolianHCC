library(gplots)

rm(list=ls())

PROJECT_DIR = "~/ACTIVE/LCS/Mongolia_HCC_Pipeline" # replace this line with your local path

destination = file.path(PROJECT_DIR,"RESULTS","class_sig")
dir.create(destination)
setwd(destination)

###########################
dataset = c("Hoshida","Lee","Roessler","TCGA","TIGER","Yamashita")
dataset_label = c("Hoshida","Lee","Roessler","TCGA","TIGER-LC","Yamashita")
n_dataset = length(dataset)
dt_col_label = vector("list",n_dataset)
dt_col_label[[1]] = paste0("S",1:3)
dt_col_label[[2]] = c("UP","DN")
dt_col_label[[3]] = c("UP","DN")
dt_col_label[[4]] = paste0("iC",1:3)
dt_col_label[[5]] = paste0("C",1:3)
dt_col_label[[6]] = c("UP","DN")

# to filter class assignments based on quality thresholds
qc_A_dev_thres = 10
#pval_thres = 0.05
pval_thres = 1
###########################

file_A = file.path(PROJECT_DIR,"RESULTS","cluster","FINAL","patient_sample_metadata_w_clustering.txt")
data_A = as.matrix(read.table(file_A,header=T,check.names=F,stringsAsFactors=F,sep="\t"))
data_A = data_A[!is.na(data_A[,"RNASeq_ID"]),] # only samples with RNA-Seq
id_A = as.numeric(data_A[,"Patient"])
cl_A = data_A[,"class"]
qc_A = as.numeric(data_A[,"class_score"]) # 0=worst quality assignment; 1=best quality assignment.
qc_A_dev = abs(qc_A-median(qc_A))/mad(qc_A)

outfile =  paste0("class_sig_comp_v2.pdf")
pdf(outfile,width=5,height=5)
for (i_dataset in 1:n_dataset) {
    file_B_label = c("cor","cos")
    file_B_label_full = c("NTP_correl","NTP_cosine")
    file_B = file.path(PROJECT_DIR,"DATA","Class_Signatures",paste0(dataset[i_dataset],"_",file_B_label,"_prediction_result.xls"))
    n_B = length(file_B)
    data_B = vector("list",n_B)
    id_B = vector("list",n_B)
    cl_B = vector("list",n_B)
    qc_B = vector("list",n_B)
    for (i_B in 1:n_B) {
        data_B[[i_B]] = as.matrix(read.table(file_B[i_B],header=T,check.names=F,stringsAsFactors=F,sep="\t"))
        id_B[[i_B]] = data_B[[i_B]][,"sample.names"]
        cl_B[[i_B]] = data_B[[i_B]][,"predict.label"]
        qc_B[[i_B]] = data_B[[i_B]][,"BH.FDR"]
    }
    
    # cor/cos consensus
    index_sel = which((qc_B[[1]]<pval_thres)&(qc_B[[2]]<pval_thres)&(cl_B[[1]]==cl_B[[2]])&(qc_A_dev<qc_A_dev_thres))
    if (length(index_sel)>10) {
        cl_A_consensus = cl_A[index_sel]
        cl_B_consensus = cl_B[[1]][index_sel]
        dt = table(cl_A_consensus,cl_B_consensus)
        colnames(dt) = dt_col_label[[i_dataset]]
        fisher.pval = fisher.test(dt)$p.value
        balloonplot(t(dt),main="",xlab="",ylab="",label=T,show.margins=F)
        title(paste0("Mongolia vs ",dataset_label[i_dataset]),line=2)
        #par(cex=0.85)
        title(paste0("Fisher's exact test p-value = ",signif(fisher.pval,digits=2)),line=0.5,cex.main=0.95)
        dt_consensus = dt
    }
}
dev.off()
