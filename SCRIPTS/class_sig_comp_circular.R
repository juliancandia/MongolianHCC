library(circlize)
require(gplots)

# Mongolia Data classified based on signatures from multiple studies

rm(list=ls())

PROJECT_DIR = "~/ACTIVE/LCS/Mongolia_HCC_Pipeline" # replace this line with your local path

setwd(file.path(PROJECT_DIR,"RESULTS","class_sig"))

dataset = c("Mongolia","TCGA","Hoshida","TIGER","Lee","Yamashita","Roessler")
dataset_label = c("Mongolia","TCGA","Hoshida","TIGER-LC","Lee","Yamashita","Roessler")
n_dataset = length(dataset)
class_label = vector("list",n_dataset)
class_label[[1]] = paste0("MO",1:4)
class_label[[2]] = paste0("iC",1:3)
class_label[[3]] = paste0("S",1:3)
#class_label[[4]] = paste0("HCC-C",1:3)
class_label[[4]] = paste0("C",1:3)
class_label[[5]] = c("UP","DN")
class_label[[6]] = c("UP","DN")
class_label[[7]] = c("UP","DN")

class_col = vector("list",n_dataset)
class_col[[1]] = c("#4DAF4A","#377EB8","orange","#E41A1C")
class_col[[2]] = c("magenta","#4DAF4A","#E41A1C")
class_col[[3]] = c("cyan","#E41A1C","#4DAF4A")
class_col[[4]] = c("#E41A1C","purple","#4DAF4A")
class_col[[5]] = c("#4DAF4A","#E41A1C")
class_col[[6]] = c("#E41A1C","#4DAF4A")
class_col[[7]] = c("#E41A1C","#4DAF4A")

# to filter class assignments based on quality thresholds
qc_A_dev_thres = 10
#pval_thres = 0.05
pval_thres = 1 # to avoid blanks
###########################

n_cl = rep(NA,n_dataset)
for (i_dataset in 1:n_dataset) {
    n_cl[i_dataset] = length(class_label[[i_dataset]])
}

file_A = file.path(PROJECT_DIR,"RESULTS","cluster","FINAL","patient_sample_metadata_w_clustering.txt")
data_A = as.matrix(read.table(file_A,header=T,check.names=F,stringsAsFactors=F,sep="\t"))
data_A = data_A[!is.na(data_A[,"RNASeq_ID"]),] # only samples with RNA-Seq
id_A = as.numeric(data_A[,"Patient"])
cl_A = data_A[,"class"]
qc_A = as.numeric(data_A[,"class_score"]) # 0=worst quality assignment; 1=best quality assignment.
qc_A_dev = abs(qc_A-median(qc_A))/mad(qc_A)

sample_id = id_A
n_sample = length(sample_id)
class = rep(NA,n_sample)
class[cl_A=="MO1"] = 1
class[cl_A=="MO2"] = 2
class[cl_A=="MO3"] = 3
class[cl_A=="MO4"] = 4
class[qc_A_dev>qc_A_dev_thres] = 0

#type = "cor"
type = "cos"
for (i_dataset in 2:n_dataset) {
    file_B = file.path(PROJECT_DIR,"DATA","Class_Signatures",paste0(dataset[i_dataset],"_",type,"_prediction_result.xls"))
    data_B = as.matrix(read.table(file_B,header=T,check.names=F,stringsAsFactors=F,sep="\t"))
    
    check = sum(data_B[,"sample.names"]==sample_id)==n_sample
    if (!check) {stop("Sample name mismatch!")}
    
    class_B = data_B[,"predict.label"]
    class_B[data_B[,"BH.FDR"]>pval_thres] = 0
    
    class = cbind(class,class_B)
}

# we reorder the samples based on Mongolia classes
o = order(class[,1])
class = class[o,]
sample_id = sample_id[o]

# we add a dummy group of sectors
n_sample_dummy = 23
sample_id = c(sample_id,paste0("d",1:n_sample_dummy))

pdf(paste0("class_sig_comp_circular_",type,".pdf"),width=5,height=5)

circos.clear()

circos.par(start.degree=90, clock.wise=F)
circos.par(cell.padding = c(0.02, 0.0, 0.0, 0.01))
circos.initialize(sample_id, xlim = c(0, 1))
for (i_dataset in 1:n_dataset) {
    circos.track(ylim = c(0, 1), track.height = 0.1)
}
circos.info(plot=F)

for (i_dataset in 1:n_dataset) {
    for (i_cl in 1:n_cl[i_dataset]) {
        index = which(class[,i_dataset]==i_cl)
        if (length(index)>0) {
            for (j in 1:length(index)) {
                ini = as.character(sample_id[index[j]])
                draw.sector(get.cell.meta.data("cell.start.degree", sector.index = ini),
                get.cell.meta.data("cell.end.degree", sector.index = ini),
                rou1 = get.cell.meta.data("cell.top.radius", track.index = i_dataset),
                rou2 = get.cell.meta.data("cell.bottom.radius", track.index = i_dataset),
                col = class_col[[i_dataset]][i_cl])
            }
        }
    }
}

dev.off()
