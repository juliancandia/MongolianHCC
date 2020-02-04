rm(list=ls())

PROJECT_DIR = "~/ACTIVE/LCS/Mongolia_HCC_Pipeline" # replace this line with your local path

# We read in table to select driver genes.
infile = file.path(PROJECT_DIR,"RESULTS","mut_driver","genes_for_oncoplot.txt")
data = read.table(infile, sep="\t", stringsAsFactors=F, header=T)
gene = data[,"gene"]
n_gene = length(gene)

path_infile = file.path(PROJECT_DIR,"DATA","Driver_Genes")
file = list.files(path=path_infile)
n_file = length(file)
driver_genes = list("vector",n_file)
for (i_file in 1:n_file) {
    tmp = as.matrix(read.table(paste0(path_infile,"/",file[i_file]),header=F,stringsAsFactors=F,sep="\t"))[,1]
    driver_genes[[i_file]] = trimws(tmp[tmp!=""])
}

mut = matrix(rep(NA,n_gene*n_file),ncol=n_file)
for (i_gene in 1:n_gene) {
    for (i_file in 1:n_file) {
        mut[i_gene,i_file] = (gene[i_gene]%in%driver_genes[[i_file]])*1
    }
}
not_reported = rep(NA,n_gene)
not_reported[which(apply(mut,1,sum)==0)] = 1

outfile = file.path(PROJECT_DIR,"RESULTS","mut_driver","Table_driver_genes.txt")
output = rbind(c("Driver Gene",file,"not reported"),cbind(gene,mut,not_reported))
write(t(output),ncol=ncol(output),file=outfile,sep="\t")
