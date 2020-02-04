require(limma) # voom normalization

rm(list=ls())

PROJECT_DIR = "~/ACTIVE/LCS/Mongolia_HCC_Pipeline" # replace this line with your local path

destination = file.path(PROJECT_DIR,"DATA","PROCESSED")
dir.create(destination)
setwd(destination)

infile = file.path(PROJECT_DIR,"DATA","ORIGINAL","patient_sample_metadata.txt")
sample_info = read.table(infile,header=T,check.names=F,stringsAsFactors=F,sep="\t")
infile = file.path(PROJECT_DIR,"DATA","ORIGINAL","RawCountFile_RSEM_genes_filtered.txt")
expr.raw = read.table(infile,header=T,check.names=F,stringsAsFactors=F,sep="\t")

# we save the list of unique gene names in the assay
gene = unique(matrix(unlist(strsplit(expr.raw[,1],'|',fixed=T)),byrow=T,ncol=2)[,2])
write(gene, ncol=1, file="genes_in_assay.txt", sep="\t")

# The first two columns contain gene symbols/ids, so use entrez ID as row names and remove the columns
rownames(expr.raw) <- expr.raw[,1]
expr.raw <- expr.raw[,-1]
colnames(expr.raw) <- unlist(lapply(strsplit(colnames(expr.raw),"\\."),"[[",7))# This parses the column names and keeps only the relevant ID string

sample_info_T <- sample_info[,"RNASeq_T"]
expr.raw_T <- expr.raw[, match(sample_info_T, colnames(expr.raw), nomatch = 0)]
sample_info_NT <- sample_info[,"RNASeq_NT"]
expr.raw_NT <- expr.raw[, match(sample_info_NT, colnames(expr.raw), nomatch = 0)]
# Gene filtering criteria:
# select = (rowMeans(expr.raw_T)>1)|(rowMeans(expr.raw_NT)>1) # alternative (more lenient) filtering criteria 
select = rowMeans(expr.raw_T)>1
expr.raw <- expr.raw[select,]

##### Normalize expression data: Voom models the mean-variance trend and adjusts count values to it, does quantile normalization, then transforms to log-2
expr <- voom(expr.raw,normalize.method = "quantile")$E

##### Write normalized expression data to file
output.df <- data.frame(entrez_id=rownames(expr), expr)
write.table(output.df, "gene_expr_voom.txt", sep="\t", quote=F, row.names=F)
