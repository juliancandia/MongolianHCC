library(circlize)
library(ComplexHeatmap)

rm(list=ls())

PROJECT_DIR = "~/ACTIVE/LCS/Mongolia_HCC_Pipeline" # replace this line with your local path

source(file.path(PROJECT_DIR,"SCRIPTS","helper_functions.oncoplot.R"))

#### Define locations of input/output
mydata_file = file.path(PROJECT_DIR,"RESULTS","cluster","FINAL","gene_expr_sig.txt")
sample_info_file = file.path(PROJECT_DIR,"DATA","PROCESSED","patient_sample_metadata_w_clustering_risk.txt")
out_dir = file.path(PROJECT_DIR,"RESULTS","cluster","FINAL")

if (! dir.exists(out_dir)) {
  dir.create(out_dir, recursive = T)
}

#### Read in data
mydata <- read.table(mydata_file,header=TRUE,stringsAsFactors = FALSE, check.names=FALSE,row.names = 1)
# gene_clusters <- data.frame(class=paste0("MO",mydata[,1]), row.names = rownames(mydata))
gene_clusters <- data.frame(class=mydata[,1], row.names = rownames(mydata), stringsAsFactors = F)
mydata <- mydata[,-1] ## First column contains cluster assignment

sample_info <- read.table(sample_info_file, header=T, sep="\t", quote="")
sample_info$Patient = paste0("PID_",sample_info$Patient) # to fix mismatch
rownames(sample_info) <- sample_info$Patient
## Match the two dfs by the Patient ID
hm_anno_info <- as.data.frame(sample_info[match(colnames(mydata),sample_info$Patient),])
## Fix sex to M/F to match with color annotation
hm_anno_info$sex <- ifelse(is.na(hm_anno_info$sex), "NA", ifelse(hm_anno_info$sex, "M", "F"))

## Define which variable to display on top/bottom; names are col names in hm_anno_info; values are "pretty" labels for printing
top_anno_names <- c(class="Class",hdv="HDV",hbv="HBV",hcv="HCV")
bot_anno_names <- c(sex="Sex",afp="AFP",cirrhosis="Cirrhosis",tumor_size="Tumor Size")
#top_anno_names <- c(class="Class")
#bot_anno_names <- c(sex="Sex",stage="Stage",tumor_size="Tumor Size",cirrhosis="Cirrhosis",hcv="HCV",hbv="HBV",hdv="HDV")
hm_anno_info <- hm_anno_info[,c(names(top_anno_names), names(bot_anno_names))]

## Now make colors
my_colors <- my_oncoplot_colors(hm_anno_info)

## Organize mydata for plotting... 
hm_data <- mydata
samples_by_cluster <- sapply(split(hm_anno_info,hm_anno_info$class), rownames)
# Match sample order to semi-supervised clustering
semi_hc <- hclust_semisupervised(t(hm_data), samples_by_cluster, dist_method = "euclidean", hclust_method = "average")
hm_data <- hm_data[,match(rownames(semi_hc$data), colnames(hm_data))]
hm_anno_info <- hm_anno_info[match(rownames(semi_hc$data), rownames(hm_anno_info)),]  ## This is prob not necessary, ComplexHeatmap should take care of matching the annotations to the data

## Make HeatmapAnnotation objects for top/bottom annotations
top_anno_data <- data.frame(hm_anno_info[,names(top_anno_names)], row.names=rownames(hm_anno_info))
colnames(top_anno_data) <- names(top_anno_names)
top_anno <- HeatmapAnnotation(df=top_anno_data, name="top_anno",col=my_colors,show_annotation_name=F,na_col="grey70", show_legend = T)
names(top_anno) <- top_anno_names[names(top_anno)]

bot_anno_data <- data.frame(hm_anno_info[,names(bot_anno_names)], row.names=rownames(hm_anno_info))
colnames(bot_anno_data) <- names(bot_anno_names)
bot_anno <- HeatmapAnnotation(df=bot_anno_data, name="bot_anno",col=my_colors,show_annotation_name=T,na_col="white", show_legend = T)
names(bot_anno) <- bot_anno_names[names(bot_anno)]

lef_anno_data <- gene_clusters
# rownames(lef_anno_data) <- names(bot_anno_names)
lef_anno <- rowAnnotation(df=lef_anno_data, name="Gene Cluster",col=my_colors,
                          show_annotation_name=F,na_col="white", show_legend = F,annotation_width=unit(1, "mm"),width=unit(1, "mm"))
names(lef_anno) <- "gene_cluster"

## Colors for the heatmap values
colorbreak = seq(-2,2,length.out = 100)
colorVals = colorRampPalette(c("blue","yellow"))(length(colorbreak))

## Plot ComplexHeatmap
hm_data <- as.matrix(hm_data)
rnaHM <- Heatmap(hm_data, name="Normalized Expression",
                cluster_rows=T,
                cluster_columns=semi_hc$hclust,
                show_row_names = F,
                show_column_names = F,
                rect_gp = gpar(col = NA, lwd = 1),
                # row_split=split_idx,
                col = colorRamp2(breaks = colorbreak, colors=colorVals),
                clustering_method_columns="centroid",
                top_annotation=top_anno,
                bottom_annotation = bot_anno,
                left_annotation = lef_anno,
                heatmap_legend_param = list(
                  title_position = "leftcenter-rot",
                  # legend_direction = "horizontal"
                  title = "Normalized Expression\n(Z-Score)"
                ),
                # heatmap_width = unit(0.3, "npc"),
                width = unit(0.8, "npc"),
                show_row_dend = F)

## Print to file
pdf(file = paste0(out_dir,"/cluster_heatmap.pdf"))
draw(rnaHM,heatmap_legend_side = "left")
## This adds the class labels to the annotation directly
decorate_annotation("Class", {
  fracs=table(top_anno_data$class)/length(top_anno_data$class)
  pos_vec = c(0,cumsum(fracs),1)
  midpoints = diff(pos_vec)/2 + pos_vec[-length(pos_vec)]
  midpoints = midpoints[1:(length(midpoints)-1)]
  for (i in 1:length(midpoints)) {
    # print(paste0(names(midpoints)[i], ": ", midpoints[i]))
    grid.text(names(midpoints)[i], midpoints[i], 0.5, default.units = "npc",gp=gpar(fontsize=9, fontface="bold",col="grey10"))
  }
})

dev.off()


