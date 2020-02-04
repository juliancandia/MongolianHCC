library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)

rm(list=ls())

PROJECT_DIR = "~/ACTIVE/LCS/Mongolia_HCC_Pipeline" # replace this line with your local path

infile = file.path(PROJECT_DIR,"RESULTS","mut_sig","weights_HDV_annotated.txt")
data = read.table(infile, header=T, sep="\t", quote="")
sig = data[,1:7]
weights = data[,-(1:7)]
infile = file.path(PROJECT_DIR,"DATA","PROCESSED","patient_sample_metadata_w_clustering_risk.txt")
metadata = read.table(infile, header=T, sep="\t", quote="")
hm_anno_info <- as.data.frame(metadata[match(colnames(weights),metadata$WES_T),])
rownames(hm_anno_info) <- hm_anno_info$WES_T

## Define which variable to display on top/bottom; names are col names in hm_anno_info; values are "pretty" labels for printing
top_anno_names <- c(class="Class",hdv="HDV",hbv="HBV",hcv="HCV")
bot_anno_names <- c("Sex","Age","Obesity","Smoking","Alcohol","Fam Hx LC","Cirrhosis",
"Albumin","Bilirubin","ALT","AFP","Tumor Size","Multinodular","TNM Staging","Risk (Cox)")
names(bot_anno_names) <- c("sex","age_bin","obesity","smoker","alcohol","FHx_LC","cirrhosis","alb","bil","alt","afp","tumor_size","multinodular","stage","risk_bin")
hm_anno_info <- hm_anno_info[,c(names(top_anno_names), names(bot_anno_names))]

## Makes the colors used for annotations
my_oncoplot_colors <- function(hm_anno_info) {
  cluster_colors = c(MO1="#4DAF4A",MO2="#377EB8",MO3="orange",MO4="#E41A1C","NA"="white")
  sex_colors <- c("1"="blue","0"="hotpink","NA"="white")
  tnm_colors <- c(brewer.pal(4, "Reds"),"grey70")
  names(tnm_colors) <- c("1","2","3","4","NA")
  
  def_bin_colors <- c("0" = "gray70","1" = "grey20","NA"="white")
  
  all_colors <- rep(list(def_bin_colors),ncol(hm_anno_info))
  names(all_colors) <- colnames(hm_anno_info)
  all_colors[[grep("sex",colnames(hm_anno_info), ignore.case = T)]] <- sex_colors
  all_colors[[grep("class",colnames(hm_anno_info), ignore.case = T)]] <- cluster_colors
  
  return(all_colors)
}
my_colors <- my_oncoplot_colors(hm_anno_info)

## Make HeatmapAnnotation objects for top/bottom annotations
top_anno_data <- data.frame(hm_anno_info[,names(top_anno_names)], row.names=rownames(hm_anno_info))
colnames(top_anno_data) <- names(top_anno_names)
top_anno <- HeatmapAnnotation(df=top_anno_data, name="top_anno",col=my_colors,show_annotation_name=T,
                              na_col="grey70", show_legend = T)
names(top_anno) <- top_anno_names[names(top_anno)]

bot_anno_data <- data.frame(hm_anno_info[,names(bot_anno_names)], row.names=rownames(hm_anno_info))
colnames(bot_anno_data) <- names(bot_anno_names)
bot_anno <- HeatmapAnnotation(df=bot_anno_data, name="bot_anno",col=my_colors,show_annotation_name=T,
                              na_col="white", show_legend = F,
                              annotation_name_gp=gpar(cex=0.8))
names(bot_anno) <- bot_anno_names[names(bot_anno)]

signif = rep(NA,nrow(sig))
signif[sig[,"p.value"]<0.05] = "0.01 < p < 0.05"
signif[sig[,"p.value"]<0.01] = "0.005 < p < 0.01"
signif[sig[,"p.value"]<0.005] = "p < 0.005"
hdv_prev = as.character(sig[,"Prevalence"])
hdv_prev[hdv_prev=="HDV_pos"] = "HDV+"
hdv_prev[hdv_prev=="HDV_neg"] = "HDV-"
left_anno_data = data.frame("Prevalence" = hdv_prev,"Significance" = signif)
hdv_prev_col = rev(brewer.pal(4, "Reds")) [c(1,4)] #c("red","blue")
names(hdv_prev_col) = c("HDV+","HDV-")
signif_col = rev(brewer.pal(3,"Blues"))
names(signif_col) = rev(c("0.01 < p < 0.05","0.005 < p < 0.01","p < 0.005"))
left_anno_col = list("Prevalence" = hdv_prev_col, "Significance" = signif_col)
left_anno <- HeatmapAnnotation(df=left_anno_data, name="left_anno",col=left_anno_col,show_annotation_name=T,na_col="grey70", show_legend = T,which = "row")

right_anno_data = data.frame("Reference" = sig[,"Reference"])
right_anno_col = list("Reference" = c("COSMIC_v2"="#00DAE0","COSMIC_v3"="#59FF4A","EnvirAg"="#FF9289"))
right_anno <- HeatmapAnnotation(df=right_anno_data, name="right_anno",col=right_anno_col,show_annotation_name=F,na_col="grey70", show_legend = T,which = "row")

## Colors for the heatmap values
colorbreak = seq(0,quantile(unlist(weights),probs=0.95),length.out=50)
colorVals = colorRampPalette(c("blue","yellow"))(length(colorbreak))

hm_data <- as.matrix(weights)
rownames(hm_data) = sig[,3]
myHM <- Heatmap(hm_data, name="Subject/Signature Weights",
                cluster_rows=F,
                cluster_columns=T,
                show_row_names = T,
                show_column_names = F,
                rect_gp = gpar(col = NA, lwd = 1),
                row_names_gp = gpar(cex=0.8),
                col = colorRamp2(breaks = colorbreak, colors=colorVals),
                clustering_method_columns="complete",
                top_annotation=top_anno,
                bottom_annotation = bot_anno,
                left_annotation = left_anno,
                right_annotation = right_anno,
                show_row_dend = F)

## Print to file
outfile = file.path(PROJECT_DIR,"RESULTS","mut_sig","weights_HDV.pdf")
pdf(file = outfile)
draw(myHM,heatmap_legend_side = "left")
dev.off()
