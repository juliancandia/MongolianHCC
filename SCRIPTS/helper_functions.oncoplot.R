
make_TCGA_comparison_heatmap <- function(onco_genes, mong_maf, split_at=NULL) {
  #### Load TCGA-LIHC maf data
  tcga_load(study = "LIHC")  ## Loads maf object named 'tcga_lihc_mc3'
  tcga_clindata <- getClinicalData(tcga_lihc_mc3)
  
  # Asian only
  asian_samples <- levels(droplevels(tcga_clindata$Tumor_Sample_Barcode[tcga_clindata$race_list=="ASIAN"]))
  tcga.asian <- subsetMaf(tcga_lihc_mc3, tsb = asian_samples,mafObj = TRUE) 
  
  # Caucasian only
  white_samples <- levels(droplevels(tcga_clindata$Tumor_Sample_Barcode[tcga_clindata$race_list=="WHITE"]))
  tcga.white <- subsetMaf(tcga_lihc_mc3, tsb = white_samples,mafObj = TRUE) 
  
  # Summarize altered frac for each dataset
  my_gene_ind.mong <- match(onco_genes,mong_maf@gene.summary$Hugo_Symbol)
  my_gene_ind.tcga <- match(onco_genes,tcga_lihc_mc3@gene.summary$Hugo_Symbol)
  my_gene_ind.asian <- match(onco_genes,tcga.asian@gene.summary$Hugo_Symbol)
  my_gene_ind.white <- match(onco_genes,tcga.white@gene.summary$Hugo_Symbol)
  altered_fracs <- data.frame(
    MONG=(mong_maf@gene.summary$AlteredSamples/as.numeric(mong_maf@summary$summary[3]))[my_gene_ind.mong],
    TCGA_All=(tcga_lihc_mc3@gene.summary$AlteredSamples/as.numeric(tcga_lihc_mc3@summary$summary[3]))[my_gene_ind.tcga],
    TCGA_Asian=(tcga.asian@gene.summary$AlteredSamples/as.numeric(tcga.asian@summary$summary[3]))[my_gene_ind.asian],
    TCGA_White=(tcga.white@gene.summary$AlteredSamples/as.numeric(tcga.white@summary$summary[3]))[my_gene_ind.white],
    row.names = onco_genes)
  
  
  hm_data <- as.matrix(altered_fracs*100)
  colnames(hm_data) <- c("Mongolian","TCGA (All)", "TCGA (Asian)", "TCGA (Cauc)")
  
  ## Define breaks for colors
  frac_breaks <- c(seq(round(min(range(hm_data, na.rm=T)),3),
                       round(mean(unlist(hm_data), na.rm=T), 3), by=0.01),
                   seq(round(mean(unlist(hm_data), na.rm=T), 3),
                       round(max(range(hm_data, na.rm=T)),3), by=0.01))
  frac_colors <- colorRampPalette(c("#fffadb","red"))(length(frac_breaks))
  
  hm_data[is.na(hm_data)] <- 0
  writeValue <- function(j, i, x, y, width, height, fill) {
    mytext <- ifelse(hm_data[i,j] > 0 & hm_data[i,j] < 1,paste0(round(hm_data[i, j],1),"%"),
                     paste0(round(hm_data[i, j],0),"%"))
    
    grid.text(mytext, x, y, gp = gpar(fontsize = 8, col="black"))
  }
  tcga_altered_frac_HM <- Heatmap(hm_data, name="% Mutated Samples",
                                  cluster_rows=F,
                                  # row_title_rot = 30,
                                  row_split=split_at,
                                  col = colorRamp2(frac_breaks,frac_colors),
                                  clustering_method_columns="centroid",
                                  cell_fun = writeValue,
                                  column_dend_reorder = c(1,0,0,0),
                                  # width = unit(0.4, "npc"))
                                  width = 0.5)
  return(list(frac_hm=tcga_altered_frac_HM, tcga_white=tcga.white, tcga_asian=tcga.asian, tcga_all=tcga_lihc_mc3))
  
}



make_mongolia_oncoplot_grouped <- function(maf_file, sample_info_file, driver_res_file=NULL, out_dir="mongolia_oncoplot_out",ignoreGenes = "default",
                                           filter_only=F,myfilters=NULL, adjust_q_val=F,
                                           driver_sig_col = "q", driver_sig_thresh = 0.1, driver_sig_freq = 0, cohort_freq_thresh = 0.095,
                                           include_all=T,annotate_empty = "", onco_width=8, onco_height=NULL, annotation_height_frac = 0.3, annotation_font_size = 9,
                                           ...) {
          
  if (! file.exists(maf_file)) {
    stop(paste0("Can't find MAF file: ", maf_file))
  }
  
  if (! file.exists(sample_info_file)) {
    stop(paste0("Can't find sample info file: ", sample_info_file))
  }
  
  if (! dir.exists(out_dir)) {
    dir.create(out_dir, recursive = T)
  }
  
  if (length(ignoreGenes)==0) {
    ignoreGenes <- c()
  } else if (ignoreGenes[1]=="default") {
    ignoreGenes <- c("TTN","MUC16","OBSCN","AHNAK2","SYNE1","FLG","MUC5B","DNAH17","PLEC","DST","SYNE2","NEB","HSPG2","LAMA5","AHNAK","HMCN1","USH2A","DNAH11","MACF1","MUC17","DNAH5","GPR98","FAT1","PKD1","MDN1","RNF213","RYR1","DNAH2","DNAH3","DNAH8","DNAH1","DNAH9","ABCA13","SRRM2","CUBN","SPTBN5","PKHD1","LRP2","FBN3","CDH23","DNAH10","FAT4","RYR3","PKHD1L1","FAT2","CSMD1","PCNT","COL6A3","FRAS1","FCGBP","RYR2","HYDIN","XIRP2","LAMA1")
  }
  
  # browser()
  maf_df <- filter_maf(maf_file, flag_genes = ignoreGenes, ...)
  
  sample_info <- read.table(sample_info_file, header=T, sep="\t", quote="")
  sample_info$Patient <- paste0("Patient_",sample_info$Patient)
  sample_info.exome <- sample_info
  colnames(sample_info.exome)[colnames(sample_info.exome)=="WES_T"] <- "Tumor_Sample_Barcode"
  
  mafObj <- read.maf(maf_df,clinicalData = sample_info.exome)
  
  maf.filter <- mafObj
  if (length(myfilters)>0) {
    print("Applying filters to MAF object...")
    for (filter_str in myfilters) {
      maf.filter<-subsetMaf(maf=maf.filter,query = filter_str,mafObj = TRUE,restrictTo = "mutations")
    }
  } else {
    print("Skipping filtering...")
  }
  
  # save
  filtered_maf_path<-paste0(out_dir,"/filtered.maf")
  print(paste0("Saving filtered MAF to ", filtered_maf_path))
  write.table(rbind(maf.filter@data, maf.filter@maf.silent), sep="\t", quote=F, file = filtered_maf_path, row.names = F)
  
  if (!filter_only) {
    driver_res <- read.table(driver_res_file, sep="\t",header=T, quote="", stringsAsFactors = F)
    colnames(driver_res)[colnames(driver_res)=="gene"] <- "Hugo_Symbol"
    driver_res$FLAG_gene <- driver_res$Hugo_Symbol %in% ignoreGenes
    
    if (adjust_q_val) {
      driver_res <- driver_res[driver_res$Hugo_Symbol %in% maf.filter@gene.summary$Hugo_Symbol,]
      driver_res[,driver_sig_col] <- p.adjust(driver_res$p,method="fdr")
    }
    
    frac_mut <- data.frame(MONG_frac_mut=(maf.filter@gene.summary$MutatedSamples/as.numeric(maf.filter@summary$summary[3])), 
                           Hugo_Symbol=maf.filter@gene.summary$Hugo_Symbol, stringsAsFactors = F)
    
    driver_res_plus <- merge.data.frame(driver_res, frac_mut)
    driver_res_plus <- driver_res_plus[order(driver_res_plus[,driver_sig_col], decreasing = F),]
    driver_genes <- driver_res_plus$Hugo_Symbol[driver_res_plus[,driver_sig_col] < driver_sig_thresh & driver_res_plus$MONG_frac_mut > driver_sig_freq]
    freq_genes <- setdiff(driver_res_plus$Hugo_Symbol[driver_res_plus$MONG_frac_mut > cohort_freq_thresh], driver_genes)
    
    gene_list <- list(driver_genes, freq_genes)
    reasons <- c(paste0("Driver Gene, ", driver_sig_col, " < ",driver_sig_thresh),
                 paste0("Cohort Freq > ",cohort_freq_thresh))
    
    genes_for_oncoplot <- data.frame(Hugo_Symbol=c(), reason=c())
    for (i in 1:length(gene_list)) {
      if (is.na(gene_list[[i]][1])) {
        next
      }
      genes_for_oncoplot <- rbind(genes_for_oncoplot,
                                  data.frame(Hugo_Symbol=gene_list[[i]],
                                             reason=reasons[i]))
    }
    
    genes_for_oncoplot <- cbind(genes_for_oncoplot,
                                q=driver_res_plus[match(genes_for_oncoplot$Hugo_Symbol, driver_res_plus$Hugo_Symbol),driver_sig_col],
                                frac=driver_res_plus$MONG_frac_mut[match(genes_for_oncoplot$Hugo_Symbol, driver_res_plus$Hugo_Symbol)])
    # browser()
    genes_for_oncoplot <- genes_for_oncoplot[order(genes_for_oncoplot$reason, -genes_for_oncoplot$frac),]
    split_idx=genes_for_oncoplot$reason
    split_colors <- rainbow(length(levels(split_idx)))
    names(split_colors) <- as.character(genes_for_oncoplot$reason[!duplicated(genes_for_oncoplot$reason)])
    split_colors <- list(Reason=split_colors)
    
    # source("scripts/helper_functions.R")
    oncomat <- createOncoMatrix(maf.filter, g=genes_for_oncoplot$Hugo_Symbol, add_missing = F)$oncoMatrix
    # write.table(oncomat,file = paste0(figures_folder,"/oncomat.txt"),sep = "\t")
    
    # maf.oncoplot <- subsetMaf(maf.filter, genes = genes_for_oncoplot$Hugo_Symbol)
    # save(maf.oncoplot, tcga_lihc_mc3, tcga.asian, tcga.white, file = paste0(data_folder,"/maf/oncoplot.",suffix,".maf.Rdata"))
    # save(maf.oncoplot, file = paste0(out_dir,"/oncoplot.mafs.Rdata"))
    write.table(genes_for_oncoplot, file = paste0(out_dir,"/genes_for_oncoplot.txt"), sep="\t", quote=F, row.names=F)
    
    
    if (include_all) {
      ### createOncoMatrix drops empty samples, so this adds them back in
      all_wes_samples <- as.character(sample_info.exome$Tumor_Sample_Barcode[!is.na(sample_info.exome$Tumor_Sample_Barcode)])
      extra_samples <- setdiff(all_wes_samples, colnames(oncomat) )
      empty_data <- matrix(data = "", nrow=nrow(oncomat), ncol=length(extra_samples), dimnames=list(rownames(oncomat), extra_samples))
      oncomat <- cbind(oncomat, empty_data)
    }
    
    ## Order rows by selected gene order
    oncomat <- oncomat[match(genes_for_oncoplot$Hugo_Symbol,rownames(oncomat)), ]
    onco_genes <- rownames(oncomat)
    
    
    tcga_comparison_results <- make_TCGA_comparison_heatmap(onco_genes, maf.filter, split_at = split_idx)
    tcga_comparison_hm <- tcga_comparison_results[[1]]
    
    
    # browser()
    #### MONGOLIA oncoplot
    oncomat.plot <- oncomat
    colnames(oncomat.plot) <- sample_info.exome$Patient[match(colnames(oncomat.plot), sample_info.exome$Tumor_Sample_Barcode)]
    
    if (is.null(onco_height)) {
      onco_height=max(round(0.3*nrow(oncomat.plot),0),5)
    }
    
    hm_anno_info <- as.data.frame(sample_info.exome[match(colnames(oncomat.plot),sample_info.exome$Patient),])
    rownames(hm_anno_info) <- hm_anno_info$Tumor_Sample_Barcode
    hm_anno_info <- hm_anno_info[,c(10:ncol(hm_anno_info))]
    hm_anno_info$sex <- ifelse(is.na(hm_anno_info$sex), "NA", ifelse(hm_anno_info$sex==1, "M", "F"))
    
    top_anno_names <- c(class="Class")
    bot_anno_names <- c("Sex","Age","HCV","HBV","HDV","Stage","Tumor Size","Cirrhosis","Obesity","Alcohol","Smoker","Heredit Liver Dis","Risk Group") 
    names(bot_anno_names) <- c("sex","age_bin","hcv","hbv","hdv","stage","tumor_size","cirrhosis","obesity","alcohol","smoker","hereditary.liver.disease","risk_bin") # column names
    # bot_anno_names <- c(sex="Sex",stage="Stage",tumor_size="Tumor Size",cirrhosis="Cirrhosis",
    #                     hcv="HCV",hbv="HBV",hdv="HDV",risk_bin="Risk Group")
    # 
    hm_anno_info <- hm_anno_info[,c(names(top_anno_names), names(bot_anno_names))]
    
    annocolors <- my_oncoplot_colors(hm_anno_info)
    
    
    my_types <- unique(unlist(apply(oncomat.plot,2,unique)))
    my_types <- my_types[!my_types %in% c(NA,"",0)]
    
    col <- c(Nonsense_Mutation="#ad7aff",Missense_Mutation="#377EB8",Frame_Shift_Del="#4DAF4A",
             In_Frame_Ins="#ff008c",Splice_Site="#FF7F00",Multi_Hit="#FFFF33",Frame_Shift_Ins="#A65628",
             In_Frame_Del="#f781bf",Translation_Start_Site="#400085",Nonstop_Mutation="#b68dfc",
             no_variants="#d6d6d6")
    
    top_anno_data <- data.frame(class=hm_anno_info$class, row.names=rownames(hm_anno_info))
    bot_anno_data <- data.frame(hm_anno_info[,! colnames(hm_anno_info) %in% "class"])
    
    order_classes <- factor(ifelse(is.na(as.character(top_anno_data$class)),"Unknown",as.character(top_anno_data$class)))
    plot_order <- orderByGroup(oncomat.plot, order_classes)
    # oncomat.plot <- oncomat.plot[,plot_order]
    
    unmutated_annodata <- ifelse(colSums(nchar(oncomat.plot))==0,annotate_empty,"")
    # oncomat.plot[,colSums(nchar(oncomat.plot))==0] <- "no_variants"
    
    # browser()
    annocolors$empty <- c("TRUE"="black","FALSE"="white")
    top_anno <- HeatmapAnnotation(empty=anno_text(unmutated_annodata, gp=gpar(fontsize=6, fontface="bold",col="grey10"),location=unit(0.5, 'npc'),which="column"),
                                  df=top_anno_data, name="top_anno",col=annocolors,show_annotation_name=F,na_col="grey70", show_legend = F, height=unit(8,'mm'))
    # draw(top_anno)
    bot_anno <- HeatmapAnnotation(df=bot_anno_data, name="bot_anno",col=annocolors,show_annotation_name=T,na_col="white", show_legend = F, 
                                  simple_anno_size_adjust=T, height=unit(onco_height*annotation_height_frac,'inches'),
                                  annotation_name_gp = gpar(fontsize = annotation_font_size))
    names(bot_anno) <- bot_anno_names[names(bot_anno)]
    
    # browser()
    names(col) <- gsub("_"," ",names(col))
    oncomat.plot <- gsub("_"," ",oncomat.plot)
    
    onco_base_default <- oncoPrint(oncomat.plot, alter_fun = alter_fun, col=col, row_order=1:nrow(oncomat.plot),
                                   name="oncoplot",
                                   show_pct = F,
                                   top_annotation=top_anno,
                                   bottom_annotation = bot_anno,
                                   row_split=split_idx,
                                   left_annotation = rowAnnotation(Reason = split_idx, col=split_colors, annotation_width = unit(0.5, "mm")),
                                   row_title = NULL,
                                   column_title=NULL,
                                   column_order = plot_order,
                                   column_gap = unit(0.001,"npc"),
                                   column_split=top_anno_data$class,
                                   width = unit(0.75, "npc"))
    
    # browser()
    save_name <- paste0(out_dir,"/oncoplot.pdf")
    
    
    pdf(file = save_name,height=onco_height,width=onco_width)
    draw(tcga_comparison_hm + onco_base_default, main_heatmap = 2)
    
    myclasses <- levels(addNA(top_anno_data$class))
    myclasses[is.na(myclasses)] <- ""
    for (i in 1:length(myclasses)) {
      decorate_annotation("class", slice=i, {
        grid.text(myclasses[i], 0.5,
                  0.5, default.units = "npc",gp=gpar(fontsize=9, fontface="bold",col="grey10"))
      })
      for (j in 1:(length(levels(split_idx))-1)) {
        decorate_heatmap_body("oncoplot", row_slice = j, column_slice=i, {
          grid.lines(unit(c(1, 0), "npc"), unit(c(0,0), "npc"), gp = gpar(lty = 1, lwd = 4))
        })
      }
      
    }
    
    dev.off()
    
    
    
    
    ####  TABLE containing data for each variant
    output_data <- maf.filter@data[maf.filter@data$Hugo_Symbol %in% onco_genes,]
    output_data$tumor_genotype <- apply(output_data[,c("Tumor_Seq_Allele1","Tumor_Seq_Allele2")], 1, paste, collapse="/")
    output_data$normal_genotype <- apply(output_data[,c("Match_Norm_Seq_Allele1","Match_Norm_Seq_Allele2")], 1, paste, collapse="/")
    
    pheno_info <- sample_info.exome[match(output_data$Tumor_Sample_Barcode, sample_info.exome$Tumor_Sample_Barcode),]
    pheno_info <- cbind(pheno_info[,"Tumor_Sample_Barcode"],pheno_info[,-c("Tumor_Sample_Barcode")])
    pheno_columns <- colnames(pheno_info)
    names(pheno_columns) <- make.names(pheno_columns, unique = T)
    
    output_data <- cbind(output_data,pheno_info)
    cols_for_table <- c("Hugo Symbol" = "Hugo_Symbol",
                        "Variant Classification"="Variant_Classification",
                        "Variant Type"="Variant_Type",
                        "Consequence"="Consequence",
                        "Chromosome"="Chromosome","Start Position" ="Start_Position","End Position"="End_Position","Strand"="Strand",
                        "Reference Allele"="Reference_Allele",
                        "Tumor Genotype"="tumor_genotype",
                        "Normal Genotype"="normal_genotype",
                        "Transcript Change"="HGVSc",
                        "Protein Change"="HGVSp_Short",
                        "Normal Depth"="n_depth",
                        "Normal Ref Depth"="n_ref_count",
                        "Normal Alt Depth"="n_alt_count",
                        "Tumor Depth"="t_depth",
                        "Tumor Ref Depth"="t_ref_count",
                        "Tumor Alt Depth"="t_alt_count",
                        "Existing Annotation"="Existing_variation",
                        "gnomAD Frequency"="gnomAD_AF",
                        "ExAC Frequency"="ExAC_AF",
                        "1000Genomes Frequency"="AF",
                        "Current Cohort Frequency"="tumor_freq",
                        pheno_columns
    )
    
    variant_info <- as.data.frame(output_data)[,cols_for_table]
    colnames(variant_info) <- names(cols_for_table)
    write.xlsx(variant_info,file = paste0(out_dir,"/Table_2.xlsx"))
    
    tcga.white <- tcga_comparison_results[[2]]
    tcga.asian <- tcga_comparison_results[[3]]
    tcga_lihc_mc3 <- tcga_comparison_results[[4]]
    all_dfs <- list(
      driver_res_plus,
      data.frame(MONG_frac_mut=(maf.filter@gene.summary$MutatedSamples/as.numeric(maf.filter@summary$summary[3])), 
                 Hugo_Symbol=maf.filter@gene.summary$Hugo_Symbol, stringsAsFactors = F),
      data.frame(TCGA_All_frac_mut=(tcga_lihc_mc3@gene.summary$MutatedSamples/as.numeric(tcga_lihc_mc3@summary$summary[3])), 
                 Hugo_Symbol=tcga_lihc_mc3@gene.summary$Hugo_Symbol, stringsAsFactors = F),
      data.frame(TCGA_Asian_frac_mut=(tcga.asian@gene.summary$MutatedSamples/as.numeric(tcga.asian@summary$summary[3])), 
                 Hugo_Symbol=tcga.asian@gene.summary$Hugo_Symbol, stringsAsFactors = F),
      data.frame(TCGA_White_frac_mut=(tcga.white@gene.summary$MutatedSamples/as.numeric(tcga.white@summary$summary[3])), 
                 Hugo_Symbol=tcga.white@gene.summary$Hugo_Symbol, stringsAsFactors = F),
      data.frame(maf.filter@gene.summary)
    )
    
    
    all_dfs_merged <- all_dfs %>%
      Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by="Hugo_Symbol"), .)
    
    write.xlsx(all_dfs_merged, file=paste0(out_dir,"/gene_mutsig_info.xlsx"))#,asTable = T)
  }
  return(maf.filter)
  
}




### Makes commands for ALVIEW to run in batch mode
### ALVIEW will creates a snapshot of alignments at each variant in each sample
make_alview_batch_file <- function(maf_data, gene_list=NULL, top_n=NULL,
                                   alview_bin="./alview_cmdline_linux",
                                   bam_file_loc="/data/CCBR/projects/ccbr921/exomeseq_hg38/",
                                   genome_build="hg38",
                                   pairs_file = "pairs", out_folder="alview_cmd_files") {
 
  # print("HELLOW!")
  if (class(maf_data) == "MAF") {
    mutect2maf = maf_data
  } else if (class(maf_data) == "character") {
    if (file.exists(maf_data)) {
      print("Reading MAF data from file...")
      mutect2maf <- read.maf(maf_file)
      maf<-subsetMaf(maf=mutect2maf,query = "gnomAD_AF<0.001",mafObj = TRUE,restrictTo = "mutations")
      maf<-subsetMaf(maf=maf,query = "ExAC_AF<0.001",mafObj = TRUE,restrictTo = "mutations")
      maf<-subsetMaf(maf=maf,query = "AF<0.01",mafObj = TRUE,restrictTo = "mutations")
      maf<-subsetMaf(maf=maf,query = "tumor_freq>0.05",mafObj = TRUE,restrictTo = "mutations")
    } else {
      stop(paste0("Can't find file: ", maf_data))
    }
  } else {
    stop(paste0("Don't know what to do with oncoplot_data of class : ", class(maf_data)))
    
  }
  
  
  if (!dir.exists(out_folder)) {
    dir.create(out_folder,recursive = T)
  }
  
  # browser()
  if (is.null(top_n)) {
    top_n = nrow(mutect2maf@gene.summary)
  }
  top_genes <- mutect2maf@gene.summary$Hugo_Symbol[1:top_n]
  
  
  
  mafdata <- mutect2maf@data
  curr_maf <- mafdata[mafdata$Hugo_Symbol %in% top_genes,]
  curr_maf <- curr_maf[order(curr_maf$Hugo_Symbol),]
  
  # If top_n not specified, use all genes in MAF
  cmd_data <- curr_maf[,c("Chromosome","Start_Position")]
  
  my_pairs <- lapply(readLines(pairs_file), function(x){unlist(strsplit(x, "\t"))})
  
  file_name <- file.path(out_folder,"alview.cmds.txt")
  if (file.exists(file_name)) {
    file.remove(file_name)
  }
  print(paste0("outputing commands to file: ",file_name))
  # browser()
  alview_files <- lapply(my_pairs, function(curr_pair, my_maf) {
    
    # curr_pair=my_pairs[[3]]
    # my_maf=curr_maf
    sample_name <- gsub("A_S.*$","",curr_pair[2])
    
    my_maf <- my_maf[my_maf$Tumor_Sample_Barcode %in% curr_pair[2],]
    
    # ./prog inbam outpng position build [imageheight] [imagewidth]
    n=nrow(my_maf)
    if (n > 0) {
      my_data_normal <- cbind(my_maf[,c("Chromosome","Start_Position","Hugo_Symbol","Variant_Classification")], 
                              sample_name=rep(curr_pair[1],each=nrow(my_maf)))
      my_data_normal.cmd <- data.frame(
                              alview_bin=alview_bin,
                              file_loc=paste0(bam_file_loc,my_data_normal$sample_name,".recal.bam"),
                              out_png=paste0("alview_out/",my_data_normal$Hugo_Symbol,"_",my_data_normal$Start_Position,"_",my_data_normal$sample_name,"_Normal.png"),
                              position=paste0(my_data_normal$Chromosome,":",my_data_normal$Start_Position-20,"-",my_data_normal$Start_Position+20),
                              build=genome_build,
                              image_height=900,
                              image_width=1600
                              )
    
      
      my_data_tumor <- cbind(my_maf[,c("Chromosome","Start_Position","Hugo_Symbol","Variant_Classification")], 
                             sample_name=rep(curr_pair[2],each=nrow(my_maf)))
      my_data_tumor.cmd <-data.frame(
                             alview_bin=alview_bin,
                             file_loc=paste0(bam_file_loc,my_data_tumor$sample_name,".recal.bam"),
                             # bam_label=paste0(my_data_tumor$Hugo_Symbol,"_",my_data_tumor$sample_name,"_Tumor"))
                             out_png=paste0("alview_out/",my_data_tumor$Hugo_Symbol,"_",my_data_tumor$Start_Position,"_",my_data_tumor$sample_name,"_Tumor.png"),
                             position=paste0(my_data_tumor$Chromosome,":",my_data_tumor$Start_Position-20,"-",my_data_tumor$Start_Position+20),
                             build=genome_build,
                             image_height=4000,
                             image_width=6000
                             )
      
      my_data <- rbind(my_data_normal.cmd, my_data_tumor.cmd)
      s <- rep(1:n, each = 2) + (0:1) * n
      my_data <- my_data[s,]
      
      print(my_data)
      
      # file_name <- paste0(out_folder,"/",sample_name,".cmds.txt")
      # cmd_data <- my_data[,c("Chromosome","Start_Position","file_loc","bam_label")]
      write.table(my_data,file_name,sep=" ",quote=FALSE, row.names = F, col.names = F, append = T)
      return(file_name)
    } else {
      print(paste0("No Variants found for ",sample_name))
      return(NA)
    }
  },curr_maf)
  
}




make_single_ribbon_plot <- function(maf, onco_genes=NULL, save_name=NULL, ribbon_color=NULL, 
                                    pval_low=0.05, pval_high=0.1,
                                    plot_frac_mut_axis=TRUE,
                                    shrink_factor=1.3, # Higher = more shrinkage; use to control whitespace (or lack thereof) around figure 
                                    scale_ribbon_to_fracmut=TRUE) {
  # pval_low <- 0.05
  # browser()
  if (!is.null(save_name)) {
    if (! dir.exists(dirname(save_name))) {
      dir.create(dirname(save_name))
    }
    plot_file <- gsub(".pdf",".interactions.pdf",save_name)
    pdf(file = plot_file,height=5,width=5)
  }
  
  # if (is.null(onco_genes)) {
  #   onco_genes = maf@gene.summary$Hugo_Symbol
  # }
  
  som_int <-  somaticInteractions(maf = maf, genes=onco_genes, pvalue = c(pval_low, pval_high), kMax=5)
  
  if (!is.null(save_name)) {
    dev.off()
  }
  # browser()
  cooccur_data <- som_int$pairs
  cooccur_data$pair_string <- apply(cooccur_data[,1:2], 1, function(x) {paste0(sort(x), collapse="_")})
  # cooccur_data$popfrac <- cooccur_data[,"11"]/rowSums(cooccur_data[,c("00","11","01","10")], na.rm = T)  ## Frac of observed mutations
  # cooccur_data$popfrac <- cooccur_data[,"11"]/as.numeric(my_maf@summary$summary[3])   ## Frac of all samples in MAF
  cooccur_data$popfrac <- NA
  cooccur_idx <- cooccur_data$Event=="Co_Occurence"
  mut_excl_idx = cooccur_data$Event=="Mutually_Exclusive"
  
  if (scale_ribbon_to_fracmut) {
    cooccur_data$popfrac1[cooccur_idx] <- unlist(cooccur_data[cooccur_idx,"11"]/as.numeric(my_maf@summary$summary[3]))
    cooccur_data$popfrac2[cooccur_idx] <- cooccur_data$popfrac1[cooccur_idx]
    cooccur_data$popfrac1[mut_excl_idx] <- unlist(cooccur_data[mut_excl_idx,"10"]/as.numeric(my_maf@summary$summary[3]))
    cooccur_data$popfrac2[mut_excl_idx] <- unlist(cooccur_data[mut_excl_idx,"01"]/as.numeric(my_maf@summary$summary[3]))
  } else {
    cooccur_data$popfrac1 <- 1
    cooccur_data$popfrac2 <- 1
  }
  # chord_data <- cooccur_data[cooccur_data$Event=="Co_Occurance",c("gene1","gene2","popfrac","pValue")]
  chord_data <- cooccur_data[,c("gene1","gene2","popfrac1","popfrac2","pValue","Event")]
  chord_data[which(is.na(chord_data[,3])),3] <- 0
  
  require(RColorBrewer)
  # browser()
  interacting_genes <- unique(unlist(chord_data[,1:2]))
  # genecols <- colorRampPalette(brewer.pal(8,"Accent"))(length(interacting_genes))
  # genecols <- colorRampPalette(brewer.pal(8,"Dark2"))(length(interacting_genes))
  genecols <- colorRampPalette(brewer.pal(8,"Set1"))(length(interacting_genes))
  # genecols <- rainbow((length(interacting_genes)))
  names(genecols) <- interacting_genes
  
  sig_line_type=2
  link_border <- ifelse(chord_data$pValue < pval_low, sig_line_type,0)
  # color_legend <- Legend(labels=gsub("_"," ",names(cluster_sig_colors)), 
  #                        legend_gp = gpar(fill = cluster_sig_colors),background = "white",size=unit(0.08,"npc"),pch=22,
  #                        type="points",direction="horizontal",nrow=1)
  line_legend <- Legend(labels=paste0("p < ",pval_low), 
                        labels_gp = gpar(fontsize = 16, fontface = "bold"),
                        legend_gp = gpar(lty=sig_line_type,lwd = 1, fontsize=16),background = "white",
                        type="lines",direction="horizontal",nrow=1)
  # mylegend <- packLegend(color_legend, line_legend)
  if (!is.null(save_name)) {
    if (! dir.exists(dirname(save_name))) {
      dir.create(dirname(save_name))
    }
    pdf(file = save_name,height=7,width=7)
  }
  
  if (is.null(ribbon_color)) {
    # ribbon_color = genecols[match(names(genecols),c(chord_data$gene1,chord_data$gene2))]
    ribbon_color = genecols[chord_data$gene1]
    # ribbon_color[chord_data$Event == "Mutually_Exclusive"] = "grey90"
    ribbon_color[chord_data$Event == "Mutually_Exclusive"] = "white"
  }
  link_border[chord_data$Event == "Mutually_Exclusive"]=1
  # shrink_factor=1.3 # Higher = more shrinkage; use to control whitespace (or lack thereof) around figure 
  circos.clear()
  circos.par(canvas.xlim=c(-shrink_factor,shrink_factor),canvas.ylim=c(-shrink_factor,shrink_factor))
  chordDiagram(chord_data[,1:4],grid.col = genecols,
               annotationTrack = c("grid",ifelse(plot_frac_mut_axis, "axis", "")),
               col=ribbon_color,
               # transparency = link_alpha, 
               link.lty = link_border, 
               link.border = "black",
               link.sort = TRUE)
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                facing = "clockwise", niceFacing = TRUE, adj = c(-0.5, 0.5))
  }, bg.border = NA) # here set bg.border to NA is important
  # draw(color_legend, x = unit(0.5, "npc"), y = unit(0.02, "npc"), just = c("center"))
  draw(line_legend, x = unit(0.5, "npc"), y = unit(0.97, "npc"), just = c("center"))
  
  if (!is.null(save_name)) {
    dev.off()
  }
  
  
}


filter_maf <- function(maf_file, flag_genes="default",save_name=NULL,no_filter=F,
                       norm_alt_max=2,t_alt_min=4,t_depth_min=20, n_callers=2,
                       gnomAD_AF_max=0.001, AF_max=0.001, ExAC_AF_max=0.001, 
                       tumor_freq_min=0.05, norm_freq_max=0.02,
                       variant_caller=NULL) {
  # browser()
  if (length(flag_genes)==0) {
    flag_genes <- c()
  } else if (flag_genes[1]=="default") {
    flag_genes <- c("TTN","MUC16","OBSCN","AHNAK2","SYNE1","FLG","MUC5B","DNAH17","PLEC","DST","SYNE2","NEB","HSPG2","LAMA5","AHNAK","HMCN1","USH2A","DNAH11","MACF1","MUC17","DNAH5","GPR98","FAT1","PKD1","MDN1","RNF213","RYR1","DNAH2","DNAH3","DNAH8","DNAH1","DNAH9","ABCA13","SRRM2","CUBN","SPTBN5","PKHD1","LRP2","FBN3","CDH23","DNAH10","FAT4","RYR3","PKHD1L1","FAT2","CSMD1","PCNT","COL6A3","FRAS1","FCGBP","RYR2","HYDIN","XIRP2","LAMA1")
  }
  maf_df.raw <- read.table(maf_file, sep="\t", header=T, fill = T, quote="\"", stringsAsFactors = F)
  maf_df.raw <- maf_df.raw[maf_df.raw$Hugo_Symbol != "Hugo_Symbol",]
  
  if (!"tumor_freq" %in% colnames(maf_df.raw)) {
    maf_df.raw$tumor_freq <- as.numeric(maf_df.raw$t_alt_count)/as.numeric(maf_df.raw$t_depth)
  }  
  if (!"norm_freq" %in% colnames(maf_df.raw)) {
    maf_df.raw$norm_freq <- as.numeric(maf_df.raw$n_alt_count)/as.numeric(maf_df.raw$n_depth)
  }
  
  if (no_filter) {
    filter_tumor_depth=rep(TRUE,nrow(maf_df.raw))
    filter_norm_alt=rep(TRUE,nrow(maf_df.raw))
    filter_tumor_alt=rep(TRUE,nrow(maf_df.raw))
    filter_genes=rep(TRUE,nrow(maf_df.raw))
    filter_pop_freq=rep(TRUE,nrow(maf_df.raw))
  } else {
    # filter_genes=!maf_df.raw$Hugo_Symbol %in% flag_genes
    # filter_norm_alt=maf_df.raw$n_alt_count=="-" | as.numeric(maf_df.raw$n_alt_count) < norm_alt_max
    # filter_tumor_alt=as.numeric(maf_df.raw$t_alt_count) > t_alt_min & as.numeric(maf_df.raw$t_depth) > t_depth_min
    # filter_pop_freq=(maf_df.raw$gnomAD_AF=="-" | as.numeric(maf_df.raw$gnomAD_AF) < gnomAD_AF_max) & 
    #   (maf_df.raw$AF=="-" | as.numeric(maf_df.raw$AF) < AF_max) &
    #   (maf_df.raw$ExAC_AF=="-" | as.numeric(maf_df.raw$ExAC_AF) < ExAC_AF_max) &
    #   (maf_df.raw$tumor_freq > tumor_freq_min)
    options(warn=-1)
    filter_tumor_depth=as.numeric(maf_df.raw$t_depth) > t_depth_min
    filter_norm_alt=maf_df.raw$norm_freq < norm_freq_max
    filter_tumor_alt=maf_df.raw$tumor_freq > tumor_freq_min
    filter_genes=!maf_df.raw$Hugo_Symbol %in% flag_genes
    filter_pop_freq=(maf_df.raw$gnomAD_AF %in% c("-","") | is.na(maf_df.raw$gnomAD_AF) | as.numeric(maf_df.raw$gnomAD_AF) < min(gnomAD_AF_max,1)) & 
      (maf_df.raw$AF %in% c("-","") | is.na(maf_df.raw$AF)  | as.numeric(maf_df.raw$AF) < min(AF_max,1)) &
      (maf_df.raw$ExAC_AF %in% c("-","") | is.na(maf_df.raw$ExAC_AF) | as.numeric(maf_df.raw$ExAC_AF) < min(ExAC_AF_max,1))
    options(warn=0)
  }
  filter_caller=rep(TRUE,nrow(maf_df.raw))
  if (! is.null(variant_caller)) {       ### Set 'variant_caller' to NULL to skip any filtering based on caller
    maf_df.raw$set[maf_df.raw$set=="" & maf_df.raw$Hugo_Symbol=="Hugo_Symbol"] <- "set"
    maf_df.raw$set[maf_df.raw$set==""] <- "N.A."
    if (variant_caller == "consensus") {   ### Set 'variant_caller' to 'consensus' to keep variants by two or more callers
      # filter_caller <- grepl("-|Intersection", maf_df.raw$set)
      filter_caller <- unlist(lapply(strsplit(maf_df.raw$set,"-"), function(x) {length(x)>=n_callers | "Intersection" %in% x}))
    } else {                             ### Set 'variant_caller' to one of the callers (mutect, mutect2, vardict, or strelka) to get only that caller
      # filter_caller <- grepl(paste0(variant_caller,"[|-]|Intersection"), maf_df.raw$set)
      filter_caller <- unlist(lapply(strsplit(maf_df.raw$set,"-"), function(x) {any(c(variant_caller,"Intersection") %in% x)}))
    }
  }
  
  # maf_df.raw <- maf_df.raw[filter_genes &filter_norm_alt & filter_tumor_alt & filter_pop_freq & filter_caller,]
  maf_df.rawest <- maf_df.raw
  # maf_df.raw <- maf_df.rawest
  maf_df.raw <- maf_df.raw[filter_tumor_depth & filter_norm_alt & filter_tumor_alt & filter_genes & filter_pop_freq & filter_caller,]
  # maf_df.raw <- maf_df.raw[complete.cases(maf_df.raw),]
  
  if (!is.null(save_name)) {
    if (! dir.exists(dirname(save_name))) {
      dir.create(dirname(save_name))
    }
    write.table(maf_df.raw, sep="\t", quote=F, file = save_name, row.names = F)
    print(paste0("Saving filtered maf to ",save_name))
    return(save_name)
  } else {
    return(maf_df.raw)
  }
}


make_single_ribbon_plot_alt <- function(maf, onco_genes, save_name=NULL, ribbon_color=NULL, 
                                    pval_low=0.05, pval_high=0.1,
                                    plot_frac_mut_axis=TRUE,
                                    scale_ribbon_to_fracmut=TRUE) {
  # pval_low <- 0.05
  # browser()
  if (!is.null(save_name)) {
    if (! dir.exists(dirname(save_name))) {
      dir.create(dirname(save_name))
    }
    plot_file <- gsub(".pdf",".interactions.pdf",save_name)
    pdf(file = plot_file,height=5,width=5)
  }
  som_int <-  somaticInteractions(maf = maf, genes=onco_genes, pvalue = c(pval_low, pval_high), kMax=5)
  
  if (!is.null(save_name)) {
    dev.off()
  }
  # browser()
  cooccur_data <- som_int$pairs
  cooccur_data$pair_string <- apply(cooccur_data[,1:2], 1, function(x) {paste0(sort(x), collapse="_")})
  # cooccur_data$popfrac <- cooccur_data[,"11"]/rowSums(cooccur_data[,c("00","11","01","10")], na.rm = T)  ## Frac of observed mutations
  # cooccur_data$popfrac <- cooccur_data[,"11"]/as.numeric(my_maf@summary$summary[3])   ## Frac of all samples in MAF
  cooccur_data$popfrac <- NA
  cooccur_idx <- cooccur_data$Event=="Co_Occurence"
  mut_excl_idx = cooccur_data$Event=="Mutually_Exclusive"
  
  if (scale_ribbon_to_fracmut) {
    cooccur_data$popfrac1[cooccur_idx] <- unlist(cooccur_data[cooccur_idx,"11"]/as.numeric(my_maf@summary$summary[3]))
    cooccur_data$popfrac2[cooccur_idx] <- cooccur_data$popfrac1[cooccur_idx]
    cooccur_data$popfrac1[mut_excl_idx] <- unlist(cooccur_data[mut_excl_idx,"10"]/as.numeric(my_maf@summary$summary[3]))
    cooccur_data$popfrac2[mut_excl_idx] <- unlist(cooccur_data[mut_excl_idx,"01"]/as.numeric(my_maf@summary$summary[3]))
  } else {
    cooccur_data$popfrac1 <- 1
    cooccur_data$popfrac2 <- 1
  }
  # chord_data <- cooccur_data[cooccur_data$Event=="Co_Occurance",c("gene1","gene2","popfrac","pValue")]
  chord_data <- cooccur_data[,c("gene1","gene2","popfrac1","popfrac2","pValue","Event")]
  chord_data[which(is.na(chord_data[,3])),3] <- 0
  
  require(RColorBrewer)
  # browser()
  interacting_genes <- unique(unlist(chord_data[,1:2]))
  # genecols <- colorRampPalette(brewer.pal(8,"Accent"))(length(interacting_genes))
  # genecols <- colorRampPalette(brewer.pal(8,"Dark2"))(length(interacting_genes))
  genecols <- colorRampPalette(brewer.pal(8,"Set1"))(length(interacting_genes))
  # genecols <- rainbow((length(interacting_genes)))
  names(genecols) <- interacting_genes
  
  sig_line_type=2
  link_border <- ifelse(chord_data$pValue < pval_low, sig_line_type,0)
  # color_legend <- Legend(labels=gsub("_"," ",names(cluster_sig_colors)), 
  #                        legend_gp = gpar(fill = cluster_sig_colors),background = "white",size=unit(0.08,"npc"),pch=22,
  #                        type="points",direction="horizontal",nrow=1)
  line_legend <- Legend(labels=paste0("p < ",pval_low), 
                        labels_gp = gpar(fontsize = 16, fontface = "bold"),
                        legend_gp = gpar(lty=sig_line_type,lwd = 1, fontsize=16),background = "white",
                        type="lines",direction="horizontal",nrow=1)
  # mylegend <- packLegend(color_legend, line_legend)
  if (!is.null(save_name)) {
    if (! dir.exists(dirname(save_name))) {
      dir.create(dirname(save_name))
    }
    pdf(file = save_name,height=7,width=7)
  }
  
  if (is.null(ribbon_color)) {
    # ribbon_color = genecols[match(names(genecols),c(chord_data$gene1,chord_data$gene2))]
    ribbon_color = genecols[chord_data$gene1]
    # ribbon_color[chord_data$Event == "Mutually_Exclusive"] = "grey90"
    ribbon_color[chord_data$Event == "Mutually_Exclusive"] = "white"
  }
  link_border[chord_data$Event == "Mutually_Exclusive"]=1
  shrink_factor=1.3 # Higher = more shrinkage; use to control whitespace (or lack thereof) around figure 
  circos.clear()
  circos.par(canvas.xlim=c(-shrink_factor,shrink_factor),canvas.ylim=c(-shrink_factor,shrink_factor))
  chordDiagram(chord_data[,1:4],grid.col = genecols,
               annotationTrack = c("grid",ifelse(plot_frac_mut_axis, "axis", "")),
               col=ribbon_color,
               # transparency = link_alpha, 
               link.lty = link_border, 
               link.border = "black",
               link.sort = TRUE)
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                facing = "clockwise", niceFacing = TRUE, adj = c(-0.5, 0.5))
  }, bg.border = NA) # here set bg.border to NA is important
  # draw(color_legend, x = unit(0.5, "npc"), y = unit(0.02, "npc"), just = c("center"))
  draw(line_legend, x = unit(0.5, "npc"), y = unit(0.97, "npc"), just = c("center"))
  
  if (!is.null(save_name)) {
    dev.off()
  }
  
  
}

## Makes the colors used for annotations
my_oncoplot_colors <- function(hm_anno_info) {
  require(circlize)
  require(RColorBrewer)
  cluster_colors = c(MO1="#4DAF4A",MO2="#377EB8",MO3="orange",MO4="#E41A1C","NA"="grey70")
  sex_colors <- c("M"="blue","F"="hotpink","NA"="white")
  tnm_colors <- c(brewer.pal(4, "Reds"),"grey70")
  names(tnm_colors) <- c("1","2","3","4","NA")
  
  def_bin_colors <- c("0" = "gray70","1" = "grey20","NA"="white")
  
  all_colors <- rep(list(def_bin_colors),ncol(hm_anno_info))
  names(all_colors) <- colnames(hm_anno_info)
  all_colors[[grep("sex",colnames(hm_anno_info), ignore.case = T)]] <- sex_colors
  all_colors[[grep("class",colnames(hm_anno_info), ignore.case = T)]] <- cluster_colors
  # all_colors[[grep("TNM_staging",colnames(hm_anno_info), ignore.case = T)]] <- tnm_colors
  
  return(all_colors)
}


### Cretaes matrix for oncoplot from maf file
### Adapted from maftools: https://github.com/PoisonAlien/maftools/blob/master/R/oncomatrix.R
createOncoMatrix = function(m, g = NULL, chatty = TRUE, add_missing = FALSE){
  
  if(is.null(g)){
    stop("Please provde atleast two genes!")
  }
  
  subMaf = subsetMaf(maf = m, genes = g, includeSyn = FALSE, mafObj = FALSE)
  
  if(nrow(subMaf) == 0){
    if(add_missing){
      numericMatrix = matrix(data = 0, nrow = length(g), ncol = length(levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])))
      rownames(numericMatrix) = g
      colnames(numericMatrix) = levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])
      
      oncoMatrix = matrix(data = "", nrow = length(g), ncol = length(levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])))
      rownames(oncoMatrix) = g
      colnames(oncoMatrix) = levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])
      
      vc = c("")
      names(vc) = 0
      
      return(list(oncoMatrix = oncoMatrix, numericMatrix = numericMatrix, vc = vc))
    }else{
      return(NULL)
    }
  }
  
  if(add_missing){
    subMaf[, Hugo_Symbol := factor(x = Hugo_Symbol, levels = g)]
  }
  
  oncomat = data.table::dcast(data = subMaf[,.(Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode)], formula = Hugo_Symbol ~ Tumor_Sample_Barcode,
                              fun.aggregate = function(x){
                                x = unique(as.character(x))
                                xad = x[x %in% c('Amp', 'Del')]
                                xvc = x[!x %in% c('Amp', 'Del')]
                                
                                if(length(xvc)>0){
                                  xvc = ifelse(test = length(xvc) > 1, yes = 'Multi_Hit', no = xvc)
                                }
                                
                                x = ifelse(test = length(xad) > 0, yes = paste(xad, xvc, sep = ';'), no = xvc)
                                x = gsub(pattern = ';$', replacement = '', x = x)
                                x = gsub(pattern = '^;', replacement = '', x = x)
                                return(x)
                              } , value.var = 'Variant_Classification', fill = '', drop = FALSE)
  
  #convert to matrix
  data.table::setDF(oncomat)
  rownames(oncomat) = oncomat$Hugo_Symbol
  oncomat = as.matrix(oncomat[,-1, drop = FALSE])
  
  variant.classes = as.character(unique(subMaf[,Variant_Classification]))
  variant.classes = c('',variant.classes, 'Multi_Hit')
  names(variant.classes) = 0:(length(variant.classes)-1)
  
  #Complex variant classes will be assigned a single integer.
  vc.onc = unique(unlist(apply(oncomat, 2, unique)))
  vc.onc = vc.onc[!vc.onc %in% names(variant.classes)]
  names(vc.onc) = rep(as.character(as.numeric(names(variant.classes)[length(variant.classes)])+1), length(vc.onc))
  variant.classes2 = c(variant.classes, vc.onc)
  
  oncomat.copy <- oncomat
  #Make a numeric coded matrix
  for(i in 1:length(variant.classes2)){
    oncomat[oncomat == variant.classes2[i]] = names(variant.classes2)[i]
  }
  
  #If maf has only one gene
  if(nrow(oncomat) == 1){
    mdf  = t(matrix(as.numeric(oncomat)))
    rownames(mdf) = rownames(oncomat)
    colnames(mdf) = colnames(oncomat)
    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  }
  
  #convert from character to numeric
  mdf = as.matrix(apply(oncomat, 2, function(x) as.numeric(as.character(x))))
  rownames(mdf) = rownames(oncomat.copy)
  
  
  #If MAF file contains a single sample, simple sorting is enuf.
  if(ncol(mdf) == 1){
    sampleId = colnames(mdf)
    mdf = as.matrix(mdf[order(mdf, decreasing = TRUE),])
    colnames(mdf) = sampleId
    
    oncomat.copy = as.matrix(oncomat.copy[rownames(mdf),])
    colnames(oncomat.copy) = sampleId
    
    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  } else{
    #Sort by rows as well columns if >1 samples present in MAF
    #Add total variants per gene
    mdf = cbind(mdf, variants = apply(mdf, 1, function(x) {
      length(x[x != "0"])
    }))
    #Sort by total variants
    mdf = mdf[order(mdf[, ncol(mdf)], decreasing = TRUE), ]
    #colnames(mdf) = gsub(pattern = "^X", replacement = "", colnames(mdf))
    nMut = mdf[, ncol(mdf)]
    
    mdf = mdf[, -ncol(mdf)]
    
    mdf.temp.copy = mdf #temp copy of original unsorted numeric coded matrix
    
    mdf[mdf != 0] = 1 #replacing all non-zero integers with 1 improves sorting (& grouping)
    tmdf = t(mdf) #transposematrix
    mdf = t(tmdf[do.call(order, c(as.list(as.data.frame(tmdf)), decreasing = TRUE)), ]) #sort
    
    mdf.temp.copy = mdf.temp.copy[rownames(mdf),] #organise original matrix into sorted matrix
    mdf.temp.copy = mdf.temp.copy[,colnames(mdf)]
    mdf = mdf.temp.copy
    
    #organise original character matrix into sorted matrix
    oncomat.copy <- oncomat.copy[,colnames(mdf)]
    oncomat.copy <- oncomat.copy[rownames(mdf),]
    
    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  }
}

### Ordering function to make rainfall pattern
orderByGroup <- function(my_oncomat, info_vec) {
  class2count <- function(class_str) {
    returnval = nchar(as.character(class_str))
    if (is.na(returnval)) {
      returnval=0
    } else {
      returnval=ifelse(returnval > 2, 1, 0)
    }
    return(returnval)
  }
  oncoprint_column_order = function(count_matrix) {
    scoreCol = function(x) {
      score = 0
      for(i in 1:length(x)) {
        if(x[i]) {
          score = score + 2^(length(x)-i*1/x[i])
          score=min(score,1e16)
        }
      }
      return(score)
    }
    scores = apply(count_matrix, 2, scoreCol)
    new_order = order(scores, decreasing=TRUE)
    return(new_order)
  }
  
  # browser()
  oncomat_bin <- apply(my_oncomat, 1:2, class2count)
  new_names_order <- c()
  prev_max_idx=0
  for (curr_level in levels(info_vec)) {
    subset_ind=which(as.character(info_vec)==curr_level)
    if (length(subset_ind) < 2) {
      curr_new_names_order <- colnames(oncomat_bin)[subset_ind]
    } else {
      oncdat=oncomat_bin[,subset_ind]
      curr_new_names_order <- colnames(oncdat)[oncoprint_column_order(oncdat)]
    }
    new_names_order <- c(new_names_order, curr_new_names_order)
  }
  new_order <- match(new_names_order, colnames(my_oncomat))
  return(new_order)
}



make_single_ribbon_plot <- function(maf, onco_genes, save_name=NULL, ribbon_color=NULL, 
                                    pval_low=0.05, pval_high=0.1,
                                    plot_frac_mut_axis=TRUE,
                                    shrink_factor=1.3, # Higher = more shrinkage; use to control whitespace (or lack thereof) around figure 
                                    scale_ribbon_to_fracmut=TRUE) {
  # pval_low <- 0.05
  # browser()
  if (!is.null(save_name)) {
    if (! dir.exists(dirname(save_name))) {
      dir.create(dirname(save_name))
    }
    plot_file <- gsub(".pdf",".interactions.pdf",save_name)
    pdf(file = plot_file,height=5,width=5)
  }
  som_int <-  somaticInteractions(maf = maf, genes=onco_genes, pvalue = c(pval_low, pval_high), kMax=5)
  
  if (!is.null(save_name)) {
    dev.off()
  }
  # browser()
  cooccur_data <- som_int$pairs
  cooccur_data$pair_string <- apply(cooccur_data[,1:2], 1, function(x) {paste0(sort(x), collapse="_")})
  # cooccur_data$popfrac <- cooccur_data[,"11"]/rowSums(cooccur_data[,c("00","11","01","10")], na.rm = T)  ## Frac of observed mutations
  # cooccur_data$popfrac <- cooccur_data[,"11"]/as.numeric(my_maf@summary$summary[3])   ## Frac of all samples in MAF
  cooccur_data$popfrac <- NA
  cooccur_idx <- cooccur_data$Event=="Co_Occurence"
  mut_excl_idx = cooccur_data$Event=="Mutually_Exclusive"
  
  if (scale_ribbon_to_fracmut) {
    cooccur_data$popfrac1[cooccur_idx] <- unlist(cooccur_data[cooccur_idx,"11"]/as.numeric(my_maf@summary$summary[3]))
    cooccur_data$popfrac2[cooccur_idx] <- cooccur_data$popfrac1[cooccur_idx]
    cooccur_data$popfrac1[mut_excl_idx] <- unlist(cooccur_data[mut_excl_idx,"10"]/as.numeric(my_maf@summary$summary[3]))
    cooccur_data$popfrac2[mut_excl_idx] <- unlist(cooccur_data[mut_excl_idx,"01"]/as.numeric(my_maf@summary$summary[3]))
  } else {
    cooccur_data$popfrac1 <- 1
    cooccur_data$popfrac2 <- 1
  }
  # chord_data <- cooccur_data[cooccur_data$Event=="Co_Occurance",c("gene1","gene2","popfrac","pValue")]
  chord_data <- cooccur_data[,c("gene1","gene2","popfrac1","popfrac2","pValue","Event")]
  chord_data[which(is.na(chord_data[,3])),3] <- 0
  
  require(RColorBrewer)
  # browser()
  interacting_genes <- unique(unlist(chord_data[,1:2]))
  # genecols <- colorRampPalette(brewer.pal(8,"Accent"))(length(interacting_genes))
  # genecols <- colorRampPalette(brewer.pal(8,"Dark2"))(length(interacting_genes))
  genecols <- colorRampPalette(brewer.pal(8,"Set1"))(length(interacting_genes))
  # genecols <- rainbow((length(interacting_genes)))
  names(genecols) <- interacting_genes
  
  sig_line_type=2
  link_border <- ifelse(chord_data$pValue < pval_low, sig_line_type,0)
  # color_legend <- Legend(labels=gsub("_"," ",names(cluster_sig_colors)), 
  #                        legend_gp = gpar(fill = cluster_sig_colors),background = "white",size=unit(0.08,"npc"),pch=22,
  #                        type="points",direction="horizontal",nrow=1)
  line_legend <- Legend(labels=paste0("p < ",pval_low), 
                        labels_gp = gpar(fontsize = 16, fontface = "bold"),
                        legend_gp = gpar(lty=sig_line_type,lwd = 1, fontsize=16),background = "white",
                        type="lines",direction="horizontal",nrow=1)
  # mylegend <- packLegend(color_legend, line_legend)
  if (!is.null(save_name)) {
    if (! dir.exists(dirname(save_name))) {
      dir.create(dirname(save_name))
    }
    pdf(file = save_name,height=7,width=7)
  }
  
  if (is.null(ribbon_color)) {
    # ribbon_color = genecols[match(names(genecols),c(chord_data$gene1,chord_data$gene2))]
    ribbon_color = genecols[chord_data$gene1]
    # ribbon_color[chord_data$Event == "Mutually_Exclusive"] = "grey90"
    ribbon_color[chord_data$Event == "Mutually_Exclusive"] = "white"
  }
  link_border[chord_data$Event == "Mutually_Exclusive"]=1
  # shrink_factor=1.3 # Higher = more shrinkage; use to control whitespace (or lack thereof) around figure 
  circos.clear()
  circos.par(canvas.xlim=c(-shrink_factor,shrink_factor),canvas.ylim=c(-shrink_factor,shrink_factor))
  chordDiagram(chord_data[,1:4],grid.col = genecols,
               annotationTrack = c("grid",ifelse(plot_frac_mut_axis, "axis", "")),
               col=ribbon_color,
               # transparency = link_alpha, 
               link.lty = link_border, 
               link.border = "black",
               link.sort = TRUE)
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                facing = "clockwise", niceFacing = TRUE, adj = c(-0.5, 0.5))
  }, bg.border = NA) # here set bg.border to NA is important
  # draw(color_legend, x = unit(0.5, "npc"), y = unit(0.02, "npc"), just = c("center"))
  draw(line_legend, x = unit(0.5, "npc"), y = unit(0.97, "npc"), just = c("center"))
  
  if (!is.null(save_name)) {
    dev.off()
  }
  
  
}


make_lolliplot_custom <- function(maf, onco_genes, save_dir=NULL, save_height=6, save_width=5,domain_db = "pfam") {
  library(maftools)
  library(trackViewer)
  library(biomaRt)
  library(dplyr)
  library(trackViewer)
  library(RColorBrewer)
  
  if (!is.null(save_dir)) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = T)
    }
  }
  # browser()
  ## Collect MAF data for oncogenes
  maf.oncoplot <- subsetMaf(maf, genes=onco_genes)
  mafdata <- as.data.frame(maf.oncoplot@data)
  myprot <- sort(unique(mafdata$UNIPARC))  ## Collects the protein IDs
  
  ## Ok, so first, grab the peptide sequence from biomart
  ## This is just to establish how long the plot needs to be; I couldn't find protein length information summarized anywhere
  ## Also this is from the "sequences" section of Biomart, so can't retrieve at the same time as the other stuff we need
  ensembl = useEnsembl("ensembl",dataset = "hsapiens_gene_ensembl", mirror = "uswest")
  protInfo.pep = getBM(attributes = c("chromosome_name","start_position","end_position","hgnc_symbol", "uniparc","peptide"),
                       filters = "uniparc", values = myprot, mart = ensembl, verbose=F)
  protInfo.pep$peptide_length <- nchar(protInfo.pep$peptide)
  
  ## Now we're going to do some text parsing on the protein-related MAF fields to extract protein position info
  mafdata$Chromosome <- gsub("chr","",mafdata$Chromosome)
  # This gets the amina acid positions of the variants from the MAF column; this fails sometimes, like for splice variants because the position is annotated as "-"
  # mafdata$aa_pos <- as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(mafdata$Protein_position,"/"),"[[",1)),"-"),"[[",1)))
  # This gets the amina acid positions of the variants from the HGVS descriptor string; seems to be comprehensive (i.e. always numerical)
  mafdata$hgvs_pos <- as.numeric(gsub("^p\\.[[:alpha:]]{1}([[:digit:]]+)[[:alpha:]|[:punct:]]{1,}[[:graph:]]*$","\\1",mafdata$HGVSp_Short))
  # This gets the protein length from the maf column; also fails for splice variants
  # mafdata$aa_len <- as.numeric(unlist(lapply(strsplit(mafdata$Protein_position,"/"),function(x){ifelse(length(x) < 2, NA, x[2])})))
  # browser()
  ## Prepare the data for turning it into a GRanges object
  maf_gr_data <- data.frame(gene_chr=mafdata$Chromosome, gene_start=mafdata$Start_Position, gene_end=mafdata$End_Position, gene_strand=mafdata$Strand, 
                            gene=mafdata$Hugo_Symbol, aa_change=mafdata$HGVSp_Short, protein_id=mafdata$UNIPARC,
                            exon=mafdata$Exon_Number,
                            # aa_pos=mafdata$aa_pos, aa_pos_extract=mafdata$hgvs_pos, 
                            # aa_len=mafdata$aa_len, aa_len_biomart=protInfo.pep$peptide_length[match(mafdata$UNIPARC, protInfo.pep$uniparc)]-1,
                            aa_pos_extract=mafdata$hgvs_pos, 
                            aa_len_biomart=protInfo.pep$peptide_length[match(mafdata$UNIPARC, protInfo.pep$uniparc)]-1,
                            stringsAsFactors = F)
  maf_gr_data <- cbind(maf_gr_data,
                       var_start_pos = maf_gr_data$aa_pos_extract,
                       var_end_pos = maf_gr_data$aa_pos_extract+1,
                       protein_start = 1,
                       protein_end = maf_gr_data$aa_len_biomart,
                       impact=mafdata$IMPACT,
                       variant_type=mafdata$Variant_Classification,
                       AF_1kG=ifelse(mafdata$AF=="-",0,as.numeric(mafdata$AF)),
                       AF_ExAC=ifelse(mafdata$ExAC_AF=="-",0,as.numeric(mafdata$ExAC_AF))                     )
  
  ## Calculate the number of samples in which each unique variant is mutated
  maf_gr_data <-maf_gr_data %>% dplyr::group_by(aa_change, variant_type) %>% dplyr::mutate(n_altered=as.double(dplyr::n()))
  
  ## Make GRanges object for the variants
  mafdata.gr <- makeGRangesFromDataFrame(maf_gr_data, keep.extra.columns = T,
                                         seqnames.field = "protein_id", start.field = "var_start_pos", end.field = "var_end_pos")
  ## Make GRanges object for the proteins
  proteins.gr <- unique(makeGRangesFromDataFrame(maf_gr_data, keep.extra.columns = F,
                                                 seqinfo = Seqinfo(seqnames = protInfo.pep$uniparc, seqlengths = protInfo.pep$peptide_length),
                                                 seqnames.field = "protein_id", start.field = "var_start_pos", end.field = "var_end_pos"))
  
  ### Now collect protein domain information
  if (domain_db == "pfam") {
    domain_attributes <- c("pfam", "pfam_start","pfam_end")
    library(PFAM.db)
    pfam_descriptions <- as.list(PFAMDE[mappedkeys(PFAMDE)])
  } else if (domain_db == "interpro") {
    domain_attributes <- c("interpro_short_description", "interpro_start","interpro_end")
  } else {
    stop("domain_db must be 'pfam' or 'interpro'")
  }
  
  protInfo.interpro = getBM(attributes = c("chromosome_name","start_position","end_position","hgnc_symbol", "uniparc",
                                           domain_attributes),
                            # "interpro","interpro_short_description","interpro_description","interpro_start","interpro_end"),
                            filters = "uniparc", values = seqnames(proteins.gr), mart = ensembl, verbose=F)
  
  interpro.gr <- makeGRangesFromDataFrame(protInfo.interpro,seqnames.field = "uniparc",
                                          seqinfo = seqinfo(proteins.gr),
                                          # start.field = "interpro_start",end.field = "interpro_end",keep.extra.columns = T)
                                          start.field = paste0(domain_db,"_start"),end.field = paste0(domain_db,"_end"),keep.extra.columns = T)
  
  if (domain_db == "pfam") {
    interpro.gr$domain_label <- unlist(pfam_descriptions[match(interpro.gr$pfam, names(pfam_descriptions))])
  } else if (domain_db == "interpro") {
    interpro.gr$domain_label <- mcols(interpro.gr)[, domain_attributes[1]]
  } else {
    stop("domain_db must be 'pfam' or 'interpro'")
  }
  
  ## Make some colors
  variant_type_colors <- c(Nonsense_Mutation="#ad7aff",Missense_Mutation="#377EB8",Frame_Shift_Del="#4DAF4A",
                           In_Frame_Ins="#e7d9ff",Splice_Site="#FF7F00",Multi_Hit="#FFFF33",Frame_Shift_Ins="#A65628",
                           In_Frame_Del="#f781bf",Translation_Start_Site="#400085",Nonstop_Mutation="#a5d8fa")
  names(variant_type_colors) <- gsub("_"," ",names(variant_type_colors))
  
  
  ### Finally, loop through the genes and make a lolliplot
  mygenes <- sort(unique(mafdata.gr$gene))
  for (curr_gene in mygenes) {
    # curr_gene="FKBP9"
    print(paste0(sum(mafdata.gr$gene == curr_gene), " of ", length(mafdata.gr)," variants fall within gene ", curr_gene))
    
    sample.gr <- sort(mafdata.gr[mafdata.gr$gene == curr_gene])
    sample.gr <- sample.gr[!duplicated(mcols(sample.gr)[,c("aa_change")])]
    
    names(sample.gr) <- gsub("p\\.","",paste0(sample.gr$aa_change," (",sample.gr$exon,")"))
    sample.gr$color <- variant_type_colors[match(gsub("_"," ",sample.gr$variant_type),names(variant_type_colors))]
    sample.gr$score <- sample.gr$n_altered + 1e-12
    sample.gr$label.parameter.rot <- 45
    sample.gr$dashline.col <- sample.gr$color
    sample.gr$label.parameter.gp <- gpar(cex=0.8)
    
    end(sample.gr)[width(sample.gr) > 1] <- start(sample.gr)[width(sample.gr) > 1]
    features <- interpro.gr[interpro.gr$hgnc_symbol==curr_gene]
    features <- keepSeqlevels(features, unique(as.character(seqnames(features))))
    features <- c(features, gaps(features))
    features <- features[strand(features)=="*"]
    features$domain_label[is.na(features$domain_label)] <- ""
    names(features) <- features$domain_label
    
    # domain_colors <- c(rainbow(length(unique(features$domain_label))-1, alpha = 1),"grey70")
    domain_colors <- c(rainbow(length(unique(features$domain_label))-1, alpha = 1),"grey70")
    names(domain_colors) <- c(sort(unique(features$domain_label))[-1],"")
    
    features$fill <- domain_colors[match(features$domain_label, names(domain_colors))]
    features$height <- 0.05
    
    # browser()
    legend_var_colors <- variant_type_colors[names(variant_type_colors) %in% gsub("_"," ",sample.gr$variant_type)]
    mylegend <- list(list(labels=names(legend_var_colors), fill=legend_var_colors))
    
    dist_from_end=max(end(sample.gr))/seqlengths(features)
    # padding=2^(scales::rescale(dist_from_end, to=c(5,10), from=c(0.7,1)))
    padding=log2(scales::rescale(dist_from_end, to=c(2,2.3), from=c(0,1)))
    padded_range <- GRanges(levels(seqnames(features)), IRanges(1,seqlengths(features) * padding))
    
    xaxis <- seq(round(min(start(features))-50,-2), round(max(end(features))+50,-2), length.out = 5)
    
    if (!is.null(save_dir)) {
      save_name=paste0(save_dir,"/",curr_gene,".pdf")
      pdf(save_name,width=save_width,height=save_height)
    }
    
    lolliplot(sample.gr, features, legend=mylegend,ylab="Number of Mutated Samples", 
              type="circle", ranges = padded_range,
              xaxis=xaxis)
    grid.text(paste0("Mutations in ",curr_gene," gene\n"), x=.5, y=.98, just="top", 
              gp=gpar(cex=1.5, fontface="bold"))
    
    if (!is.null(save_dir)) {
      dev.off()
    }
  }  
  
  
}

hclust_semisupervised <- function(data, groups, dist_method = "euclidean",
                                  dist_p = 2, hclust_method = "complete") {
  ## From: https://stackoverflow.com/questions/49541604/merging-multiple-hclust-objects-or-dendrograms
  hclist <- lapply(groups, function (group) {
    hclust(dist(data[group,], method = dist_method, p = dist_p), method = hclust_method)
  })
  d <- as.dendrogram(hclist[[1]])
  for (i in 2:length(hclist)) {
    d <- merge(d, as.dendrogram(hclist[[i]]))
  }
  hc <- as.hclust(d)
  data_reordered <- data[unlist(groups),]
  
  return(list(data = data_reordered, hclust = hc))
}



mutation_colors <- c(Nonsense_Mutation="#ad7aff",Missense_Mutation="#377EB8",Frame_Shift_Del="#4DAF4A",
         In_Frame_Ins="#ff008c",Splice_Site="#FF7F00",Multi_Hit="#FFFF33",Frame_Shift_Ins="#A65628",
         In_Frame_Del="#f781bf",Translation_Start_Site="#400085",Nonstop_Mutation="#b68dfc",
         no_variants="#d6d6d6")

### List defining functions for color and shape of cells in oncoplot
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # "0" = function(x, y, w, h) {
  #   grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
  #             gp = gpar(fill = "#CCCCCC", col = NA))
  # },
  "Nonsense Mutation" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Nonsense Mutation"], col = NA))
  },
  "Missense Mutation" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Missense Mutation"], col = NA))
  },
  "Frame Shift Del" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Frame Shift Del"], col = NA))
  },
  "In Frame Ins" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["In Frame Ins"], col = NA))
  },
  "Splice Site" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = col["Splice Site"], col = NA))
  },
  "Multi Hit" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = col["Multi Hit"], col = NA))
  },
  "Frame Shift Ins" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = col["Frame Shift Ins"], col = NA))
  },
  "In Frame Del" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["In Frame Del"], col = NA))
  },
  "Nonstop Mutation" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = col["Nonstop Mutation"], col = NA))
  },
  "Translation Start Site" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = col["Translation Start Site"], col = NA))
  },
  "no variants" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              # gp = gpar(fill = "#e0e0e0", col = NA))
              gp = gpar(fill = "#CCCCCC", col = NA))
  }
)
