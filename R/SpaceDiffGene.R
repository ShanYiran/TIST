#' Title
#'
#' @param SC_obj load RNA-seq use an existed Seruat object
#' @param Spot_manifest The TIST results returned by function "Meta_St_img_unsupervised"
#' @param savePath The address where the results stored.
#' @param MC_obj_file MC obj
#' @param Mc_manifest Mc_manifest
#' @param Mc_pair_SC Mc_pair_SC
#' @param methods TIST network method is "Walktrap_id"
#' @param SPARK_file SPARK_file generate by SPARK
#' @param MC_ident_file MC_ident_file
#' @param DEgeneplot DEgeneplot
#' @param Score_show T/F
#' @param spark_methods SPARK/SPARK-x
#' @param netfile The TIST-net build by function "Meta_St_img_unsupervised"
#'
#' @return MD_markers
#'
SpaceDiffGene <- function(SC_obj,
                          Spot_manifest,
                          savePath,
                          MC_obj_file=NULL,
                          Mc_manifest=NULL,
                          Mc_pair_SC,
                          methods= "walktrap",
                          SPARK_file,
                          MC_ident_file = NULL,
                          DEgeneplot = F,
                          Score_show = T,
                          spark_methods = NULL,
                          netfile = NULL){
  if(methods=="walktrap"){
    Idents(SC_obj) <- Spot_manifest$Walktrap_id
  }
  if(is.null(MC_obj_file)==F){
    MC_obj <- readRDS(MC_obj_file)
    mc_matrix <- MC_obj@assays$RNA@scale.data
  }
  MD_markers <- FindAllMarkers(SC_obj,features = rownames(SC_obj))
  saveRDS(MD_markers,file = paste0(savePath,"TIST_markers.RDS"))
  write.csv(MD_markers,file = paste0(savePath,"TIST_markers.csv"))
  #MDMarkers gene vs high var genes
  #MD_gene <- MD_markers %>% group_by(cluster) %>% top_n(n=30,wt=avg_logFC)
  #MD_gene <- MD_markers[abs(MD_markers$avg_logFC)>=1,]
  MD_gene <- MD_markers %>% top_n(n=150,wt=avg_log2FC)
  if(!is.null(spark_methods)){
    SPARK_gene <-readRDS(SPARK_file)
    if(spark_methods=="SPARK"){
      features3 <- SPARK_gene@res_mtest
      features3<- features3[sort(features3$adjusted_pvalue,index.return=TRUE)$ix,]
      #features3<- features3[sort(features3$adjusted_pvalue,index.return=TRUE)$ix,]
      features3 <- rownames(features3)
    }
    if(spark_methods=="SPARKX"){
      features3 <- SPARK_gene$res_mtest
      features3 <- features3[sort(features3$adjustedPval,index.return=TRUE)$ix,]
      #features3<- features3[sort(features3$adjusted_pvalue,index.return=TRUE)$ix,]
      features3 <- rownames(features3)
    }
    venn.plot <- venn.diagram(
      x = list(
        MD =unique(MD_gene$gene)[1:100],
        spark = features3[1:100]
      ),
      filename = paste0(savePath,"gene_vs_x.jpeg"),
      cex = 2.5,
      cat.cex = 2.5,
      cat.pos = c(-20, 20),
      ext.line.lty = "dotted",
      ext.line.lwd = 3,
      ext.pos = 12,
      ext.dist = -0.12,
      ext.length = 0.85
    )
    if(Score_show){
      SD_score <- c()
      SPARK_score <- c()
      SD_gene <- unique(MD_gene$gene)[1:100]
      Spark_gene <- features3[1:100]
      ST_net <- readRDS(file = netfile)
      KNN <- as_adjacency_matrix(ST_net)
      KNN <- as.matrix(KNN)
      SC_obj <- ScaleData(SC_obj,features = c(VariableFeatures(SC_obj),SD_gene,Spark_gene))
      colnames(SC_obj@assays$RNA@scale.data) <- ST_filter_str(colnames(SC_obj@assays$RNA@scale.data),'-')
      for(i in 1:100){
        SD_gene_expr <- SC_obj@assays$RNA@scale.data[SD_gene[i],colnames(KNN)]
        SD_score <- c(SD_score,as.numeric(t(SD_gene_expr)%*%KNN%*%(SD_gene_expr)))
        SPARK_gene_expr <- SC_obj@assays$RNA@scale.data[Spark_gene[i],colnames(KNN)]
        SPARK_score <- c(SPARK_score,as.numeric(t(SPARK_gene_expr)%*%KNN%*%(SPARK_gene_expr)))
      }
      #SD_score <- (SD_score-min(SD_score))/(max(SD_score)-min(SD_score))
      #SPARK_score <- (SPARK_score-min(SPARK_score))/(max(SPARK_score)-min(SPARK_score))
      names(SD_score) <- SD_gene
      names(SPARK_score) <- Spark_gene
      saveRDS(SD_score,paste0(savePath,"SD_score.RDS"))
      saveRDS(SPARK_score,paste0(savePath,"SPARK_score.RDS"))

      x_o <- SD_score
      x_e <- SPARK_score
      score <- c(x_o,x_e)
      lable <- c(rep("SD",length(x_o)),rep("SPARK",length(x_e)))

      df <- data.frame(score = score, lable = lable)

      p <- gghistogram(df,x = "score", add = "mean", rug = TRUE, fill = "lable",palette = c("#00AFBB","#E7B800"),bins = 30)
      ggsave(filename = paste0(savePath,"SDvsSPARK.pdf"),
             p, width = 6, height = 5, dpi = 150, limitsize = FALSE)

      gl <- intersect(SD_gene,Spark_gene)
      feature_plot(gene_plot = gl,SC_obj, Spot_manifest,savePath,savefilename = "UNI_gene_plot",height = 550,pointsize = 2.5)
      mc_marker_plot(gl = gl,MC_obj_file,Mc_manifest,Spot_manifest,savePath,savefilename = "UNI_gene_MC_plot")
    }
  }

  if(DEgeneplot){
    all.markers <- MD_markers
    Mc_pair_SC <- readRDS(file = MC_ident_file)
    names(Mc_pair_SC) <- ST_filter_str(names(Mc_pair_SC),"-")
    #DoHeatmap(SC_obj, features = top5$gene) + NoLegend()
    Spot_manifest <- cbind(Spot_manifest,Mc_pair_SC)

    colnames(Spot_manifest)[8] <- "mc"
    markers_MC <- all.markers
    top20_mc <- markers_MC %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
    top5_mc <- markers_MC %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
    top20_mc <- top20_mc$gene
    top5_mc <- top5_mc$gene

    par(mar=c(10,0.5,0.5,0.5))
    cols <- array(getDefaultColors(length(table(Spot_manifest$Walktrap_id))))
    rownames(cols) <- Mc_manifest$mc_names

    mc_matrix <- MC_obj@assays$RNA@scale.data
    annotation_col = data.frame(
      #mc_type = Mc_manifest$mc_type,
      mc_name = Mc_manifest$mc_names
    )

    #annotation_col <- annotation_col[sort(annotation_col$mc_names,index.return=TRUE)$ix,]
    #annocol <- as.data.frame(annotation_col[,"mc_type"])
    annocol <- as.data.frame(annotation_col[,"mc_name"])
    rownames(annocol) <- annotation_col$mc_name
    colnames(annocol) <- "mc_name"
    annocol_mc <- annocol
    ann_colors = list(
      mc_name = cols
    )
    #MC_corr <- rcorr(mc_matrix[intersect(rownames(mc_matrix),top20_mc),rownames(annocol)])
    MC_corr <- rcorr(mc_matrix[VariableFeatures(SC_obj),rownames(annocol)])
    mc_corr <- (MC_corr$r)
    tm <- mc_corr
    diag(tm) <- 0
    mc_max <- max(tm)
    diag(mc_corr) <- mc_max
    #head(ann_colors)
    png(filename=paste0(savePath,"MC_corr.png"), width=1000, height=900)
    pheatmap(as.matrix(mc_corr), annotation_col = annocol_mc,annotation_row = annocol_mc,cluster_rows =T,cluster_cols =T,show_colnames = F,legend = T,annotation_colors = ann_colors)
    dev.off()
    mc_heat <- pheatmap(as.matrix(mc_corr), annotation_col = annocol_mc,annotation_row = annocol_mc,cluster_rows =T,cluster_cols =T,show_colnames = F,legend = F,annotation_colors = ann_colors)
    mc_cluster <- mc_heat$tree_row$labels[mc_heat$tree_row$order]


    pdf(paste0(savePath,"MC_diffgene_inte_new.pdf"), width=9, height=10)
    pheatmap(mc_matrix[top5_mc,mc_cluster], annotation_col =annocol,cluster_rows = F,cluster_cols = F,annotation_colors = ann_colors)
    dev.off()

    #dev.off()
    #mc_cluster <- annotation_col$mc_name
    #SC_diff_gene_show
    SC_barcode_to_MC <- data.frame(barcode = colnames(SC_obj),Mc_id = Mc_pair_SC)
    scrowid <- c()
    annocol_cs <- c()
    mc_id <- c()
    for(i in 1:length(mc_cluster)){
      tmidlist <- SC_barcode_to_MC$barcode[which(SC_barcode_to_MC$Mc_id==mc_cluster[i])]
      annocol_cs <- c(annocol_cs,rep(Mc_manifest$mc_type[Mc_manifest$mc_names==mc_cluster[i]],Mc_manifest$spot_num.Freq[Mc_manifest$mc_names==mc_cluster[i]]))
      scrowid <- c(scrowid,tmidlist)
      mc_id <- c(mc_id,rep(mc_cluster[i],Mc_manifest$spot_num.Freq[Mc_manifest$mc_names==mc_cluster[i]]))
    }
    annocol_cs <- data.frame(mc_id = mc_id)
    rownames(annocol_cs) <- scrowid
    ann_colors_sc = list(
      mc_id = cols
    )
    SC_obj <- ScaleData(SC_obj,features =rownames(SC_obj))
    colnames(SC_obj@assays$RNA@scale.data) <- ST_filter_str(colnames(SC_obj@assays$RNA@scale.data),'-')

    sc_matrix <- as.matrix(SC_obj@assays$RNA@scale.data)
    rownames(annocol_cs) <- ST_filter_str(rownames(annocol_cs),'-')
    show_sc_matrix <- sc_matrix[top5_mc,rownames(annocol_cs)]
    show_sc_matrix[which(show_sc_matrix>4)] <- 4
    show_sc_matrix[which(show_sc_matrix< -4)] <- -4
    num.cluster <- as.data.frame(table(mc_id))
    rownames(num.cluster) <- num.cluster[,1]
    num.cluster <- num.cluster[mc_cluster,2]
    #num.cluster <- num.cluster[as.character(1 : length(num.cluster))]
    gaps_col <- cumsum(num.cluster)

    de.pre <- list(expr.data = show_sc_matrix, spot.cluster = Spot_manifest$mc, gaps_col = gaps_col)

    pdf(paste0(savePath,"SC_diffgene_inte_new.pdf"), width=15, height=10)
    #pheatmap(show_sc_matrix, annotation_col =annocol_cs,cluster_rows =T,cluster_cols =F,show_colnames = F,annotation_colors = ann_colors_sc,cutree_cols = length(cols))
    pheatmap(de.pre$expr.data,
             color = colorRampPalette(c("#4393C3", "white", "#D6604D"))(100),
             annotation_col = annocol_cs,
             annotation_colors = ann_colors_sc,
             fontsize = 7,
             gaps_col = de.pre$gaps_col,
             cluster_rows = F, cluster_cols = F,
             show_colnames = F)
    dev.off()

    #SC_diff_gene_show
    SC_barcode_to_MC <- data.frame(barcode = colnames(SC_obj),Mc_id = Mc_pair_SC)
    scrowid <- c()
    annocol_cs <- c()
    mc_id <- c()
    for(i in 1:length(mc_cluster)){
      tmidlist <- SC_barcode_to_MC$barcode[which(SC_barcode_to_MC$Mc_id==mc_cluster[i])]
      annocol_cs <- c(annocol_cs,rep(Mc_manifest$mc_type[Mc_manifest$mc_names==mc_cluster[i]],Mc_manifest$spot_num.Freq[Mc_manifest$mc_names==mc_cluster[i]]))
      scrowid <- c(scrowid,tmidlist)
      mc_id <- c(mc_id,rep(mc_cluster[i],Mc_manifest$spot_num.Freq[Mc_manifest$mc_names==mc_cluster[i]]))
    }
    annocol_cs <- data.frame(mc_id = mc_id)
    rownames(annocol_cs) <- scrowid
    ann_colors_sc = list(
      mc_id = cols
    )


    SC_obj <- FindVariableFeatures(SC_obj)
    scrowid <- ST_filter_str(scrowid,"-")
    SC_corr <- rcorr(sc_matrix[VariableFeatures(SC_obj),scrowid])
    sc_corr <- (SC_corr$r)
    tm <- sc_corr
    diag(tm) <- 0
    sc_max <- max(tm)
    diag(sc_corr) <- sc_max

    png(filename=paste0(savePath,"SC_corr.png"), width=1000, height=900)
    pheatmap(as.matrix(sc_corr), annotation_col = annocol_cs,annotation_row = annocol_cs,cluster_rows =F,cluster_cols =F,show_colnames = F,show_rownames = F,legend = F,annotation_colors = ann_colors_sc)
    dev.off()
  }

  return(MD_markers)
}



