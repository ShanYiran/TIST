#Meta_expr_matrix
#' Title
#'
#' @param exprPath The address where the RNA-seq file is stored. Such as "/home/data/filtered_feature_bc_matrix/".
#' @param Spot_manifest The TIST results returned by function "Meta_St_img_unsupervised"
#' @param savePath The address where the results stored.
#' @param imagefile The address where the Image file is stored. Such as "/home/data/spatial/tissue_hires_image.png".
#' @param merge_method mean or max for neighbor expression
#' @param method TIST network method is "Walktrap_id"
#'
#' @return Mc_manifest
#' @export
#'
Meta_expr_matrix <- function(exprPath,
                             Spot_manifest,
                             savePath,
                             imagefile,
                             merge_method = "max",
                             method = "Walktrap_id"){
  expr_obj <- Read10X(data.dir = exprPath)
  expr_obj <- CreateSeuratObject(counts = expr_obj,min.cells = 3,min.features = 200)
  expr_obj <- NormalizeData(object = expr_obj, normalization.method = "LogNormalize")
  expr_obj <- FindVariableFeatures(expr_obj,selection.method = "vst", nfeatures = 2000)
  expr_obj_scale <- ScaleData(expr_obj, features = VariableFeatures(expr_obj))
  expr_col_name <- ST_filter_str(colnames(expr_obj),'-')
  mc_id <- Spot_manifest
  mc_id <- mc_id[which(mc_id$barcode%in%expr_col_name),]
  barcode <- data.frame(barcode = expr_col_name)
  mc_id <- merge(barcode,mc_id,by.x = "barcode",by.y = "barcode")

  if(method=="Walktrap_id"){
    mc_id <- mc_id$Walktrap_id
  }
  mc_index <-table(mc_id)
  mc_names <- paste0("SD_",1:length(mc_index))
  Meta_matrix <- matrix(0,nrow =length(rownames(expr_obj)) ,ncol = length(mc_names))
  Mc_manifest <- data.frame(mc_names = mc_names,spot_num = mc_index)
  #mc_turn <- data.frame(mc_names = mc_names,mc_index =mc_index)
  rownames(Mc_manifest) <- Mc_manifest$spot_num.mc_id
  Idents(expr_obj) <- Mc_manifest[as.character(mc_id),"mc_names"]
  saveRDS(expr_obj,paste0(savePath,"expr_obj.RDS"))
  for( i in 1:length(mc_names)){
    id <- mc_names[i]
    this_expM <- expr_obj@assays$RNA@counts[,which(Idents(expr_obj)==id),drop=F]
    if(merge_method =="mean"){
      Meta_matrix[,i] <- round(rowMeans(this_expM))
    }
    if(merge_method =="max"){
      Meta_matrix[,i] <- max.col(this_expM)
    }
  }
  colnames(Meta_matrix) <- mc_names
  rownames(Meta_matrix) <- rownames(expr_obj)
  Mc_obj <- CreateSeuratObject(counts = Meta_matrix)
  #Mc_obj <- NormalizeData(object = Mc_obj, normalization.method = "LogNormalize")
  #feature_num <- length(colnames(Mc_obj))-2
  #Mc_obj <- FindVariableFeatures(Mc_obj,selection.method = "vst", nfeatures = feature_num)
  #Mc_obj <- ScaleData(Mc_obj, features = rownames(Mc_obj))
  #Mc_obj <- RunPCA(Mc_obj,features = VariableFeatures(Mc_obj))

  write.csv(Meta_matrix,file = paste0(savePath,"Meta_matrix.csv"))
  saveRDS(Mc_obj,paste0(savePath,"Mc_obj.RDS"))
  mc_index <- data.frame(org.id = rownames(mc_index),mc.id <- mc_names)
  saveRDS(mc_index,paste0(savePath,"mc_index.csv"))
  rownames(mc_index) <- mc_index$org.id

  saveRDS(Idents(expr_obj),file = paste0(savePath,"MC_idents.RDS"))
  saveRDS(Mc_manifest,file = paste0(savePath,"Mc_manifest.RDS"))

  return(Mc_manifest)

}


#MC and spot paired analyses
#heatmap for high var genes
Matching_analysis <- function(SC_obj_file,
                              Spot_manifest,
                              Mc_manifest,
                              MC_obj_file,
                              MC_ident_file,
                              savePath,
                              netfile,
                              res_rate = 0.5,
                              plot_res_net = F,
                              pair_ana = F,
                              marker_pair_ana = F){
  #diff_gene_heatmap
  #top5_mc <- all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
  SC_obj <- readRDS(SC_obj_file)
  MC_obj <- readRDS(MC_obj_file)

  all.markers <- FindAllMarkers(SC_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  saveRDS(all.markers,file = paste0(savePath,"Mc_markers.RDS"))

  SC_obj <- ScaleData(SC_obj, features = rownames(SC_obj))
  Mc_pair_SC <- readRDS(file = MC_ident_file)
  names(Mc_pair_SC) <- ST_filter_str(names(Mc_pair_SC),"-")
  #DoHeatmap(SC_obj, features = top5$gene) + NoLegend()
  Spot_manifest <- cbind(Spot_manifest,Mc_pair_SC)

  colnames(Spot_manifest)[8] <- "mc"
  if(marker_pair_ana == T){

    markers_MC <- all.markers
    top20_mc <- markers_MC %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
    top5_mc <- markers_MC %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
    top20_mc <- top20_mc$gene
    top5_mc <- top5_mc$gene

    par(mar=c(10,0.5,0.5,0.5))
    cols <- array(getDefaultColors(length(table(Spot_manifest$Walktrap_id))))
    rownames(cols) <- Mc_manifest$mc_names

    mc_matrix <- MC_obj@assays$RNA@scale.data
    gene_list <- top20_mc
    #MC_diff_gene_show
    #mc_matrix <- as.data.frame(MC_obj@assays$RNA@scale.data)
    gl <- gene_list
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
    sc_matrix <- as.matrix(SC_obj@assays$RNA@scale.data)

    show_sc_matrix <- sc_matrix[top5_mc,rownames(annocol_cs)]
    show_sc_matrix[which(show_sc_matrix>4)] <- 4

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
  #cluster_corr_show
  #image(mc_corr)
  #image(mc_corr, xaxt='n', yaxt='n', xlab="MC", ylab="MC")
  #mtext("#cells", 1, line=0.5, cex = mcp_heatmap_text_cex, las=2)
  #ct.list <- list(
  #  "Hepatocytes" = c("HAMP", "CYP1A2", "APOF", "CYP3A4", "C9", "ADH4",
  #                    "SLC22A1", "CYP2A6", "HGFAC", "CYP2C8", "CYP2E1", "RDH16",
  #                    "SPP2", "UROC1", "AFM", "PCK1", "F9", "CYP4A11",
  #                    "SERPINA11", "APOA5", "CYP8B1", "SLC10A1"),
  #  "Malignancy" = c("SPINK1", "GPC3", "AKR1B10", "TOP2A", "REG3A", "CAP2",
  #                   "UBE2C", "MDK", "LCN2", "AURKA"),
  #  "Bile duct cells" = c("KRT19", "KRT7"),
  #  "Endothelial cells" = c("CLDN5", "PECAM1", "CD34", "FLT1", "VWF", "ENG", "CDH5"),
  #  "Fibroblasts" = c("COL1A2", "FAP", "PDPN", "DCN", "COL3A1", "COL6A1", "COL1A1"),
  #  "T cells" = c("CD3D", "CD3E", "CD3G"),
  #  "B cells" = c("CD19", "MS4A1", "CD79A"),
  #  "NK cells" = c("NCAM1", "KLRF1", "NCR1", "KLRC1"),
  #  "Myeloid cells" = c("ITGAX", "CD33", "CEACAM8", "CD68", "CD163"),
  #  "TLS" = c("FDCSP","CR2","CXCL13","LTF","CD52","MS4A1","CCL19","LINC00926","LTB","CORO1A","CD79B","TXNIP","CD19","LIMD2","CD37","ARHGAP45","BLK","TMC8","CCL21","PTPN6","ATP2A3","IGHM","SPIB","TMSB4X","CXCR4","NCF1","CD79A","ARHGAP9","DEF6","EVL","TBC1D10C","RASAL3","INPP5D","RNASET2","RASGRP2","TNFRSF13C","RAC2","CD22","ARHGEF1","AC103591.3","TRAF3IP3","HLA-DQB1","CD53","ARHGAP4","TRBC2","POU2AF1","TRAF5","OGA","FCRL3","HLA-DQA1"),
  #  "chemokine" = c("CCL2","CCL3","CCL4","CCL5","CCL8","CCL18","CCL19","CCL21","CXCL9","CXCL10","CXCL13","CXCL11"))
  #rownames(MC_obj@assays$RNA@scale.data) <- toupper(rownames(MC_obj@assays$RNA@scale.data))
  #gl <- intersect(rownames(MC_obj@assays$RNA@scale.data),unlist(ct.list))
  if(pair_ana == T){
    mc_genepair_show(gl1 = gl,gl2 = gl,MC_obj_file,Mc_manifest,savePath)
    sc_genepair_show(gl1 = gl,gl2 = gl,SC_obj_file,Spot_manifest,savePath)
  }
  gl <- all.markers %>% group_by(cluster) %>% top_n(n=3,wt=avg_logFC)
  gl <- gl$gene
  gl <- unique(gl)
  feature_plot(gene_plot = gl,SC_obj, Spot_manifest,savePath)
  mc_marker_plot(gl,MC_obj_file,Mc_manifest,Spot_manifest,savePath)
}

mc_marker_plot <- function(gl,MC_obj_file,Mc_manifest,Spot_manifest,savePath,savefilename=NULL,TLS =F){
  MC_obj <- readRDS(MC_obj_file)
  mc_expr <- MC_obj@assays$RNA@data
  #rownames(MC_obj@assays$RNA@scale.data) <- toupper(rownames(MC_obj@assays$RNA@scale.data))
  #rownames(mc_expr) <- toupper(rownames(mc_expr))
  if(is.null(gl)){
    ct.list <- list(
      "Hepatocytes" = c("HAMP", "CYP1A2", "APOF", "CYP3A4", "C9", "ADH4",
                        "SLC22A1", "CYP2A6", "HGFAC", "CYP2C8", "CYP2E1", "RDH16",
                        "SPP2", "UROC1", "AFM", "PCK1", "F9", "CYP4A11",
                        "SERPINA11", "APOA5", "CYP8B1", "SLC10A1"),
      "Malignancy" = c("SPINK1", "GPC3", "AKR1B10", "TOP2A", "REG3A", "CAP2",
                       "UBE2C", "MDK", "LCN2", "AURKA"),
      "Bile duct cells" = c("KRT19", "KRT7"),
      "Endothelial cells" = c("CLDN5", "PECAM1", "CD34", "FLT1", "VWF", "ENG", "CDH5"),
      "Fibroblasts" = c("COL1A2", "FAP", "PDPN", "DCN", "COL3A1", "COL6A1", "COL1A1"),
      "T cells" = c("CD3D", "CD3E", "CD3G"),
      "B cells" = c("CD19", "MS4A1", "CD79A"),
      "NK cells" = c("NCAM1", "KLRF1", "NCR1", "KLRC1"),
      "Myeloid cells" = c("ITGAX", "CD33", "CEACAM8", "CD68", "CD163"),
      "TLS" = c("FDCSP","CR2","CXCL13","LTF","CD52","MS4A1","CCL19","LINC00926","LTB","CORO1A","CD79B","TXNIP","CD19","LIMD2","CD37","ARHGAP45","BLK","TMC8","CCL21","PTPN6","ATP2A3","IGHM","SPIB","TMSB4X","CXCR4","NCF1","CD79A","ARHGAP9","DEF6","EVL","TBC1D10C","RASAL3","INPP5D","RNASET2","RASGRP2","TNFRSF13C","RAC2","CD22","ARHGEF1","AC103591.3","TRAF3IP3","HLA-DQB1","CD53","ARHGAP4","TRBC2","POU2AF1","TRAF5","OGA","FCRL3","HLA-DQA1"),
      "chemokine" = c("CCL2","CCL3","CCL4","CCL5","CCL8","CCL18","CCL19","CCL21","CXCL9","CXCL10","CXCL13","CXCL11"))
    gl <- intersect(rownames(MC_obj@assays$RNA@scale.data),unlist(ct.list))
  }
  gl <- intersect(rownames(MC_obj@assays$RNA@data),gl)
  expr_gl_mc <- matrix(0,nrow = length(gl),ncol = dim(Spot_manifest)[1])
  Spot_manifest$Walktrap_id <- as.numeric(Spot_manifest$Walktrap_id)
  for(i in 1:dim(expr_gl_mc)[1]){
    for(j in 1:dim(expr_gl_mc)[2]){
      expr_gl_mc[i,j] <- mc_expr[gl[i],Mc_manifest$mc_names[which(Mc_manifest$spot_num.mc_id==Spot_manifest$Walktrap_id[j])]]
    }
  }
  colnames(expr_gl_mc) <- Spot_manifest$barcode
  rownames(expr_gl_mc) <- gl
  pltdat <- cbind.data.frame(Spot_manifest[,c("imagerow","imagecol")], t(expr_gl_mc))
  genetitle <- gl
  colnames(pltdat)[1:2] <- c("y","x")
  pltdat$y = -pltdat$y
  if(is.null(savefilename)){
    png(filename=paste0(savePath,"MC_feature_plot.png"), width=1600, height=500*(length(gl)/4))
    pp <- lapply(1:(ncol(pltdat)-2), function(x) {
      pattern_plot2(pltdat, x, main = T, titlesize = 1.5, title = genetitle[x])
    })
    grid.arrange(grobs = pp, ncol = 4)
  }

  if(!is.null(savefilename)){
    png(filename=paste0(savePath,savefilename,".png"), width=1600, height=500*(length(gl)/4))
    pp <- lapply(1:(ncol(pltdat)-2), function(x) {
      pattern_plot2(pltdat, x, main = T, titlesize = 1.5, title = genetitle[x])
    })
    grid.arrange(grobs = pp, ncol = 4)
  }
  dev.off()
  #TLS
  if(TLS==T){
    gl <- c("CCL2","CCL3","CCL4","CCL5","CCL8","CCL18","CCL19","CCL21","CXCL9","CXCL10","CXCL13","CXCL11")
    #gl <- c("FDCSP","CR2","CXCL13","LTF","CD52","MS4A1","CCL19","LINC00926","LTB","CORO1A","CD79B","TXNIP","CD19","LIMD2","CD37","ARHGAP45","BLK","TMC8","CCL21","PTPN6","ATP2A3","IGHM","SPIB","TMSB4X","CXCR4","NCF1","CD79A","ARHGAP9","DEF6","EVL","TBC1D10C","RASAL3","INPP5D","RNASET2","RASGRP2","TNFRSF13C","RAC2","CD22","ARHGEF1","AC103591.3","TRAF3IP3","HLA-DQB1","CD53","ARHGAP4","TRBC2","POU2AF1","TRAF5","OGA","FCRL3","HLA-DQA1")
    gl <- intersect(rownames(expr_gl_mc),gl)

    tm_matrix <- expr_gl_mc[gl,]

    TLS_expr <- colSums(tm_matrix)

    TLS_expr_tm <- TLS_expr*0
    TLS_expr_tm[TLS_expr>=quantile(TLS_expr,0.984)] <- 1

    pltdat <-  cbind.data.frame(Spot_manifest[,c("imagerow","imagecol")], TLS_expr_tm)
    pal <- colorRampPalette(c("mediumseagreen", "lightyellow2", "deeppink"))
    colnames(pltdat)[1:2] <- c("y","x")
    pltdat$y = -pltdat$y
    png(filename=paste0(savePath,"MC_TLS.png"), width=600, height=500)
    ggplot(pltdat, aes(x = x, y = y, color = pltdat[, 1 + 2])) + geom_point(size = 2) +
      # scale_color_gradientn(colours=pal(5))+
      scale_color_gradientn(colours = pal(5)) + coord_equal() +
      theme(panel.grid = element_blank(),
            legend.position = "top",
            panel.border = element_blank(),
            plot.title = element_text(hjust = 0.5)) +
      ggplot_config(base.size = 6)
    dev.off()

    if(FALSE){
      Spot_manifest_tm <- Spot_manifest
      Spot_manifest$mc[which(Spot_manifest$mc%nin%c("SD_7"))] <- "0"
      col <- array(getDefaultColors(length(table(Spot_manifest$mc))))
      cols <- as.array(c("#DBDBDB",col[7],col[15]))
      names(cols) <- c("0","SD_7","15")
      ST_meta_plot2 <- ggplot(Spot_manifest, aes(x = imagecol, y = -imagerow, color = mc)) + geom_point(size = 2) +
        scale_colour_manual(values = cols) + scale_x_discrete(expand = c(0.09,0.09)) + scale_y_discrete(expand = c(0.09,0.09)) + coord_equal() +  theme_bw() +
        xlab("x") + ylab("y")
      ggsave(filename = file.path(savePath, "/TLS_final.png"), ST_meta_plot2,
             width = 8, height = 7, dpi = 500)
    }

    id <-TLS_expr
    pltdat <- cbind(Spot_manifest[, c("imagerow","imagecol")],id)
    colnames(pltdat)[1:2] <- c("y","x")
    pltdat$y = -pltdat$y
    xpand = 0
    ypand = 1
    pal <- colorRampPalette(c("#ADADAD", "#EEB4B4", "#7A378B"))
    p <- ggplot(pltdat, aes(x = x, y = y, color = id)) + geom_point(size = 3)+ scale_color_gradientn(colours = pal(5)) +scale_x_discrete(expand = c(xpand, ypand)) + scale_y_discrete(expand = c(xpand, ypand)) + coord_equal() +

      # labs(title = colnames(pd)[igene+2], x = NULL, y = NULL)+
      theme_bw()
    ggsave(filename = paste0(savePath, "TLS_score.png"), p,
           width = 8, height = 7, dpi = 500)



  }
}

mc_genepair_show <- function(gl1,
                             gl2,
                             MC_obj_file,
                             Mc_manifest,
                             savePath){
  MC_obj <- readRDS(MC_obj_file)
  mc_expr <- MC_obj@assays$RNA@scale.data
  for(i in 1:length(gl1)){
    for(j in 1:length(gl2)){
      gene1 <- gl1[i]
      gene2 <- gl2[j]
      if(gene1!=gene2){
        gene1_expr <- (mc_expr[gene1,]+1)
        gene2_expr <- (mc_expr[gene2,]+1)
        mc_id <- Mc_manifest$mc_type
        mc_gp_df <- data.frame(gene1_expr = gene1_expr,gene2_expr = gene2_expr,mc_id = mc_id)
        cols <- array(getDefaultColors(10))
        rownames(cols) <- c("Hepatocytes","Malignancy","Bile duct cells","Endothelial cells","Fibroblasts","T cells","B cells","NK cells","Myeloid cells","UnKnown")
        png(filename=paste0(savePath,gene1,"_",gene2,".png"), width=500, height=500)
        plot(mc_gp_df$gene1_expr,mc_gp_df$gene2_expr, cex= 3*1, col="black", pch=21, bg=cols[mc_gp_df$mc_id],xlab = gene1,ylab = gene2)
        text(mc_gp_df$gene1_expr,mc_gp_df$gene2_expr,1:length(mc_gp_df$mc_id), cex=0.9)
        dev.off()
      }
    }
  }
}


sc_genepair_show <- function(gl1,
                             gl2,
                             SC_obj_file,
                             Spot_manifest,
                             savePath){
  SC_obj <- readRDS(SC_obj_file)
  sc_expr <- SC_obj@assays$RNA@data
  for(i in 1:length(gl1)){
    for(j in 1:length(gl2)){
      gene1 <- gl1[i]
      gene2 <- gl2[j]
      if(gene1!=gene2){
        gene1_expr <- (sc_expr[gene1,])
        gene2_expr <- (sc_expr[gene2,])
        sc_id <- Spot_manifest$cell_type_new_2
        sc_gp_df <- data.frame(gene1_expr = gene1_expr,gene2_expr = gene2_expr,sc_id = sc_id)
        cols <- array(getDefaultColors(10))
        rownames(cols) <- c("Hepatocytes","Malignancy","Bile duct cells","Endothelial cells","Fibroblasts","T cells","B cells","NK cells","Myeloid cells","UnKnown")
        png(filename=paste0(savePath,gene1,"_",gene2,"_sc.png"), width=900, height=900)
        plot(sc_gp_df$gene1_expr,sc_gp_df$gene2_expr, cex= 1, col="black", pch=21, bg=cols[sc_gp_df$sc_id],xlab = gene1,ylab = gene2)
        #text(sc_gp_df$gene1_expr,sc_gp_df$gene2_expr, 1:length(sc_gp_df$gene2_expr), cex=0.9)
        dev.off()
      }
    }
  }
}


feature_plot <- function(gene_plot, SC_obj, Spot_manifest,savefile,savefilename=NULL,width = 1600,height = 500,pointsize = 2){
  #rawcount <- SC_obj@assays$RNA@counts[which(rownames( SC_obj@assays$RNA@counts)%in%in(gene_plot,VariableFeatures(SC_obj))),]
  #rownames(SC_obj@assays$RNA@data) <- toupper(rownames(SC_obj@assays$RNA@data))
  gene_plot <- intersect(gene_plot,rownames(SC_obj@assays$RNA@data))
  rel_vst_ct <- SC_obj@assays$RNA@data[gene_plot,]
  barcode <- ST_filter_str(colnames(rel_vst_ct),'-')
  rel_vst_ct <- cbind(barcode,as.data.frame(t(rel_vst_ct)))
  Spot_manifest <- merge(Spot_manifest,rel_vst_ct,x.by = "barcode",y.by = "barcode")
  pltdat <- cbind.data.frame(Spot_manifest[, c("imagerow","imagecol",gene_plot)])
  colnames(pltdat)[1:2] <- c("y","x")
  pltdat$y = -pltdat$y
  if(is.null(savefilename)){
    png(filename=paste0(savefile,"feature_plot.png"),  width=width, height=height*(length(gene_plot)/4))
    pp <- lapply(1:(ncol(pltdat) - 2), function(x) {
      pattern_plot2(pltdat, x ,main = T, titlesize = 1.5, title = gene_plot[x],pointsize = 2)
    })
    grid.arrange(grobs = pp, ncol = 4)
  }
  if(!is.null(savefilename)){
    png(filename=paste0(savefile,savefilename,".png"),  width=width, height=height*(length(gene_plot)/4))
    pp <- lapply(1:(ncol(pltdat) - 2), function(x) {
      pattern_plot2(pltdat, x ,main = T, titlesize = 1.5, title = gene_plot[x],pointsize = 2)
    })
    grid.arrange(grobs = pp, ncol = 4)
  }
  #pal <- colorRampPalette(c("mediumseagreen", "lightyellow2", "deeppink"))
  #gpt <- ggplot(pltdat, aes(x = x, y = y, color = pltdat[, 1 + 2])) + geom_point(size = 1) +
    # scale_color_gradientn(colours=pal(5))+
  #  scale_color_gradientn(colours = pal(5)) + coord_equal() +
    # labs(title = colnames(pd)[igene+2], x = NULL, y = NULL)+
  #  theme_bw()
  #find space spatially features _Method 3 SpatialDE??use python code??
  #gpt
  dev.off()
}





