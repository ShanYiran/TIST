#' Title
#'
#' @param exprPath The address where the RNA-seq file is stored. Such as "/home/data/filtered_feature_bc_matrix/".
#' @param spaceFile The address where the local file of all spots is stored. Such as "spatial/tissue_positions_list.csv".
#' @param savePath The address where the results stored.
#'
#' @return
#' @export
#'
#' @examples
#' \donttest{
#'    Spark_methods(exprPath,
#'    spaceFile,
#'    savePath)
#' }
#'
#'
Spark_methods <- function(exprPath,
                          spaceFile,
                          savePath){
  expr_obj <- Read10X(data.dir = exprPath)
  expr_obj <- CreateSeuratObject(counts = expr_obj,min.cells = 3,min.features = 200)
  expr_obj <- NormalizeData(object = expr_obj, normalization.method = "LogNormalize")
  expr_obj <- FindVariableFeatures(expr_obj,selection.method = "vst", nfeatures = 2000)
  Spot_space <- read.csv(file = spaceFile,header = F)
  sel.cols <- c("barcode","tissue", "row", "col", "imagerow", "imagecol")
  colnames(Spot_space) <- sel.cols
  Spot_space$barcode <- ST_filter_str(Spot_space$barcode,'-')
  colnames(expr_obj@assays$RNA@data) <- ST_filter_str(colnames(expr_obj@assays$RNA@data),'-')
  Spot_space <- Spot_space[which(Spot_space$barcode%in%colnames(expr_obj@assays$RNA@data)),]
  Spot_manifest <- Spot_space

  rawcount <- expr_obj@assays$RNA@counts[which(rownames(expr_obj@assays$RNA@counts)%in%VariableFeatures(expr_obj)),]
  count <- as.matrix(rawcount)
  colnames(count) <- ST_filter_str(colnames(count),'-')
  #Spot_manifest <- Spot_manifest[which(colnames(count)%in%Spot_manifest$barcode),]
  setorder(setDT(Spot_manifest)[, barcode := colnames(count)], barcode)

  info <- Spot_manifest[,c(3,4)]
  info <- as.data.frame(info)
  colnames(info) <- c("y","x")
  info$y <- 0-info$y
  rownames(info) <- Spot_manifest$barcode
  spark <- CreateSPARKObject(counts = (count), location = info[, 1:2],
                             percentage = 0.1, min_total_counts = 10)

  spark@lib_size <- apply(spark@counts, 2, sum)
  spark <- spark.vc(spark, covariates = NULL, lib_size = spark@lib_size, verbose = T, fit.maxiter = 500)
  spark <- spark.test(spark, check_positive = T, verbose = T)
  saveRDS(spark, file = paste0(savePath,"SPARK.rds"))
  #head(spark@res_mtest[,c("combined_pvalue","adjusted_pvalue")])
  #features3 <- spark@res_mtest
  #features3<- features3[sort(features3$combined_pvalue,index.return=TRUE)$ix,]
}


Spark_X_methods <- function(exprPath,
                            spaceFile,
                            savePath){
  expr_obj <- Read10X(data.dir = exprPath)
  expr_obj <- CreateSeuratObject(counts = expr_obj,min.cells = 3,min.features = 200)
  expr_obj <- NormalizeData(object = expr_obj, normalization.method = "LogNormalize")

  sp_count <- as.matrix(expr_obj@assays$RNA@counts)

  Spot_space <- read.csv(file = spaceFile,header = F)
  sel.cols <- c("barcode","tissue", "row", "col", "imagerow", "imagecol")
  colnames(Spot_space) <- sel.cols
  #Spot_space$barcode <- ST_filter_str(Spot_space$barcode,'-')
  Spot_space <- Spot_space[which(Spot_space$barcode%in%colnames(expr_obj@assays$RNA@data)),]
  Spot_manifest <- Spot_space
  info <- Spot_manifest[,c("row","col")]
  rownames(info) <- Spot_manifest$barcode
  location <- as.matrix(info)

  mt_idx <- grep("mt-",rownames(sp_count))
  if(length(mt_idx)!=0){
    sp_count    <- sp_count[-mt_idx,]
  }

  sparkX <- sparkx(sp_count,location,numCores=1,option="mixture")

  saveRDS(sparkX, file = paste0(savePath,"sparkX.rds"))
  #head(spark@res_mtest[,c("combined_pvalue","adjusted_pvalue")])
  #features3 <- spark@res_mtest
  #features3<- features3[sort(features3$combined_pvalue,index.return=TRUE)$ix,]
}
