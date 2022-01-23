#' Title
#'
#' @param SC_obj You can also load RNA-seq use an existed Seruat object.
#' @param Spot_manifest The TIST results returned by function "Meta_St_img_unsupervised"
#' @param savePath The address where the results stored.
#' @param centre The covariance matrix of distance Gaussian kernel, we default element in covariance matrix is identical.
#' @param broad The mean of distance Gaussian kernel.
#' @param spaceFile The address where the local file of all spots is stored. Such as "spatial/tissue_positions_list.csv".
#' @param Maskfile The address where the Mask file is stored. Such as "/home/data/Imginit/mask.txt". You can use no mask to generate an all 1 matrix. The more recommended approach is generate this mask file use our python code in this package.
#' @param imagefile The address where the Image file is stored. Such as "/home/data/spatial/tissue_hires_image.png".
#' @param spot_r_min A min window to select Image features. Default is 12.
#' @param img_import The image feature weight for feature integration.
#' @param spot_r_max A max window to select Image features. Default is 20.
#' @param Step Control the step length of random walk for walkstrap method.
#' @param diffusion_rate The diffusion rate, is related to temperature, liquid viscosity, and molecular size.
#' @param netfile The TIST-net build by function "Meta_St_img_unsupervised"
#' @param cluster_num The number of init cluster for image. Default is 50.
#' @param geneenrichment Whether to perform gene enhancement analysis. Default True.
#' @param genelist if geneenrichment is True, the analysis gene list.
#'
#' @examples
#' \donttest{
#'     Diffusion_net_analysis(SC_obj=SC_obj,
#'     Spot_manifest = Spot_manifest_imgunsup,
#'     savePath= savePath,
#'     centre = 1,
#'     broad = 2,
#'     spaceFile=spaceFile,
#'     Maskfile=Maskfile,
#'     imagefile=imagefile,
#'     spot_r_min = 12,
#'     img_import = 1,
#'     spot_r_max = 20,
#'     Step = NULL,
#'     diffusion_rate = dd,
#'     netfile = netfile,
#'     geneenrichment = T,
#'     genelist = NULL)
#' }
#'
Diffusion_net_analysis <- function(SC_obj,
                            Spot_manifest,
                            savePath,
                            centre = 1,
                            broad = 1,
                            spaceFile,
                            Maskfile,
                            imagefile,
                            spot_r_min = 12,
                            img_import = 1,
                            spot_r_max = 20,
                            Step = 10,
                            diffusion_rate = 0.1,
                            enrichment_score = TRUE,
                            netfile,
                            cluster_num = 50,
                            geneenrichment = F,
                            genelist = NULL
                            ){
  img <- load.image(imagefile)
  #image_info(img)
  Mask <- read.table(file = Maskfile,sep = ',')
  Mask <- as.matrix(Mask)
  Mask <- as.cimg(t(Mask))
  img_m <- img
  R(img_m) <-R(img_m)*Mask
  G(img_m) <-G(img_m)*Mask
  B(img_m) <-B(img_m)*Mask
  #grayscale(img_m) %>% hist(main="Luminance values in boats picture")
  #R(img_m) %>% hist(main="Red channel values in boats picture")
  #bdf <- as.data.frame(img_m)
  #head(bdf,3)
  #bdf <- mutate(bdf,channel=factor(cc,labels=c('R','G','B')))
  #ggplot(bdf,aes(value,col=channel))+geom_histogram(bins=30)+facet_wrap(~ channel)
  img_org <- img_m
  im <- grayscale(img_org) %>% isoblur(2)
  im <- im * Mask
  #im <- imsharpen(im,2)
  #im <- isoblur(im,10)
  #im <- grayscale(img_org)

  ###cal marcov RF for im
  mcim <- as.matrix(im)
  SC_obj <- ScaleData(SC_obj)
  #RenameCells(SC_obj,new.names = ST_filter_str(colnames(SC_obj),'-'))
  expr_matrix <- as.matrix(SC_obj@assays$RNA@counts[VariableFeatures(SC_obj),])
  colnames(expr_matrix) <- ST_filter_str(colnames(SC_obj),'-')
  diffusion_matrix <- matrix(0,nrow = dim(expr_matrix)[1],ncol = dim(expr_matrix)[2])
  space <- Spot_manifest
  rownames(space) <- Spot_manifest$barcode
  Count <- colMeans(expr_matrix)
  for(i in 1:length(colnames(expr_matrix))){
    this_node <- colnames(expr_matrix)[i]
    dis <- abs(space$row-space[this_node,"row"])+abs(space$col-space[this_node,"col"])/2
    dis_score <- centre * exp(-(dis*dis)/(broad*broad))
    Dc <-  Count-Count[i]
    sc <- Dc * dis_score
    sc[which(sc<0)] <- 0
    if(sum(sc)==0) next
    sc <- sc/sum(sc)
    tm_diff_matrix <- apply(expr_matrix, 1, function(x){x*sc*diffusion_rate})
    tm_diff_matrix <- t(tm_diff_matrix)
    tm_diff_matrix[,i] <- expr_matrix[,i]*(1-diffusion_rate)
    diffusion_matrix <- diffusion_matrix + tm_diff_matrix
  }
  diffusion_matrix <- round(diffusion_matrix)
  if(F){
    id <-diffusion_matrix["Spink8",]
    #id <- Dc
    id <- abs(id)
    pltdat <- cbind(Spot_manifest[, c("imagerow","imagecol")],id)
    colnames(pltdat)[1:2] <- c("y","x")
    pltdat$y = -pltdat$y
    xpand = 0
    ypand = 1
    pal <- colorRampPalette(c("#ADADAD", "#EEB4B4", "#7A378B"))
    p <- ggplot(pltdat, aes(x = x, y = y, color = id)) + geom_point(size = 3)+ scale_color_gradientn(colours = pal(5)) +scale_x_discrete(expand = c(xpand, ypand)) + scale_y_discrete(expand = c(xpand, ypand)) + coord_equal() +

      # labs(title = colnames(pd)[igene+2], x = NULL, y = NULL)+
      theme_bw()
    ggsave(filename = paste0(savePath, "/Diffusion_c1_b5/","Spink8_final.png"), p,
           width = 8, height = 7, dpi = 500)


    id <-SC_obj@assays$RNA@counts["Spink8",]
    #id <- Dc
    id <- abs(id)
    pltdat <- cbind(Spot_manifest[, c("imagerow","imagecol")],id)
    colnames(pltdat)[1:2] <- c("y","x")
    pltdat$y = -pltdat$y
    xpand = 0
    ypand = 1
    pal <- colorRampPalette(c("#ADADAD", "#EEB4B4", "#7A378B"))
    p <- ggplot(pltdat, aes(x = x, y = y, color = id)) + geom_point(size = 3)+ scale_color_gradientn(colours = pal(5)) +scale_x_discrete(expand = c(xpand, ypand)) + scale_y_discrete(expand = c(xpand, ypand)) + coord_equal() +

      # labs(title = colnames(pd)[igene+2], x = NULL, y = NULL)+
      theme_bw()
    ggsave(filename = paste0(savePath, "/Diffusion_c1_b5/","Spink8_before.png"), p,
           width = 8, height = 7, dpi = 500)

  }

  if(FALSE){
    id <-tm_diff_matrix["Spink8",]
    #id <- Dc
    #id[which(id==max(id))] = mean(id)
    id <- abs(id)
    pltdat <- cbind(Spot_manifest[, c("imagerow","imagecol")],id)
    colnames(pltdat)[1:2] <- c("y","x")
    pltdat$y = -pltdat$y
    xpand = 0
    ypand = 1
    pal <- colorRampPalette(c("#ADADAD", "#EEB4B4", "#7A378B"))
    p <- ggplot(pltdat, aes(x = x, y = y, color = id)) + geom_point(size = 3)+ scale_color_gradientn(colours = pal(5)) +scale_x_discrete(expand = c(xpand, ypand)) + scale_y_discrete(expand = c(xpand, ypand)) + coord_equal() +

      # labs(title = colnames(pd)[igene+2], x = NULL, y = NULL)+
      theme_bw()
    ggsave(filename = paste0(savePath, "/Diffusion_c1_b5/",i,"Spink8.png"), p,
           width = 8, height = 7, dpi = 500)


    id <-dis_score
    #id <- Dc
    id <- abs(id)
    pltdat <- cbind(Spot_manifest[, c("imagerow","imagecol")],id)
    colnames(pltdat)[1:2] <- c("y","x")
    pltdat$y = -pltdat$y
    xpand = 0
    ypand = 1
    pal <- colorRampPalette(c("#ADADAD", "#EEB4B4", "#7A378B"))
    p <- ggplot(pltdat, aes(x = x, y = y, color = id)) + geom_point(size = 3)+ scale_color_gradientn(colours = pal(5)) +scale_x_discrete(expand = c(xpand, ypand)) + scale_y_discrete(expand = c(xpand, ypand)) + coord_equal() +

      # labs(title = colnames(pd)[igene+2], x = NULL, y = NULL)+
      theme_bw()
    ggsave(filename = paste0(savePath, "/Diffusion_c1_b5/",i,"_before.png"), p,
           width = 8, height = 7, dpi = 500)

  }

  dir.create(paste0(savePath,"Diffusion","_d",diffusion_rate,"_b",broad))
  savePath_ <- paste0(savePath,"Diffusion","_d",diffusion_rate,"_b",broad,"/")
  saveRDS(diffusion_matrix,paste0(savePath_,"diffusion_matrix.RDS"))
  #Diff_obj <- SC_obj
  #colnames(diffusion_matrix) <- colnames(expr_matrix)
  colnames(diffusion_matrix) <- ST_filter_str(colnames(diffusion_matrix),'-')
  TPM <- as.matrix(diffusion_matrix)
  TPM <- as(TPM, "sparseMatrix")
  Diff_obj <- CreateSeuratObject(diffusion_matrix)
  label <- readRDS(paste0(savePath,"McRFlabel.RDS"))
  expr_obj_scale <- Diff_obj
  expr_obj_scale <- NormalizeData(expr_obj_scale)
  #expr_obj_scale <- FindVariableFeatures(expr_obj_scale)
  expr_obj_scale <- ScaleData(expr_obj_scale)
  expr_obj_scale <- RunPCA(expr_obj_scale,features = rownames(expr_obj_scale))
  expr_obj_scale <- FindNeighbors(expr_obj_scale)
  cellcorr <- as.matrix(expr_obj_scale@graphs$RNA_snn)
  colnames(cellcorr) <- ST_filter_str(colnames(cellcorr),'-')
  rownames(cellcorr) <- ST_filter_str(rownames(cellcorr),'-')
  #expr_data <- as.matrix(expr_obj_scale@assays$RNA@data)
  expr_data <- expr_obj_scale@assays$RNA@scale.data
  Spot_space <- read.csv(file = spaceFile,header = F)
  sel.cols <- c("barcode","tissue", "row", "col", "imagerow", "imagecol")
  colnames(Spot_space) <- sel.cols
  Spot_space$barcode <- ST_filter_str(Spot_space$barcode,'-')
  colnames(expr_obj_scale@assays$RNA@scale.data) <- ST_filter_str(colnames(expr_obj_scale@assays$RNA@scale.data),'-')
  Spot_space <- Spot_space[which(Spot_space$barcode%in%colnames(expr_obj_scale@assays$RNA@scale.data)),]
  Spot_manifest <- Spot_space

  #sample_bar <- Spot_manifest$barcode[1:100]

  #Spot_manifest <- Spot_manifest[1:100,]


  #Spot_manifest <- Spot_manifest[sort(Spot_manifest$col,index.return=TRUE)$ix,]
  #Spot_manifest <- Spot_manifest[sort(Spot_manifest$row,index.return=TRUE)$ix,]
  nodes <- Spot_manifest[,c("barcode")]
  nodes <- cbind(1:length(nodes),nodes)
  nodes <- data.frame(nodes)
  colnames(nodes) <- c("id","name")
  edges <- data.frame()
  #tm_matrix <- as.matrix(grayscale(img_org))

  #cellcorr <- cellcorr[which(rownames(cellcorr)%in%sample_bar),which(colnames(cellcorr)%in%sample_bar)]

  #tm_matrix <- as.matrix(im)
  tm_matrix <- as.matrix(label)
  tm_matrix_org <- as.matrix(im)

  col <- array(c("#f5f5f5",getDefaultColors(n = 400)))
  rownames(col) <- c(0,1:400)
  len <- length(colnames(cellcorr))

  sample_id <- which(cellcorr>=0.05)
  for(td in 1:length(sample_id)){
    tt <- sample_id[td]
    i = ceiling(tt/len)
    j = tt%%len
    if(j==0) j = len
    if(i<=j) next
    idx <- which(Spot_manifest$barcode==colnames(cellcorr)[i])
    idy <- which(Spot_manifest$barcode==colnames(cellcorr)[j])
    xi <- Spot_manifest$row[idx]
    yi <- Spot_manifest$col[idx]
    xj <- Spot_manifest$row[idy]
    yj <- Spot_manifest$col[idy]
    ximgrow = round(Spot_manifest[idx,"imagerow"]*sacle_score)
    ximgcol = round(Spot_manifest[idx,"imagecol"]*sacle_score)
    yimgrow = round(Spot_manifest[idy,"imagerow"]*sacle_score)
    yimgcol = round(Spot_manifest[idy,"imagecol"]*sacle_score)

    distodis = abs(xj-xi)+abs(yj-yi)#

    M1 = tm_matrix[(ximgcol-spot_r_max):(ximgcol+spot_r_max),(ximgrow-spot_r_max):(ximgrow+spot_r_max)]
    M2 = tm_matrix[(yimgcol-spot_r_max):(yimgcol+spot_r_max),(yimgrow-spot_r_max):(yimgrow+spot_r_max)]
    if(length(which(M1>0))==0) next
    if(length(which(M2>0))==0) next

    M1_v <- hist(as.vector(M1)[which(as.vector(M1)>0)],breaks = seq(0,cluster_num,length.out = cluster_num),plot = F)$counts/sum( hist(as.vector(M1)[which(as.vector(M1)>0)],breaks = seq(0,cluster_num,length.out = cluster_num),plot = F)$counts)
    M2_v <- hist(as.vector(M2)[which(as.vector(M2)>0)],breaks = seq(0,cluster_num,length.out = cluster_num),plot = F)$counts/sum( hist(as.vector(M2)[which(as.vector(M2)>0)],breaks = seq(0,cluster_num,length.out = cluster_num),plot = F)$counts)
    weight_1 <- 1-max(min(distance(rbind(M1_v,M2_v),method = 'kullback-leibler',test.na = F,mute.message = T),1))

    M1 = tm_matrix[(ximgcol-spot_r_min):(ximgcol+spot_r_min),(ximgrow-spot_r_min):(ximgrow+spot_r_min)]
    M2 = tm_matrix[(yimgcol-spot_r_min):(yimgcol+spot_r_min),(yimgrow-spot_r_min):(yimgrow+spot_r_min)]
    M1_v <- hist(as.vector(M1)[which(as.vector(M1)>0)],breaks = seq(0,cluster_num,length.out = cluster_num),plot = F)$counts/sum( hist(as.vector(M1)[which(as.vector(M1)>0)],breaks = seq(0,cluster_num,length.out = cluster_num),plot = F)$counts)
    M2_v <- hist(as.vector(M2)[which(as.vector(M2)>0)],breaks = seq(0,cluster_num,length.out = cluster_num),plot = F)$counts/sum( hist(as.vector(M2)[which(as.vector(M2)>0)],breaks = seq(0,cluster_num,length.out = cluster_num),plot = F)$counts)
    weight_2 <-  1-max(min(distance(rbind(M1_v,M2_v),method = 'kullback-leibler',test.na = F,mute.message = T),1))

    M1 = tm_matrix_org[(ximgcol-spot_r_max):(ximgcol+spot_r_max),(ximgrow-spot_r_max):(ximgrow+spot_r_max)]
    M2 = tm_matrix_org[(yimgcol-spot_r_max):(yimgcol+spot_r_max),(yimgrow-spot_r_max):(yimgrow+spot_r_max)]
    M1_v <- as.vector(M1)[which(as.vector(M1)>0)]
    M2_v <- as.vector(M2)[which(as.vector(M2)>0)]
    mean_diff_max <- abs(mean(M1_v)-mean(M2_v))*255
    M1 = tm_matrix_org[(ximgcol-spot_r_min):(ximgcol+spot_r_min),(ximgrow-spot_r_min):(ximgrow+spot_r_min)]
    M2 = tm_matrix_org[(yimgcol-spot_r_min):(yimgcol+spot_r_min),(yimgrow-spot_r_min):(yimgrow+spot_r_min)]
    M1_v <- as.vector(M1)[which(as.vector(M1)>0)]
    M2_v <- as.vector(M2)[which(as.vector(M2)>0)]
    mean_diff_min <- abs(mean(M1_v)-mean(M2_v))*255

    weight_3 <- min(1,1/(mean_diff_max + mean_diff_min + 1))

    weight <- min((cellcorr[i,j]+img_import*0.5*(weight_1+weight_2+weight_3))/(img_import+1),1)
    weight <- min(weight*10/distodis,1)  #
    #weight <- cellcorr[i,j] + weight*3/distodis
    edges <- rbind(edges,c(nodes$id[which(nodes$name==colnames(cellcorr)[i])],nodes$id[which(nodes$name==colnames(cellcorr)[j])],weight))
  }


  colnames(edges) <- c("from","to","weight")
  edges <- unique(edges)
  edges$weight <- as.numeric(edges$weight)
  ST_net <- graph_from_data_frame(
    d = edges,
    vertices = nodes,
    directed = F)
  saveRDS(ST_net,file = paste0(savePath_,"ST_imgunsupnet.RDS"))
  expr_data <- expr_obj_scale@assays$RNA@scale.data

  #hist(E(ST_net)$weight[which(E(ST_net)$weight>0.01)])

  expr_info_matrix <- t(expr_data)
  search_map <- list(c(2,0),c(-1,-1),c(1,-1),c(-2,0),c(-1,1),c(1,1))

  if(is.null(Step)){
    step_testlist <- 1:50
    step_testSDbw <- rep(0,50)
    discrete_points_num <- rep(0,50)
    I_Score <- rep(0,50)
    step_c_dis <- rep(0,50)
    c_num <- rep(0,50)
    for(i in 1:50){
      net_walktrap_test<-walktrap.community(ST_net,weights=abs(E(ST_net)$weight),step=i)
      walktrap_id <- cbind(net_walktrap_test$names,net_walktrap_test$membership)
      walktrap_id <- as.data.frame(walktrap_id)
      colnames(walktrap_id) <- c("barcode","Walktrap_id")
      Spot_manifest_t <- merge(Spot_manifest,walktrap_id,by.x = "barcode",by.y = "barcode")
      sort(table(walktrap_id$Walktrap_id))
      Spot_manifest_t$Walktrap_id <- as.numeric(Spot_manifest_t$Walktrap_id)
      use_set <- names(table(Spot_manifest_t$Walktrap_id))[which(table(Spot_manifest_t$Walktrap_id)>5)]
      Spot_manifest_t$Walktrap_id[which(!Spot_manifest_t$Walktrap_id%in%use_set)] <- 0
      Spot_manifest_t$Walktrap_id <- Renamedomin(namelist = Spot_manifest_t$Walktrap_id)
      c_num[i] <- length(table(Spot_manifest_t$Walktrap_id))

      #Number of discrete points
      for(k in 1:dim(Spot_manifest_t)[1]){
        xrow  = Spot_manifest_t[k,"row"]
        xcol  = Spot_manifest_t[k,"col"]
        ld <- Spot_manifest_t[k,"Walktrap_id"]
        tm = 0;
        for(j in 1:length(search_map)){
          yrow = xrow + search_map[[j]][1]
          ycol = xcol + search_map[[j]][2]
          id <- which(Spot_manifest_t$row==yrow & Spot_manifest_t$col==ycol)
          if(length(id)==0) tm = tm+1
          if(length(id)!=0){
            if(as.numeric(Spot_manifest_t[id,"Walktrap_id"])!=ld) tm=tm+1
          }
        }
        if(tm==6) discrete_points_num[i] = discrete_points_num[i] + 1
      }

      step_testSDbw[i] <- modularity(net_walktrap_test)
      step_c_dis[i] <- (c_num[i] /discrete_points_num[i] )
      I_Score[i] <- step_testSDbw[i] + step_c_dis[i]

      Spot_manifest_t$Walktrap_id <- as.character(Spot_manifest_t$Walktrap_id)
      ST_meta_plot_tm <- ggplot(Spot_manifest_t, aes(x = imagecol, y = -imagerow, color = Walktrap_id)) + geom_point(size = 2) +
        scale_colour_manual(values = col) + scale_x_discrete(expand = c(0.09,0.09)) + scale_y_discrete(expand = c(0.09,0.09)) + coord_equal() +  theme_bw() +
        xlab("x") + ylab("y")+ theme( legend.position = "none" )
      ggsave(filename = paste0(savePath_, "/ST_meta_imgunsup_plot_Step_",i,"_C",c_num[i],".png"), ST_meta_plot_tm,
             width = 8, height = 7, dpi = 500)
    }

    score_df <- cbind(c_num,step_testSDbw,I_Score,discrete_points_num)
    score_df <- as.data.frame(score_df)

    saveRDS(score_df,file = paste0(savePath_,"score_df.rds"))

    score_df <- score_df[sort(score_df$c_num,index.return=TRUE,decreasing = T)$ix,]
    attach(score_df)
    tm <- aggregate(score_df, by=list(c_num = score_df$c_num), FUN = max)

    png(file =paste0(savePath_,"Iscore.png"))
    plot(tm$c_num,tm$I_Score,type = "o",col = "red", xlab = "c_num", ylab = "I_Score",
         main = "I_Score")
    dev.off()

    png(file =paste0(savePath_,"step_clusternum_imgunsup.png"))
    plot(tm$c_num,tm$c_num,type = "o",col = "blue", xlab = "c_num", ylab = "c_num",
         main = "c_num")
    dev.off()
    png(file =paste0(savePath_,"SDbw_imgunsup.png"))
    plot(tm$c_num,tm$step_testSDbw,type = "o",col = "red", xlab = "c_num", ylab = "SDbw_img",
         main = "SDbw_img")
    dev.off()
    png(file =paste0(savePath_,"discrete_points_num.png"))
    plot(tm$c_num,tm$discrete_p,type = "o",col = "red", xlab = "c_num", ylab = "discrete_p",
         main = "discrete_points_num")
    dev.off()

    Step = which.max(I_Score)
    message("Step=",Step)
  }

  net_walktrap<-walktrap.community(ST_net,weights=E(ST_net)$weight,steps = Step,merges = T)
  #dendPlot(net_walktrap)
  walktrap_id <- cbind(net_walktrap$names,net_walktrap$membership)
  walktrap_id <- as.data.frame(walktrap_id)
  colnames(walktrap_id) <- c("barcode","Walktrap_id")
  Spot_manifest <- merge(Spot_manifest,walktrap_id,by.x = "barcode",by.y = "barcode")
  sort(table(walktrap_id$Walktrap_id))
  Spot_manifest$Walktrap_id <- as.numeric(Spot_manifest$Walktrap_id)
  use_set <- names(table(Spot_manifest$Walktrap_id))[which(table(Spot_manifest$Walktrap_id)>2)]
  Spot_manifest$Walktrap_id[which(!Spot_manifest$Walktrap_id%in%use_set)] <- 0
  use_set <- names(table(Spot_manifest$Walktrap_id))[which(table(Spot_manifest$Walktrap_id)>2)]
  Spot_manifest$Walktrap_id[which(!Spot_manifest$Walktrap_id%in%use_set)] <- 1
  Spot_manifest$Walktrap_id <- Renamedomin(namelist = Spot_manifest$Walktrap_id)
  Spot_manifest$Walktrap_id <- as.character(Spot_manifest$Walktrap_id)
  ST_meta_plot <- ggplot(Spot_manifest, aes(x = imagecol, y = -imagerow, color = Walktrap_id)) + geom_point(size = 2) +
    scale_colour_manual(values = col) + scale_x_discrete(expand = c(0.09,0.09)) + scale_y_discrete(expand = c(0.09,0.09)) + coord_equal() +  theme_bw() +
    xlab("x") + ylab("y")+ theme( legend.position = "none" )
  ggsave(filename = file.path(savePath_, "/ST_meta_imgunsup_plot.png"), ST_meta_plot,
         width = 8, height = 7, dpi = 500)


  Spot_manifest$Walktrap_id <- as.character(Spot_manifest$Walktrap_id)


  saveRDS(Spot_manifest,file = paste0(savePath_,"Spot_manifest_imgunsup.RDS"))
  write.csv(Spot_manifest,file = paste0(savePath_,"Spot_manifest_imgunsup.csv"))
  #score
  #tm <- I_score(netfile,membership(net_walktrap),savePath_,Spot_manifest,savefile = "SD_Iscore")
  tm <- SCmethod_score(label_file, Spot_manifest)

  #SNN
  Diff_obj_ <- Diff_obj
  Diff_obj_ <- FindVariableFeatures(Diff_obj_)
  Diff_obj_ <- ScaleData(Diff_obj_)
  Diff_obj_ <- RunPCA(Diff_obj_)
  Diff_obj_ <- FindNeighbors(Diff_obj_)
  snn_result <- FindClusters(Diff_obj_)
  idd <- as.double(Idents(snn_result))
  Spot_manifest[,7] <- idd
  id <- as.character(Idents(snn_result))
  pltdat <- cbind(Spot_manifest[, c("imagerow","imagecol")],id)
  colnames(pltdat)[1:2] <- c("y","x")
  pltdat$y = -pltdat$y
  xpand = 0
  ypand = 1
  pal <- colorRampPalette(c("#ADADAD", "#EEB4B4", "#7A378B"))
  cols_ <- array(getDefaultColors(length(table(id))))
  ST_meta_plot_tm <- ggplot(pltdat, aes(x = x, y = y, color = id)) + geom_point(size = 2) +
    scale_colour_manual(values = cols_) + scale_x_discrete(expand = c(0.09,0.09)) + scale_y_discrete(expand = c(0.09,0.09)) + coord_equal() +  theme_bw() +
    xlab("x") + ylab("y")+ theme( legend.position = "none" )
  ggsave(filename = paste0(savePath_, "SNN.png"), ST_meta_plot_tm,
         width = 8, height = 7, dpi = 500)

  names(idd) <- ST_filter_str(colnames(snn_result),'-')
  #tm <- rbind(tm,I_score(netfile,idd,savePath_,Spot_manifest,savefile = "SNN_Iscore"))

  tm <- rbind(tm,SCmethod_score(label_file, Spot_manifest))
  #gene_plot
  if(geneenrichment){
    if(is.null(genelist)){
      genelist = c('Spink8',
                   'Slc17a6',
                   'Kcne2',
                   'Prkcd',
                   'Crlf1',
                   'Slc6a11',
                   'Ctxn3',
                   'Baiap3',
                   'Cbln4',
                   'Mfge8',
                   'Gpx3',
                   'Fam163b',
                   'Hap1',
                   'Rora',
                   '1110008P14Rik',
                   'Meig1',
                   'Arpp21',
                   'Ndn',
                   'Cabp7',
                   'Adarb1',
                   'Ttr',
                   'Mbp',
                   'Rarres2',
                   'Ptgds',
                   'Lrrc10b',
                   'Dsp',
                   'Mef2c',
                   'Adcy1',
                   'Atp2b1',
                   'Neurod6',
                   '6330403K07Rik',
                   'Amotl1',
                   'Snap25',
                   'Dynlrb2',
                   'Lamp5',
                   'Tbr1',
                   'Ramp3',
                   'Resp18',
                   'C1ql2',
                   'Agt',
                   'Plp1',
                   'Folr1',
                   'Camk2n1',
                   'Cfh',
                   'Camk2a',
                   'Rasgrf2',
                   'Mgp',
                   'Ccdc153',
                   'Mobp',
                   'Ahi1',
                   'Cbln1',
                   'Nptxr',
                   'Arpc5',
                   'Tmem212',
                   'Col1a2',
                   'Hpcal1',
                   'Klk8',
                   'Ddn',
                   'Lypd1',
                   '3110035E14Rik',
                   'Hpca',
                   'Kl',
                   'Icam5',
                   'Neurl1a',
                   'Dclk1',
                   'Mal',
                   'Epop',
                   'Enpp2',
                   'Tcf7l2',
                   'Slc6a13',
                   'Cnp',
                   'Zcchc12',
                   'Sparc')
    }

    genelist <- intersect(rownames(Diff_obj),genelist)
    Plot_all_resnet(SC_obj,
                    netfile = paste0(savePath_,"ST_imgunsupnet.RDS"),
                    Diff_obj = expr_obj_scale,
                    Spot_manifest,
                    savePath_,
                    gene_name_list = genelist,
                    gene_plot = T)

  }

  return(Spot_manifest)
}


I_score <- function(netfile, idents, savePath, Spot_manifest,savefile = "I_score"){
  search_map <- list(c(2,0),c(-1,-1),c(1,-1),c(-2,0),c(-1,1),c(1,1))
  net <- readRDS(netfile)
  cluster_score <- modularity(net,membership = idents)
  tid <- data.frame(barcode = names(idents),id = as.vector(idents))

  #Number of discrete points
  Spot_manifest <- merge(Spot_manifest,tid,by.x = "barcode",by.y = "barcode")
  cc <- length(table(Spot_manifest$id))
  discrete_points_num <- 0
  for(k in 1:dim(Spot_manifest)[1]){
    xrow  = Spot_manifest[k,"row"]
    xcol  = Spot_manifest[k,"col"]
    ld <- Spot_manifest[k,"id"]
    tm = 0
    for(j in 1:length(search_map)){
      yrow = xrow + search_map[[j]][1]
      ycol = xcol + search_map[[j]][2]
      id <- which(Spot_manifest$row==yrow & Spot_manifest$col==ycol)
      if(length(id)==0) tm = tm+1
      if(length(id)!=0){
        if(as.numeric(Spot_manifest[id,"id"])!=ld) tm=tm+1
      }
    }
    if(tm==6) discrete_points_num = discrete_points_num + 1
  }
  I_score <- cluster_score + (cc/discrete_points_num )
  tm <- data.frame(cluster_num = cc, discrete_points_num = discrete_points_num, modulatrity_score = cluster_score,Iscore = I_score)
  write.csv(tm,paste0(savePath,savefile,".csv"))
  saveRDS(tm,paste0(savePath,savefile,".RDS"))
  return(tm)
}



gene_recover_diffusion <- function(SC_obj,netfile,Diff_obj,Spot_manifest,savePath,gene_name,gene_plot = F){
  ST_net <- readRDS(file = netfile)
  SC_obj <- NormalizeData(SC_obj)
  colnames(SC_obj@assays$RNA@data) <- ST_filter_str(colnames(SC_obj),'-')

  expr_data <- SC_obj@assays$RNA@data[gene_name,]
  names(expr_data) <- colnames(SC_obj@assays$RNA@data)
  #x_o <- as.vector(expr_data+0.1)
  x_o <- as.vector(expr_data[which(expr_data>0.1)])
  res_o <- fitdistr(x_o,"lognormal")
  observed_hist <- hist(x_o,breaks = 10,plot = F)
  xfit <-seq(min(x_o), max(x_o), by=(max(x_o)-min(x_o))/length(x_o))
  yfit_o <-dlnorm(xfit,res_o[[1]][1], res_o[[1]][2])
  yfit_o <- yfit_o*diff(observed_hist$mids[1:2])*length(xfit)
  #lines(xfit, yfit_o, col="blue", lwd=2)
  l_o <- yfit_o - 1*sqrt(vcov(res_o)[1,1]+0.1)
  u_o <- yfit_o + 1*sqrt(vcov(res_o)[1,1]+0.1)
  #calculate expected plot
  #cfg <- cluster_fast_greedy(ST_net)
  nodes <- V(ST_net)[Spot_manifest$barcode]
  #g <- induced_subgraph(ST_net, nodes)
  #png(filename=paste0(savePath,Mc_name,"circle_plot.png"), width=500, height=500)
  #plot(g,layout = layout_in_circle)
  #dev.off()
  g <- ST_net
  x_e <- Diff_obj@assays$RNA@data[gene_name,]
  names(x_e) <- ST_filter_str(colnames(Diff_obj@assays$RNA@data),'-')
  x_e <- x_e+0.1
  #x_e <- x_e[which(x_e>0.1)]
  drop_matrix <- x_e
  space <- Spot_manifest
  rownames(space) <- Spot_manifest$barcode
  x_o <- drop_matrix
  x_tm <- drop_matrix
  nodes <- V(ST_net)[Spot_manifest$barcode]
  space <- Spot_manifest
  rownames(space) <- Spot_manifest$barcode

  for(i in 1:dim(space)[1]){
    this_node <- names(x_o)[i]
    #this_n <- neighbors(ST_net,this_node)
    this_e <-  E(g)[inc(this_node)]
    gs <- subgraph.edges(g,this_e)
    this_enode <- gs %>% ends(E(gs))
    if(dim(this_enode)[1]>1){
      for(t in 1:length(this_e)){
        if(this_enode[t,1]==this_node) x_o[i] <- x_o[i] + this_e[t]$weight*x_tm[this_enode[t,2]]
        if(this_enode[t,2]==this_node) x_o[i] <- x_o[i] + this_e[t]$weight*x_tm[this_enode[t,1]]
      }
    }
    x_o[i] <- x_o[i]/sum(this_e$weight)
  }
  message(gene_name)
  saveRDS(x_o,paste0(savePath,gene_name,"_enhance.RDS"))

  x_e <- as.vector(x_o)
  res_e <- fitdistr(x_e,"lognormal")
  expected_hist <- hist(x_e,breaks = 10,plot = F)
  yfit_e <-dlnorm(xfit,res_e[[1]][1], res_e[[1]][2])
  yfit_e <- yfit_e*diff(expected_hist$mids[1:2])*length(xfit)
  #lines(xfit, yfit_e, col="blue", lwd=2)
  l_e <- yfit_e - 1*sqrt(vcov(res_e)[1,1]+0.1)
  u_e <- yfit_e + 1*sqrt(vcov(res_e)[1,1]+0.1)
  mashi_dis <- mashi(a = yfit_e,b = yfit_o)
  KL_dis <- kl.dist(yfit_e,yfit_o)$D
  dat_plot_spot<- data.frame(xfit, yfit_o, l_o, u_o,yfit_e,l_e,u_e)
  if(gene_plot==T){
    tm_e <- data.frame(barcode = names(x_o),x_e = x_e)
    pltdat <- merge(Spot_manifest,tm_e, by.x = "barcode", by.y = "barcode")
    pltdat <- pltdat[, c("imagerow","imagecol","x_e")]
    colnames(pltdat)[1:2] <- c("y","x")
    pltdat$y = -pltdat$y
    pal <- colorRampPalette(c("mediumseagreen", "lightyellow2", "deeppink"))
    xpand = 0
    ypand = 1
    gpt <- ggplot(pltdat, aes(x = x, y = y, color = x_e)) + geom_point(size = 3) +
      # scale_color_gradientn(colours=pal(5))+
      scale_color_gradientn(colours = pal(5)) + scale_x_discrete(expand = c(xpand, ypand)) + scale_y_discrete(expand = c(xpand, ypand)) + coord_equal() +

      # labs(title = colnames(pd)[igene+2], x = NULL, y = NULL)+
      theme_bw()+guides(fill=guide_legend(title=" "))
    ggsave(filename = paste0(savePath,gene_name,"_after.png"),
           gpt, width = 7, height = 10, dpi = 150, limitsize = FALSE)

    tm_e <- data.frame(barcode = names(x_o),d_matrix = drop_matrix)
    pltdat <- merge(Spot_manifest,tm_e, by.x = "barcode", by.y = "barcode")
    pltdat <- pltdat[, c("imagerow","imagecol","d_matrix")]
    colnames(pltdat)[1:2] <- c("y","x")
    pltdat$y = -pltdat$y
    pal <- colorRampPalette(c("mediumseagreen", "lightyellow2", "deeppink"))
    xpand = 0
    ypand = 1
    gpt <- ggplot(pltdat, aes(x = x, y = y, color = d_matrix)) + geom_point(size = 3) +
      # scale_color_gradientn(colours=pal(5))+
      scale_color_gradientn(colours = pal(5)) + scale_x_discrete(expand = c(xpand, ypand)) + scale_y_discrete(expand = c(xpand, ypand)) + coord_equal() +

      # labs(title = colnames(pd)[igene+2], x = NULL, y = NULL)+
      theme_bw()+guides(fill=guide_legend(title=" "))
    ggsave(filename = paste0(savePath,gene_name,"_before.png"),
           gpt, width = 7, height = 10, dpi = 150, limitsize = FALSE)

  }
  return(list(Mahal = round(mashi_dis, 2),KL = round(KL_dis, 2)))
}



Plot_all_resnet <- function(SC_obj,
                            netfile,
                            Diff_obj,
                            Spot_manifest,
                            savePath,
                            gene_name_list,
                            gene_plot = F){
  #gene_name_list <- intersect(gene_name_list,VariableFeatures(SC_obj))
  #KL_matrix <- array(0,dim = c(length(gene_name_list), 1))
  #Mahal_matirx <- array(0,dim = c(length(gene_name_list), 1))
  KL_matrix <- c()
  Mahal_matirx <- c()

  for(i in 1:length(gene_name_list)){
    gene_name <- gene_name_list[i]
    ans <- gene_recover_diffusion(SC_obj,netfile,Diff_obj,Spot_manifest,savePath,gene_name,gene_plot)
    KL_matrix <- c(KL_matrix,ans$KL)
    Mahal_matirx <- c(Mahal_matirx,ans$Mahal)
  }
  names(KL_matrix) <- gene_name_list
  names(Mahal_matirx) <- gene_name_list
  KL_matrix[which(KL_matrix>1)] = 1

  saveRDS(Mahal_matirx,paste0(savePath,"Mahal_matirx.RDS"))
  saveRDS(KL_matrix,paste0(savePath,"KL_matrix.RDS"))
  write.csv(Mahal_matirx,paste0(savePath,"Mahal_matirx.csv"))
  write.csv(KL_matrix,paste0(savePath,"KL_matrix.csv"))
}


mashi <-function(a,b){
  return (((a-b)%*% t(t(a-b))) / cov(a,b))
}


