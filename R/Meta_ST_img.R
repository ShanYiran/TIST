#' Title
#'
#' @param Maskfile The address where the Mask file is stored. Such as "/home/data/Imginit/mask.txt". You can use no mask to generate an all 1 matrix. The more recommended approach is generate this mask file use our python code in this package.
#' @param imagefile The address where the Image file is stored. Such as "/home/data/spatial/tissue_hires_image.png".
#' @param exprPath The address where the RNA-seq file is stored. Such as "/home/data/filtered_feature_bc_matrix/".
#' @param SC_obj If exprPath is NULL, you can also load RNA-seq use an existed Seruat object. Default is NULL.
#' @param spaceFile The address where the local file of all spots is stored. Such as "spatial/tissue_positions_list.csv".
#' @param colors Color bar for visualization. We will give a default color scheme if NULL.
#' @param savePath The address where the results stored.
#' @param Method The Network-based methods to find clusters. Default param is walktrap,
#' @param imgMethod The image feather select method. Default is Marcov, Other optional: Kmeans. We will use original image feature if NULL.
#' @param imginitMethod If imgMethod is Marcov, this param is to select label init method, such as "Kmeans" default or "Random".
#' @param Step Control the step length of random walk for walkstrap method.
#' @param spot_r_min A min window to select Image features. Default is 12.
#' @param spot_r_max A max window to select Image features. Default is 20.
#' @param cluster_num The number of init cluster for image. Default is 50.
#' @param maxiter The number of iterations for marcov. Default is 10.
#' @param sacle_score The scale to map the spots coordinate in spaceFile with image.
#' @param img_import The image feature weight for feature integration.
#'
#' @return TIST results.
#'
#' @examples
#' \donttest{
#'   Spot_manifest_imgunsup <- Meta_St_img_unsupervised(Maskfile = Maskfile,
#'   imagefile = imagefile,
#'   spaceFile = spaceFile,
#'   exprPath = exprPath,
#'   colors = NULL,
#'   savePath = savePath,
#'   Method = "walktrap",
#'   sacle_score = sacle_score)
#' }
#'
#'
Meta_St_img_unsupervised <- function(Maskfile,
                                     imagefile,
                                     exprPath,
                                     SC_obj = NULL,
                                     spaceFile,
                                     colors = NULL,
                                     savePath,
                                     Method = "walktrap",
                                     imgMethod = "Marcov",
                                     imginitMethod = "Kmeans",
                                     Step = NULL,
                                     spot_r_min = 12,
                                     spot_r_max = 20,
                                     cluster_num = 50,
                                     maxiter = 10,
                                     sacle_score = 0.09920635,
                                     img_import = 1){
  #image part
  img <- load.image(imagefile)
  #image_info(img)
  Mask <- read.table(file = Maskfile,sep = ',')
  Mask <- as.matrix(Mask)
  Mask <- as.cimg(Mask)
  img_m <- img
  R(img_m) <-R(img_m)*Mask
  G(img_m) <-G(img_m)*Mask
  B(img_m) <-B(img_m)*Mask
  plot(img_m)
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
  plot(im)


  ###cal marcov RF for im
  mcim <- as.matrix(im)


  if(imginitMethod == "Kmeans"){
    label = kmeans(as.vector(mcim),centers = cluster_num,nstart = 50)# time 1min
    label = array(label$cluster,dim = c(size(mcim)[1],size(mcim)[2]))
    label = label * as.matrix(Mask)
    plot(as.cimg(label*10))
    saveRDS(label,file = paste0(savePath,"Kmeans.RDS"))
  }
  if(imginitMethod == "Random"){
    label = randi(imax = cluster_num,n = size(mcim)[1],m = size(mcim)[2])
    label = label * as.matrix(Mask)
    plot(as.cimg(label*10))
  }


  if(imgMethod=="Marcov"){
    for( iter in 1:maxiter){
      label_img <- as.cimg(label)
      label_u <- as.matrix(convolve( label_img, filter = as.cimg(matrix(c(0,1,0,0,0,0,0,0,0),ncol = 3,byrow = T))))
      label_d <- as.matrix(convolve( label_img, filter = as.cimg(matrix(c(0,0,0,0,0,0,0,1,0),ncol = 3, byrow = T))))
      label_l <- as.matrix(convolve( label_img, filter = as.cimg(matrix(c(0,0,0,1,0,0,0,0,0),ncol = 3, byrow = T))))
      label_r <- as.matrix(convolve( label_img, filter = as.cimg(matrix(c(0,0,0,0,0,1,0,0,0),ncol = 3, byrow = T))))
      label_ul <- as.matrix(convolve( label_img, filter = as.cimg(matrix(c(1,0,0,0,0,0,0,0,0),ncol = 3, byrow = T))))
      label_ur <- as.matrix(convolve( label_img, filter = as.cimg(matrix(c(0,0,1,0,0,0,0,0,0),ncol = 3, byrow = T))))
      label_dl <- as.matrix(convolve( label_img, filter = as.cimg(matrix(c(0,0,0,0,0,0,1,0,0),ncol = 3, byrow = T))))
      label_dr <- as.matrix(convolve( label_img, filter = as.cimg(matrix(c(0,0,0,0,0,0,0,0,1),ncol = 3, byrow = T))))

      p_c <- array(0,dim = c(cluster_num, dim(label)[1]*dim(label)[2]))
      for(i in 1:cluster_num){
        label_i <- i * matrix(1,nrow = dim(label)[1],ncol = dim(label)[2])
        temp <- as.numeric(!(label_i - label_u)) + as.numeric(!(label_i - label_d)) +
          as.numeric(!(label_i - label_l)) + as.numeric(!(label_i - label_r)) +
          as.numeric(!(label_i - label_ul)) + as.numeric(!(label_i - label_ur)) +
          as.numeric(!(label_i - label_dl)) + as.numeric(!(label_i - label_dr))
        p_c[i,] <-as.vector(temp)/8
      }
      p_c[which(p_c==0)] <- 0.001;
      mu <- array(0,dim=c(1,cluster_num))
      sigma <- array(0,dim = c(1,cluster_num))
      for(i in 1:cluster_num){
        data_c <- mcim[which(label==i)]
        mu[i] = mean(data_c)
        sigma[i] = var(data_c)
      }
      tmmcim <- as.vector(mcim)
      p_sc <- array(0,dim=c(cluster_num,dim(label)[1]*dim(label)[2]))
      for(j in 1:cluster_num){
        MU = array(mu[j],dim = c(dim(mcim)[1]*dim(mcim)[2],1))
        p_sc[j,] = 1/sqrt(2*pi*sigma[j])*exp(-(tmmcim-MU)**2/2/sigma[j])
      }
      p_sc[which(p_sc==0)] <- 0.001;
      label_n <- array(0,dim = dim(label)[1]*dim(label)[2])
      tmlist <- log(p_c) + log(p_sc)
      tmlist[which(is.na(tmlist))] <- 0
      label_n <- apply(tmlist,2,which.max)
      #label_n <- apply(tmlist,2,function(t) which.max(t))
      # for(i in 1:(dim(label)[1]*dim(label)[2])){
      #    tmlist <- log(p_c[,i])+log(p_sc[,i])
      #    tmlist[which(is.na(tmlist))] <- 0
      #    label_n[i] = which(tmlist==max(tmlist))
      #  }
      label_new = array(label_n,dim = c(dim(label)[1],dim(label)[2]))
      label = label_new
      label = label * as.matrix(Mask)
    }

    plot(as.cimg(label))
    saveRDS(label,file = paste0(savePath,"McRFlabel.RDS"))
  }
  if(!is.null(SC_obj)){
    expr_obj_scale <- SC_obj
  }
  else{
    expr_obj <- Read10X(data.dir = exprPath)
    expr_obj <- CreateSeuratObject(counts = expr_obj,min.cells = 3,min.features = 200)
    expr_obj <- NormalizeData(object = expr_obj)
    expr_obj <- FindVariableFeatures(expr_obj,selection.method = "vst", nfeatures = 2000)
    expr_obj_scale <- ScaleData(expr_obj, features = VariableFeatures(expr_obj))
  }
  expr_obj_scale <- RunPCA(expr_obj_scale)
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

  if(!is.null(colors)){
    col <- colors
    rownames(col) <- c(1:length(colors))
  }
  else{
    col <- array(c("#f5f5f5",getDefaultColors(n = 400)))
    rownames(col) <- c(0,1:400)
  }

  len <- length(colnames(cellcorr))

  sample_id <- which(cellcorr!=0)
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

    M1_v <- hist(as.vector(M1)[which(as.vector(M1)>0)],breaks = seq(0,cluster_num,length.out = cluster_num),plot = F)$counts/(sum( hist(as.vector(M1)[which(as.vector(M1)>0)],breaks = seq(0,cluster_num,length.out = cluster_num),plot = F)$counts))
    M2_v <- hist(as.vector(M2)[which(as.vector(M2)>0)],breaks = seq(0,cluster_num,length.out = cluster_num),plot = F)$counts/(sum( hist(as.vector(M2)[which(as.vector(M2)>0)],breaks = seq(0,cluster_num,length.out = cluster_num),plot = F)$counts))
    weight_1 <- 1-max(min(distance(rbind(M1_v,M2_v),method = 'kullback-leibler',test.na = F,mute.message = T),1))

    M1 = tm_matrix[(ximgcol-spot_r_min):(ximgcol+spot_r_min),(ximgrow-spot_r_min):(ximgrow+spot_r_min)]
    M2 = tm_matrix[(yimgcol-spot_r_min):(yimgcol+spot_r_min),(yimgrow-spot_r_min):(yimgrow+spot_r_min)]
    if(min(M1)==max(M1)) next
    if(min(M2)==max(M2)) next
    M1_v <- hist(as.vector(M1)[which(as.vector(M1)>0)],breaks = seq(0,cluster_num,length.out = cluster_num),plot = F)$counts/(sum( hist(as.vector(M1)[which(as.vector(M1)>0)],breaks = seq(0,cluster_num,length.out = cluster_num),plot = F)$counts))
    M2_v <- hist(as.vector(M2)[which(as.vector(M2)>0)],breaks = seq(0,cluster_num,length.out = cluster_num),plot = F)$counts/(sum( hist(as.vector(M2)[which(as.vector(M2)>0)],breaks = seq(0,cluster_num,length.out = cluster_num),plot = F)$counts))
    weight_2 <-  1-max(min(distance(rbind(M1_v,M2_v),method = 'kullback-leibler',test.na = F,mute.message = T),1))

    M1 = tm_matrix_org[(ximgcol-spot_r_max):(ximgcol+spot_r_max),(ximgrow-spot_r_max):(ximgrow+spot_r_max)]
    M2 = tm_matrix_org[(yimgcol-spot_r_max):(yimgcol+spot_r_max),(yimgrow-spot_r_max):(yimgrow+spot_r_max)]
    if(min(M1)==max(M1)) next
    if(min(M2)==max(M2)) next
    M1_v <- as.vector(M1)[which(as.vector(M1)>0)]
    M2_v <- as.vector(M2)[which(as.vector(M2)>0)]
    mean_diff_max <- abs(mean(M1_v)-mean(M2_v))*255
    M1 = tm_matrix_org[(ximgcol-spot_r_min):(ximgcol+spot_r_min),(ximgrow-spot_r_min):(ximgrow+spot_r_min)]
    M2 = tm_matrix_org[(yimgcol-spot_r_min):(yimgcol+spot_r_min),(yimgrow-spot_r_min):(yimgrow+spot_r_min)]
    if(min(M1)==max(M1)) next
    if(min(M2)==max(M2)) next
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
  saveRDS(ST_net,file = paste0(savePath,"ST_imgunsupnet.RDS"))
  expr_data <- expr_obj_scale@assays$RNA@scale.data

  hist(E(ST_net)$weight[which(E(ST_net)$weight>0.01)])

  expr_info_matrix <- t(expr_data)
  search_map <- list(c(2,0),c(-1,-1),c(1,-1),c(-2,0),c(-1,1),c(1,1))
  if(Method=='walktrap'){
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
        ggsave(filename = paste0(savePath, "/ST_meta_imgunsup_plot_Step_",i,"_C",c_num[i],".png"), ST_meta_plot_tm,
               width = 8, height = 7, dpi = 500)
      }

      score_df <- cbind(c_num,step_testSDbw,I_Score,discrete_points_num)
      score_df <- as.data.frame(score_df)
      saveRDS(score_df,file = paste0(savePath,"score_df.rds"))

      score_df <- score_df[sort(score_df$c_num,index.return=TRUE,decreasing = T)$ix,]
      attach(score_df)
      tm <- aggregate(score_df, by=list(c_num = score_df$c_num), FUN = max)

      png(file =paste0(savePath,"Iscore.png"))
      plot(tm$c_num,tm$I_Score,type = "o",col = "red", xlab = "c_num", ylab = "I_Score",
           main = "I_Score")
      dev.off()

      png(file =paste0(savePath,"step_clusternum_imgunsup.png"))
      plot(tm$c_num,tm$c_num,type = "o",col = "blue", xlab = "c_num", ylab = "c_num",
           main = "c_num")
      dev.off()
      png(file =paste0(savePath,"SDbw_imgunsup.png"))
      plot(tm$c_num,tm$step_testSDbw,type = "o",col = "red", xlab = "c_num", ylab = "SDbw_img",
           main = "SDbw_img")
      dev.off()
      png(file =paste0(savePath,"discrete_points_num.png"))
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
    ggsave(filename = file.path(savePath, "/ST_meta_imgunsup_plot.png"), ST_meta_plot,
           width = 8, height = 7, dpi = 500)


    Spot_manifest$Walktrap_id <- as.character(Spot_manifest$Walktrap_id)
  }

  saveRDS(expr_obj_scale,file = paste0(savePath,"expr_obj.RDS"))

  saveRDS(Spot_manifest,file = paste0(savePath,"Spot_manifest_imgunsup.RDS"))
  write.csv(Spot_manifest,file = paste0(savePath,"Spot_manifest_imgunsup.csv"))
  return(Spot_manifest)
}


mirrorindex <- function(num,length){
  if(num<0){
    num = -num-1
  }
  if(num>=length){
    num = 2*length-num-1
  }
  return(num)
}

R_imageconv <- function(img,kernel){
  r = nrow(img);
  c = ncol(img);
  img_rst <- matrix(nrow=r,ncol=c);
  for(i in 1:r){
    for(j in 1:c){
      rst = 0;
      for(h in 0:2){
        for(k in 0:2){
          imagI = i+h;
          imagJ = j+k;
          imagI = mirrorindex(imagI,r);
          imagJ = mirrorindex(imagJ,c);
          rst = rst+ img[imagI,imagJ]*kernel[h+1,k+1];
        }
        if(rst <= 0){
          rst = 0;
        }
        if(rst>=1){
          rst = 1;
        }
        img_rst[i,j] = rst;
      }
    }
  }
  return(img_rst)
}


generate_mask <- function(imagefile,
                          savePath,
                          sacle_score,
                          spaceFile,
                          spot_r_max = 20){
  img <- load.image(imagefile)
  im <- grayscale(img) %>% isoblur(2)
  plot(im)

  Spot_space <- read.csv(file = spaceFile,header = F)
  sel.cols <- c("barcode","tissue", "row", "col", "imagerow", "imagecol")
  colnames(Spot_space) <- sel.cols
  Spot_space$barcode <- ST_filter_str(Spot_space$barcode,'-')
  mask <- zeros(dim(im)[1],dim(im)[2])
  for(i in 1:length(Spot_space$barcode)){
    ximgrow = round(Spot_space[i,"imagerow"]*sacle_score)
    ximgcol = round(Spot_space[i,"imagecol"]*sacle_score)
    mask[(ximgcol-spot_r_max):(ximgcol+spot_r_max),(ximgrow-spot_r_max):(ximgrow+spot_r_max)] = 1
  }
  write.table(mask,paste0(savePath,"mask1.txt"),quote = F,row.names = F,col.names = F,sep = ',')
}
