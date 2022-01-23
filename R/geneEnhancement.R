geneEnhancement <- function(SC_obj,netfile,Spot_manifest,savePath,gene_name,gene_plot = F,width=1600, height=350,pointsize = 3,res = 100){
  ST_net <- readRDS(file = netfile)
  colnames(SC_obj@assays$RNA@data) <- ST_filter_str(colnames(SC_obj),'-')

  png(filename=paste0(savePath,"geneEnhancement_plot.png"),  width=width, height=height*ceiling(length(gene_name)/4),res = res*ceiling(length(gene_name)/4))
  pp <- lapply(1:length(gene_name), function(x) {
    gn <- gene_name[x]
    expr_data <- SC_obj@assays$RNA@data[gn,]
    names(expr_data) <- colnames(SC_obj@assays$RNA@data)
    x_o <- (expr_data+0.1)
    x_tm <- x_o
    nodes <- V(ST_net)[Spot_manifest$barcode]
    g <- ST_net
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
    message(gn)
    pltdat <- cbind(Spot_manifest[, c("imagerow","imagecol")],x_o)
    colnames(pltdat)[1:2] <- c("y","x")
    pltdat$y = -pltdat$y
    pal <- colorRampPalette(c("mediumseagreen", "lightyellow2", "deeppink"))
    xpand = 0
    ypand = 1
    ggplot(pltdat, aes(x = x, y = y, color = x_o)) + geom_point(size = pointsize) +
      # scale_color_gradientn(colours=pal(5))+
      scale_color_gradientn(colours = pal(5)) + scale_x_discrete(expand = c(xpand, ypand)) + scale_y_discrete(expand = c(xpand, ypand)) + coord_equal() +
      labs(title = gn, x = NULL, y = NULL) +
      # labs(title = colnames(pd)[igene+2], x = NULL, y = NULL)+
      theme_bw()+guides(fill=guide_legend(title=" "))
  })
  grid.arrange(grobs = pp, ncol = 4)
  dev.off()


  png(filename=paste0(savePath,"no_enhance.png"),  width=width, height=height*ceiling(length(gene_name)/4),res = res*ceiling(length(gene_name)/4))
  pp <- lapply(1:length(gene_name), function(x) {
    gn <- gene_name[x]
    x_o <- SC_obj@assays$RNA@data[gn,]
    pltdat <- cbind(Spot_manifest[, c("imagerow","imagecol")],x_o)
    colnames(pltdat)[1:2] <- c("y","x")
    pltdat$y = -pltdat$y
    pal <- colorRampPalette(c("mediumseagreen", "lightyellow2", "deeppink"))
    xpand = 0
    ypand = 1
    ggplot(pltdat, aes(x = x, y = y, color = x_o)) + geom_point(size = pointsize) +
      # scale_color_gradientn(colours=pal(5))+
      scale_color_gradientn(colours = pal(5)) + scale_x_discrete(expand = c(xpand, ypand)) + scale_y_discrete(expand = c(xpand, ypand)) + coord_equal() +
      labs(title = gn, x = NULL, y = NULL) +
      # labs(title = colnames(pd)[igene+2], x = NULL, y = NULL)+
      theme_bw()+guides(fill=guide_legend(title=" "))
  })
  grid.arrange(grobs = pp, ncol = 4)
  dev.off()


}


