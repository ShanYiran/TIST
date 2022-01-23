

#' ggplot_config
#'
#' @param base.size The size of text.
#'
#' @return A theme.
#' @export
#'
ggplot_config <- function(base.size = 8){
  p <- theme_classic() +
    theme(plot.title = element_text(size = 2 * base.size),
          axis.title.x = element_text(size = 2 * base.size, vjust = -0.2),
          axis.title.y = element_text(size = 2 * base.size, vjust = 0.2),
          axis.text.x = element_text(size = 2 * base.size),
          axis.text.y = element_text(size = 2 * base.size),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_text(size = 2 * base.size - 2),
          legend.text = element_text(size = 1.5 * base.size)
    )
  return(p)
}

##############reference################
getDefaultColors <- function(n = NULL, type = 1){
  if(type == 1){
    colors <- c("#cb7c77", "#68d359", "#6a7dc9", "#c9d73d", "#c555cb",
                "#d7652d", "#7cd5c8", "#c49a3f", "#507d41", "#5d8d9c",
                "#90353b", "#674c2a", "#1B9E77", "#c5383c", "#0081d1",
                "#ffd900", "#502e71", "#c8b693", "#aed688", "#f6a97a",
                "#c6a5cc", "#798234", "#6b42c8", "#cf4c8b", "#666666",
                "#feb308", "#ff1a1a", "#1aff1a", "#1a1aff", "#ffff1a")
  }else if(type == 2){
    if(n <= 8){
      colors <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
                  "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")
    }else if(n <= 14){
      colors <- c("#437BFE", "#FEC643", "#43FE69", "#FE6943", "#C643FE",
                  "#43D9FE", "#B87A3D", "#679966", "#993333", "#7F6699",
                  "#E78AC3", "#333399", "#A6D854", "#E5C494")
    }
    else if(n <= 20){
      colors <- c("#87b3d4", "#d5492f", "#6bd155", "#683ec2", "#c9d754",
                  "#d04dc7", "#81d8ae", "#d34a76", "#607d3a", "#6d76cb",
                  "#ce9d3f", "#81357a", "#d3c3a4", "#3c2f5a", "#b96f49",
                  "#4e857e", "#6e282c", "#d293c8", "#393a2a", "#997579")
    }else if(n <= 30){
      colors <- c("#628bac", "#ceda3f", "#7e39c9", "#72d852", "#d849cc",
                  "#5e8f37", "#5956c8", "#cfa53f", "#392766", "#c7da8b",
                  "#8d378c", "#68d9a3", "#dd3e34", "#8ed4d5", "#d84787",
                  "#498770", "#c581d3", "#d27333", "#6680cb", "#83662e",
                  "#cab7da", "#364627", "#d16263", "#2d384d", "#e0b495",
                  "#4b272a", "#919071", "#7b3860", "#843028", "#bb7d91")
    }else{
      colors <- c("#982f29", "#5ddb53", "#8b35d6", "#a9e047", "#4836be",
                  "#e0dc33", "#d248d5", "#61a338", "#9765e5", "#69df96",
                  "#7f3095", "#d0d56a", "#371c6b", "#cfa738", "#5066d1",
                  "#e08930", "#6a8bd3", "#da4f1e", "#83e6d6", "#df4341",
                  "#6ebad4", "#e34c75", "#50975f", "#d548a4", "#badb97",
                  "#b377cf", "#899140", "#564d8b", "#ddb67f", "#292344",
                  "#d0cdb8", "#421b28", "#5eae99", "#a03259", "#406024",
                  "#e598d7", "#343b20", "#bbb5d9", "#975223", "#576e8b",
                  "#d97f5e", "#253e44", "#de959b", "#417265", "#712b5b",
                  "#8c6d30", "#a56c95", "#5f3121", "#8f846e", "#8f5b5c")
    }
  }else if(type == 3){
    # colors <- c("#07a2a4", "#9a7fd1", "#588dd5", "#f5994e",
    #             "#c05050", "#59678c", "#c9ab00", "#7eb00a")
    colors <- c("#c14089", "#6f5553", "#E5C494", "#738f4c", "#bb6240",
                "#66C2A5", "#2dfd29", "#0c0fdc")
  }
  if(!is.null(n)){
    if(n <= length(colors)){
      colors <- colors[1:n]
    }else{
      step <- 16777200 %/% (n - length(colors)) - 2
      add.colors <- paste0("#", as.hexmode(seq(from = sample(1:step, 1),
                                               by = step, length.out = (n-length(colors)))))
      colors <- c(colors, add.colors)
    }
  }
  return(colors)
}

getColors <- function(n){
  if(n==6){
    colors <- c("#cb7c77", "#68d359", "#6a7dc9", "#c9d73d", "#c555cb","#f5f5f5")
  }
  if(n==5){
    colors <- c("#cb7c77", "#68d359", "#6a7dc9", "#c9d73d", "#f5f5f5")
  }
  if(n==7){
    colors <- c("#cb7c77", "#68d359", "#6a7dc9", "#c9d73d", "#d7652d", "#7cd5c8", "#f5f5f5")
  }
  if(n>7){
    colors <- c("#cb7c77", "#68d359", "#6a7dc9", "#c9d73d", "#f5f5f5","#c555cb",
                "#d7652d", "#7cd5c8", "#c49a3f", "#507d41", "#5d8d9c",
                "#90353b", "#674c2a", "#1B9E77", "#c5383c", "#0081d1",
                "#ffd900", "#502e71", "#c8b693", "#aed688", "#f6a97a",
                "#c6a5cc", "#798234", "#6b42c8", "#cf4c8b", "#666666",
                "#feb308", "#ff1a1a", "#1aff1a", "#1a1aff", "#ffff1a") 
    
  }
  return(colors)
}

ST_filter_str <- function(data,strsep = '_',filter_id = 1){
  filter_barcode <- strsplit(data,strsep)
  filter_barcode_ <- rep(1,length(filter_barcode))
  for(i in 1:length(filter_barcode_)){
    filter_barcode_[i] <- filter_barcode[[i]][filter_id]
  }
  return(filter_barcode_)
}
Renamedomin <- function(namelist){
  vis <- rep(0,5000)
  cnt = 1
  namelist <- namelist+1
  for(i in 1:length(namelist)){
    if(vis[namelist[i]]==0){
      vis[namelist[i]] = cnt
      cnt = cnt + 1
    }
    namelist[i] = vis[namelist[i]]
  }
  return(namelist)
}
Renamedomin0 <- function(namelist){
  vis <- rep(0,5000)
  vis[1] <- 1
  cnt = 2
  namelist <- namelist+1
  for(i in 1:length(namelist)){
    if(vis[namelist[i]]==0){
      if(namelist[[i]]!=1){
        vis[namelist[i]] = cnt
        cnt = cnt + 1
      }
    }
    namelist[i] = vis[namelist[i]]
  }
  return(namelist)
}