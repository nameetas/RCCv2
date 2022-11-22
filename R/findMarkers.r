#' @title Find markers from specific clusters
#' @description Find markers using normalize data in each cluster
#' @param expr log2 CPM matrix with rows as genes and columns as cells
#' @param info meta data
#' @param fields meta data column names to be used as annotations for marker identification
#' @return markers file
#' @example
#' findMarkersAll('log2cpm.rds', 'meta.csv', c('Major.celltype', 'celltype', 'sample'))
#'
#'

findMarkersAll = function(expr, info, fields){

  infoFile = info
  data = readRDS(expr)
  library(scran)

  info = read.csv(infoFile)
  rownames(info) = info[,1]
  info = info[,-1]

  if(identical(colnames(data), rownames(info))){
    print("data is ordered")
  } else {
    data = data[, order(colnames(data))]
    info = info[order(rownames(info)), ]
  }

  num = which(colnames(info) %in% fields)

  for(i in 1:length(num)){
    lev1 = colnames(info)[num[i]]

    a = findMarkers(as.matrix(data), info[, lev1], direction = "up", lfc = 1)
    mat = NULL

    for(i in 1:length(a)){
      print(i)
      clusterInfo = a[[i]]
      clusterInfo = subset(clusterInfo, clusterInfo$FDR < 0.01)
      #lfc_mat = as.data.frame(clusterInfo[,4:ncol(clusterInfo)])
      #lfc = apply(lfc_mat, 1, min)
      #clusterInfo = cbind(clusterInfo, lfc)
      #clusterInfo = subset(clusterInfo, clusterInfo$lfc > 1)
      #clusterInfo = subset(clusterInfo, clusterInfo$lfc > 1)

      if(nrow(clusterInfo) == 0){
        mat = qpcR:::cbind.na(mat, rep("0", 5))
      } else {
        mat = qpcR:::cbind.na(mat, rownames(clusterInfo))
      }
    }

    mat = mat[,-1]
    colnames(mat) = names(a)
    write.csv(mat, file = paste0(lev,"_markers.csv"))
  }

  files = list.files(".", pattern = "*markers.csv")
  mat = NULL

  for(i in 1:length(files)){
    data = read.csv(files[i])
    data = data[,-1]
    colnames(data) = gsub("X", "K", colnames(data))
    if(i == 1){
      mat = data
    } else {
      mat = qpcR:::cbind.na(mat, data)
    }
  }
  num = grep("_0$", colnames(mat))
  mat = mat[,-c(num)]
  write.csv(mat, file = "combineMarkers.csv")
}
