#' @title Cluster annotation
#' @description This function allows the users to compare different annotation classes to see the overlap between two meta fields
#' @param file meta data file
#' @param attribute Field one from meta data
#' @param clus Field two from meta data
#' @param outputDir output directory for cluster annotation file
#' @return PDF file with overlap between two annotation fields
#' @example
#' clusterAnnotation('output.csv', 'Major.celltype', 'celltype', '.')
#'

clusterAnnotation = function(file, attribute, clus, outputDir){
  library(matrixStats)
  library(ComplexHeatmap)
  library(circlize)
  library(clue)

  pMatrix.min <- function(A, B){
    norm
    n <- nrow(A)
    D <- matrix(NA, n, n)
    for (i in 1:n) {
      for (j in 1:n) {
        D[j, i] <- (sum((B[j, ] - A[i, ])^2))
      }
    }
    vec <- c(solve_LSAP(D))
    list(A = A[vec,], pvec=vec)
  }

  setwd(outputDir)
  sampleInfo = read.csv(file)
  colNum = which(colnames(sampleInfo) == attribute)

  mat = matrix(0, nrow = length(unique(sampleInfo[, clus])), ncol = length(unique(sampleInfo[,colNum])))
  colnames(mat) = unique(sampleInfo[,colNum])
  rownames(mat) = unique(sampleInfo[,clus])

  attr_df = as.data.frame(table(sampleInfo[,colNum]))
  rownames(attr_df) = attr_df[,1]

  for(i in 1:nrow(mat)){
    b = subset(sampleInfo, sampleInfo[,clus] == rownames(mat)[i])
    df = as.data.frame(table(b[,colNum]))
    rownames(df) = df[,1]
    for(j in 1:nrow(df)){
      x = rownames(df)[j]
      mat[i,x] = (df[j,2]/attr_df[x,2])*100
    }
  }
  B <- diag(1, nrow(mat))
  X <- pMatrix.min(mat,B)
  d = X$A

  pdf(paste0(attribute, "_vs_", clus,".pdf"),height = 15, width = 10)
  h1 = Heatmap(d, col = colorRamp2(c(0,50,100), c("white","orangered","red")),
               na_col = "black", row_names_gp = gpar(fontsize = 12),
               cluster_rows=F, cluster_columns=F)
  draw(h1)
  dev.off()
}
