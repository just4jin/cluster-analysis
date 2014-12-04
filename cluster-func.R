# Jin Li
# cluster analysis functions

## function remove outliers from data sets
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

## function runs the K-Means algorithm on data matrix X
runkMeans <- function(X, initial_centroids, max_iters) {
  # Initialize values
  m <- nrow(X)
  n <- ncol(X)
  K <- nrow(initial_centroids)
  centroids <- initial_centroids
  idx = matrix(0, m, 1)
  
## Run K-Means
  for (i in 1:max_iters){
    idx <- findClosestCentroids(X, centroids)
    centroids = computeCentroids(X, idx, K)
  }
  
  return(list(centroids = centroids, idx = idx))
}

## function computes the centroid memberships for every example
findClosestCentroids <- function(X, centroids){
  K <- nrow(centroids)
  m <- nrow(X)
  idx <- matrix(0, m, 1)
  tmp_d <- matrix(0, K, 1);
  for (i in 1:m){
    x <- X[i,]
    for (j in 1:K){
      tmp_d[j] = sqrt(sum((x - centroids[j,]) ^ 2))
    }
    idx[i] <- which.min(tmp_d)
  }
  idx
}

## function returs the new centroids by computing the means of the data points assigned to each centroid
computeCentroids <- function(X, idx, K){
  m <- nrow(X)
  n <- ncol(X)
  centroids <- matrix(0, K, n)
  for (j in 1:K){
    points_set <- 1 * (idx == j)
    centroids[j,] = colSums(X * points_set[,rep(1,n)]) / sum(points_set)
  }
  centroids
}
