cluster-analysis
================
cluster analysis on yeast data sets from UCI machine learning repository https://archive.ics.uci.edu/ml/datasets/Yeast


X - matrix with 8 variables and attribute values 

initialize K to 10 and randomly selected K rows without replacement out of X and stored in initial_centroids as basis for running K-means iterations



function runMeans - takes data matrix X, initial centroids randomly generated and the maximum of iterations assigned as parameters, and used result generated from findClosetCentroids and computeCentroids to return a list of computed centroids

function findClosestCentroids - takes data matrix X and centroids to return the centroid memberships for every example

function computeCentroids - takes data matrix X, the result of findClosetCentroids and K to generate the new centroids by computing the means of the data points assigned to each centroid
