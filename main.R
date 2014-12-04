##*****************************************************************
# Jin Li
# cluster analysis on yeast data sets
##*****************************************************************

# Source functions from the func file
source("./cluster-func.R")

# Read yeast.data 
yeast <- read.table("yeast.data", header=FALSE)

# Add column names
colnames(yeast) <- c("seq_name", "mcg", "gvh", "alm", "mit", "erl", "pox", "vac", "nuc", "class")

# Plot histogram and boxplot for each feature, to examine normality and locate outliers
# erl and pox are binary attributes, which will be droped out from the following analysis.
attach(yeast)
par(mfrow=c(4,4))
hist(mcg, main = "Histogram of mcg")
hist(gvh, main = "Histogram of gvh")
hist(alm, main = "Histogram of alm")
hist(mit, main = "Histogram of mit")
hist(erl, main = "Histogram of erl") #binary attribute (0.5 and 1)
hist(pox, main = "Histogram of pox") #attribute with discrete values
hist(vac, main = "Histogram of vac")
hist(nuc, main = "Histogram of nuc")

boxplot(mcg, main="Boxplot of mcg")
boxplot(gvh, main="Boxplot of gvh")
boxplot(alm, main="Boxplot of alm")
boxplot(mit, main="Boxplot of mit")
boxplot(erl, main="Boxplot of erl") #no outliers
boxplot(pox, main="Boxplot of pox") #no outliers
boxplot(vac, main="Boxplot of vac")
boxplot(nuc, main="Boxplot of nuc")

detach(yeast)

# Try 3 different transformation methods to get the most normalised features
# Statistical transformation
stat_yeast <- yeast
stat_yeast$mcg<-(stat_yeast$mcg-mean(stat_yeast$mcg))/sd(stat_yeast$mcg)
stat_yeast$gvh<-(stat_yeast$gvh-mean(stat_yeast$gvh))/sd(stat_yeast$gvh)
stat_yeast$alm<-(stat_yeast$alm-mean(stat_yeast$alm))/sd(stat_yeast$alm)
stat_yeast$mit<-(stat_yeast$mit-mean(stat_yeast$mit))/sd(stat_yeast$mit)
stat_yeast$vac<-(stat_yeast$vac-mean(stat_yeast$vac))/sd(stat_yeast$vac)
stat_yeast$nuc<-(stat_yeast$nuc-mean(stat_yeast$nuc))/sd(stat_yeast$nuc)

# Log transformation y=ln(xk-c),c=floor(min(xk))
log_yeast <- yeast
log_yeast$mcg<-log(log_yeast$mcg-floor(min(log_yeast$mcg)))
log_yeast$gvh<-log(log_yeast$gvh-floor(min(log_yeast$gvh)))
log_yeast$alm<-log(log_yeast$alm-floor(min(log_yeast$alm)))
log_yeast$mit<-log(log_yeast$mit-floor(min(log_yeast$mit)))
log_yeast$vac<-log(log_yeast$vac-floor(min(log_yeast$vac)))
log_yeast$nuc<-log(log_yeast$nuc-floor(min(log_yeast$nuc)))

# Standardized transformation
stan_yeast <- yeast
stan_yeast$mcg<-(stan_yeast$mcg-min(stan_yeast$mcg))/(max(stan_yeast$mcg)-min(stan_yeast$mcg))
stan_yeast$gvh<-(stan_yeast$gvh-min(stan_yeast$gvh))/(max(stan_yeast$gvh)-min(stan_yeast$gvh))
stan_yeast$alm<-(stan_yeast$alm-min(stan_yeast$alm))/(max(stan_yeast$alm)-min(stan_yeast$alm))
stan_yeast$mit<-(stan_yeast$mit-min(stan_yeast$mit))/(max(stan_yeast$mit)-min(stan_yeast$mit))
stan_yeast$vac<-(stan_yeast$vac-min(stan_yeast$vac))/(max(stan_yeast$vac)-min(stan_yeast$vac))
stan_yeast$nuc<-(stan_yeast$nuc-min(stan_yeast$nuc))/(max(stan_yeast$nuc)-min(stan_yeast$nuc))

# Display and compare histograms and boxplots for features transformed by different methods. Decide which method to use for each feature.
par(mfrow=c(4,4))
hist(yeast$mcg,main="Histogram of mcg")
hist(stat_yeast$mcg, main="Histogram of statistical transformed mcg")
hist(log_yeast$mcg, main="Histogram of log transformed mcg")
hist(stan_yeast$mcg, main="Histogram of standardized mcg")
boxplot(yeast$mcg, main="Boxplot of mcg")
boxplot(stat_yeast$mcg, main="Boxplot of statistical transformed mcg")
boxplot(log_yeast$mcg, main="Boxplot of log transformed mcg")
boxplot(stan_yeast$mcg, main="Boxplot of standardized mcg")
# Log transformation worked best for mcg

hist(yeast$gvh,main="Histogram of gvh")
hist(stat_yeast$gvh, main="Histogram of statistical transformed gvh")
hist(log_yeast$gvh, main="Histogram of log transformed gvh")
hist(stan_yeast$gvh, main="Histogram of standardized gvh")
boxplot(yeast$gvh, main="Boxplot of gvh")
boxplot(stat_yeast$gvh, main="Boxplot of statistical transformed gvh")
boxplot(log_yeast$gvh, main="Boxplot of log transformed gvh")
boxplot(stan_yeast$gvh, main="Boxplot of standardized gvh")
# Log transformation worked best for gvh

hist(yeast$alm,main="Histogram of alm")
hist(stat_yeast$alm, main="Histogram of statistical transformed alm")
hist(log_yeast$alm, main="Histogram of log transformed alm")
hist(stan_yeast$alm, main="Histogram of standardized alm")
boxplot(yeast$alm, main="Boxplot of alm")
boxplot(stat_yeast$alm, main="Boxplot of statistical transformed alm")
boxplot(log_yeast$alm, main="Boxplot of log transformed alm")
boxplot(stan_yeast$alm, main="Boxplot of standardized alm")
# Statistical transformation worked best for alm

hist(yeast$mit,main="Histogram of mit")
hist(stat_yeast$mit, main="Histogram of statistical transformed mit")
hist(log_yeast$mit, main="Histogram of log transformed mit")
hist(stan_yeast$mit, main="Histogram of standardized mit")
boxplot(yeast$mit, main="Boxplot of mit")
boxplot(stat_yeast$mit, main="Boxplot of statistical transformed mit")
boxplot(log_yeast$mit, main="Boxplot of log transformed mit")
boxplot(stan_yeast$mit, main="Boxplot of standardized mit")
# Log transformation worked best for mit

hist(yeast$vac,main="Histogram of vac")
hist(stat_yeast$vac, main="Histogram of statistical transformed vac")
hist(log_yeast$vac, main="Histogram of log transformed vac")
hist(stan_yeast$vac, main="Histogram of standardized vac")
boxplot(yeast$vac, main="Boxplot of vac")
boxplot(stat_yeast$vac, main="Boxplot of statistical transformed vac")
boxplot(log_yeast$vac, main="Boxplot of log transformed vac")
boxplot(stan_yeast$vac, main="Boxplot of standardized vac")
# Log transformation worked best for vac

hist(yeast$nuc,main="Histogram of nuc")
hist(stat_yeast$nuc, main="Histogram of statistical transformed nuc")
hist(log_yeast$nuc, main="Histogram of log transformed nuc")
hist(stan_yeast$nuc, main="Histogram of standardized nuc")
boxplot(yeast$nuc, main="Boxplot of nuc")
boxplot(stat_yeast$nuc, main="Boxplot of statistical transformed nuc")
boxplot(log_yeast$nuc, main="Boxplot of log transformed nuc")
boxplot(stan_yeast$nuc, main="Boxplot of standardized nuc")
# Log transformation worked best for nuc

trans_yeast <- log_yeast
trans_yeast$alm <- stat_yeast$alm

# Remove outliers
trans_yeast$mcg<-remove_outliers(trans_yeast$mcg)
trans_yeast$gvh<-remove_outliers(trans_yeast$gvh)
trans_yeast$alm<-remove_outliers(trans_yeast$alm)
trans_yeast$mit<-remove_outliers(trans_yeast$mit)
trans_yeast$vac<-remove_outliers(trans_yeast$vac)
trans_yeast$nuc<-remove_outliers(trans_yeast$nuc)
clean_yeast<-na.omit(trans_yeast)

# Determine linear correlations between features
cor(clean_yeast[,2:9])
# There is no strong linear correlation between features in our data set.

# Clustering analysis. We are going to use K-means clustering method, one type of prototype based clustering method.
# Randomely initialize 10 centroids
X <- clean_yeast[,2:9]
K <- 10
	
initial_centroids <- X[sample(nrow(X), size = K, replace = FALSE),]
# Run K-Means over 10 iterations
max_iters <- 10
result <- runkMeans(X, initial_centroids, max_iters)
# Graphical display of clusters
par(mfrow=c(1,1))
plot(X, col=result$idx) #scatterplot matrix
table(result$idx) #how many points belong to each cluster
result
pdf("output.pdf")
plot(X[c(3,4)],col=result$idx)
result
cat("Centroid 1 is marked in red on the scatter plot, which has", sum(1*(result$idx==1)), "points centered around it.\n")
cat("Centroid 2 is marked in blue on the scatter plot, which has", sum(1*(result$idx==2)), "points centered around it.\n")
cat("Centroid 3 is marked in green on the scatter plot, which has", sum(1*(result$idx==3)), "points centered around it.\n")
cat("Centroid 4 is marked in yellow on the scatter plot, which has", sum(1*(result$idx==4)), "points centered around it.\n")
cat("Centroid 5 is marked in purple on the scatter plot, which has", sum(1*(result$idx==5)), "points centered around it.\n")
cat("Centroid 6 is marked in red on the scatter plot, which has", sum(1*(result$idx==6)), "points centered around it.\n")
cat("Centroid 7 is marked in blue on the scatter plot, which has", sum(1*(result$idx==7)), "points centered around it.\n")
cat("Centroid 8 is marked in green on the scatter plot, which has", sum(1*(result$idx==8)), "points centered around it.\n")
cat("Centroid 9 is marked in yellow on the scatter plot, which has", sum(1*(result$idx==9)), "points centered around it.\n")
cat("Centroid 10 is marked in purple on the scatter plot, which has", sum(1*(result$idx==10)), "points centered around it.\n")
result$centroids
sink("centro.pdf",append=FALSE, split=FALSE)
