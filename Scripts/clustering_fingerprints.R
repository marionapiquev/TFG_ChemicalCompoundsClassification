########################################
#INSTALLING R PACKAGES (ONLY RON ONCE!!!)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("ChemmineR")
#BiocManager::install("Rcpi")

#library(devtools)
#install_github("CDK-R/rpubchem", dependencies=TRUE)
#devtools::install_github('rajarshi/chemblr/package')
########################################
#LOADING R PACKAGES
library(rcdk)
library(ChemmineR) 
library(rpubchem)
library(chemblr)
library(fingerprint)
library(Rcpi)
library(cluster)
library("factoextra")
library(fpc)
library(scatterplot3d)
library(gplots)
library(ggcorrplot)
library(dbscan)
########################################
#INCREASE THE SIZE OF JAVA RUNTIME HEAP
options(java.parameters = "-Xmx4g")
########################################
#IMPORT MOLECULES IN R 

#Set the directory
setwd("/Users/MarionaP/Desktop/TFG/DrugBank")

#Load molecules using a Database
mols <- load.molecules("DrugBank_3Dstructures.sdf")
sdfset <- read.SDFset("DrugBank_3Dstructures.sdf") # 9213 molecules

########################################
#GET FINGERPRINTS
#We extract the fingerprints of the compunds
fps <- lapply(mols, get.fingerprint, type='extended')

#Compute the intermolecular distances by the Tanimoto index
fp.sim <- fp.sim.matrix(fps, method='tanimoto')
fp.dist <- 1 - fp.sim

########################################
#K-MEANS
#Perform the clustering
km.res <- kmeans(fp.dist, 10)
km.res

#fviz_cluster(km.res, data = fp.dist, ellipse.type = "convex")

########################################
#HIERARCHICAL CLUSTERING
#Prepare the data for HC
d <- dist(fp.dist) 
distance <- as.dist(d) 

#Perform the clustering
res.hc <- hclust(distance, method = "ward.D2")
plot(res.hc,main = paste('Dendrogram'), xlab = 'Fingerprints', ylab = 'Distances')

#Cut tree in n clusters
grp <- cutree(res.hc, k = 4)
plot(res.hc, cex = 0.6) 
rect.hclust(res.hc, k = 4, border = 2:5)  

########################################
#STRUCTURAL CLUSTERING USING DBSCAN
result <- dbscan(fp.dist, eps = 0.5)
