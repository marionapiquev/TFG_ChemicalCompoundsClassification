########################################
#INSTALLING R PACKAGES (ONLY RUN ONCE!!!)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("ChemmineR")
#BiocManager::install("Rcpi")
#BiocManager::install("ChemmineOB")

#library(devtools)
#install_github("CDK-R/rpubchem", dependencies=TRUE)
#devtools::install_github('rajarshi/chemblr/package')
########################################
#LOADING R PACKAGES
library(rcdk)
library(ChemmineR) 
library(ChemmineOB) 
library(rpubchem)
library(chemblr)
library(fingerprint)
library(Rcpi)
library(cluster)
library("factoextra")
library(fpc)
library(ggcorrplot)
library(plotly)
library(ggplot2)
library(dbscan)
########################################
#LOAD THE DATAFRAME
#Set the directory 
setwd("/Users/MarionaP/Desktop/TFG/DrugBank")
DF <- read.csv("df_molecular_descriptors_scaled.csv")

########################################
#K-MEANS
#Prepare the data for clustering
cnames <- colnames(DF)
df <- DF[,cnames[!cnames %in% c("X","id")]]

#Perform the clustering
km.res <- kmeans(pca_DF, 10)
km.res

#Add clustering results to a dataframe
df_km <- pca_DF
df_km['cluster'] <- km.res['cluster']

#Visualize the clustering
fig <- plot_ly(df_km, x = ~Comp.1, y = ~Comp.2, z = ~Comp.3, color = ~cluster)
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'MLogP'),yaxis = list(title = 'OB_MW'),zaxis = list(title = 'OB_MR')))
fig

#PCA
pca.DF <- princomp(df) 
summary(pca.DF)
#Show the percentage of variances explained by each principal component.
fviz_eig(pca.DF)
#Graph of individuals (Individuals with a similar profile are grouped together)
fviz_pca_ind(pca.DF,
             col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
pca_DF <- as.data.frame(pca.DF$scores)
p <- plot_ly(pca_DF, x=pca_DF$Comp.1, y=pca_DF$Comp.2, 
             z=pca_DF$Comp.3) %>%
  add_markers(size=1.5)
print(p)

########################################
#HIERARCHICAL CLUSTERING
#Prepare the data for HC
cnames <- colnames(DF)
df <- DF[,cnames[!cnames %in% c("X","id","X.1")]]
dist_matrix <- daisy(df, metric = "euclidean")

#Perform the clustering
hc <- hclust(dist_matrix, method = 'ward.D')
plot(hc,main = paste('Dendrogram'), xlab = 'Molecules', ylab = 'Euclidean distances')

#Cut tree in n clusters
grp <- cutree(hc, k = 4)
plot(hc,main = paste('Dendrogram'),xlab = 'Molecules', ylab = 'Euclidean distances') 
rect.hclust(hc, k = 4, border = 2:5)

#Add clustering results to a dataframe
df_hc <- DF
df_hc['cluster'] <- grp

#Visualize the clustering
fig <- plot_ly(df_hc, x = ~MLogP, y = ~OB_MW, z = ~OB_logP, color = ~cluster)
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'MLogP'),yaxis = list(title = 'OB_MW'),zaxis = list(title = 'OB_logP')))
fig

ggplot(data = df_hc, mapping = aes(x = MLogP, y = OB_MW,color = cluster))+geom_point() 
########################################
#DBSCAN
#Prepare the data for clustering
cnames <- colnames(DF)
df <- DF[,cnames[!cnames %in% c("X","id")]]

#Perform the clustering
dbscan_result <- dbscan(df, eps = 2, minPts = 10)

#Add clustering results to a dataframe
df_dbscan <- DF
df_dbscan['cluster'] <- dbscan_result$cluster

#Visualize the clustering
df_clusters <- filter(df_dbscan, cluster > 0)
ggplot(data = df_clusters, mapping = aes(x = MLogP,y = OB_MW, color = cluster))+geom_point() 