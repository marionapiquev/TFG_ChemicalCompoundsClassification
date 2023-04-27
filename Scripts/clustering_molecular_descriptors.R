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
library(matrixStats)
########################################
#LOAD THE DATAFRAME
#Set the directory 
setwd("/Users/MarionaP/Desktop/TFG/DrugBank")
Df_scaled <- read.csv("df_molecular_descriptors_scaled.csv")
df_original <- read.csv("df_molecular_descriptors.csv")
smiles <- read.csv("df_smiles.csv")
df_filter <- read.csv("df_molecular_descriptors_scaled_filtered.csv")
########################################
#PCA
pca <- function(df){
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
}

pca_DF <- pca(df)

p <- plot_ly(pca_DF, x=pca_DF$Comp.1, y=pca_DF$Comp.2, 
             z=pca_DF$Comp.3) %>% add_markers(size=1.5)
print(p)
########################################
#K-MEANS
#Prepare the data for clustering
DF <- Df_scaled #Choose the dataframe you want to use
cnames <- colnames(DF)
df <- DF[,cnames[!cnames %in% c("X","id")]]

#Since we are with a high dimension dataframe we use PCA to make the clustering easier
pca_DF <- pca(df)

#Perform the clustering
km.res <- kmeans(pca_DF, 10)
km.res

#Add clustering results to a dataframe using PCA dataframe
df_km <- pca_DF
df_km['cluster'] <- km.res['cluster']

#Visualize the clustering
fig <- plot_ly(df_km, x = ~Comp.1, y = ~Comp.2, z = ~Comp.3, color = ~cluster)
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'MLogP'),yaxis = list(title = 'OB_MW'),zaxis = list(title = 'OB_MR')))
fig

#Add clustering results to a dataframe using the original dataframe
df_km <- df_original
df_km['cluster'] <- km.res['cluster']

#Visualize the clustering
fig <- plot_ly(df_km, x = ~MLogP, y = ~OB_MW, z = ~OB_logP, color = ~cluster)
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'MLogP'),yaxis = list(title = 'OB_MW'),zaxis = list(title = 'OB_logP')))
fig

########################################
#HIERARCHICAL CLUSTERING
#Prepare the data for HC
DF <- Df_scaled #Choose the dataframe you want to use
cnames <- colnames(DF)
df <- DF[,cnames[!cnames %in% c("X","id","X.1")]]
dist_matrix <- daisy(df, metric = "euclidean")

#Perform the clustering
hc <- hclust(dist_matrix, method = 'ward.D')
plot(hc,main = paste('Dendrogram'), xlab = 'Molecules', ylab = 'Euclidean distances')

#Cut tree in n clusters
grp <- cutree(hc, k = 10)
plot(hc,main = paste('Dendrogram'),xlab = 'Molecules', ylab = 'Euclidean distances') 
rect.hclust(hc, k = 10, border = 2:5)

#Add clustering results to a dataframe
df_hc <- df_original
df_hc['cluster'] <- grp

#Visualize the clustering
fig <- plot_ly(df_hc, x = ~MLogP, y = ~OB_MW, z = ~OB_logP, color = ~cluster)
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'MLogP'),yaxis = list(title = 'OB_MW'),zaxis = list(title = 'OB_logP')))
fig

#Visualize the result using the pca dataframe 
df_hc <- pca_DF
df_hc['cluster'] <- grp

fig <- plot_ly(df_km, x = ~Comp.1, y = ~Comp.2, z = ~Comp.3, color = ~cluster)
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'MLogP'),yaxis = list(title = 'OB_MW'),zaxis = list(title = 'OB_MR')))
fig
#ggplot(data = df_hc, mapping = aes(x = MLogP, y = OB_MW,color = cluster))+geom_point() 

########################################
#DBSCAN
#Prepare the data for clustering
DF <- Df_scaled
cnames <- colnames(DF)
df <- DF[,cnames[!cnames %in% c("X","id","X.1")]]

#Perform the clustering
#BO: eps = 2, minPts = 10
#df_filter (small): eps = 5, minPts = 5
#df_filter (small, aromatic): eps = 6, minPts = 5
dbscan_result <- dbscan(df,eps = 2, minPts = 10)

#Add clustering results to a dataframe
df_dbscan <- df_original
df_dbscan['cluster'] <- dbscan_result$cluster

#Visualize the clustering
df_clusters <- filter(df_dbscan, cluster > 0)
ggplot(data = df_clusters, mapping = aes(x = MLogP,y = OB_MW, color = cluster))+geom_point() 

fig <- plot_ly(df_clusters, x = ~MLogP, y = ~OB_MW, z = ~OB_logP, color = ~cluster)
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'MLogP'),yaxis = list(title = 'OB_MW'),zaxis = list(title = 'OB_logP')))
fig

########################################
#HDBSCAN
#Prepare the data for clustering
DF <- Df_scaled
cnames <- colnames(DF)
df <- DF[,cnames[!cnames %in% c("X","id","X.1")]]

#Perform the clustering
h_dbscan_result <- hdbscan(df, minPts = 3)

#Visualize the Simplified Tree
plot(h_dbscan_result, show_flat = T)

#Add clustering results to a dataframe
df_hdbscan <- df_original
df_hdbscan['cluster'] <- h_dbscan_result$cluster

#Visualize the clustering
df_clusters <- filter(df_hdbscan, cluster > 0)
ggplot(data = df_clusters, mapping = aes(x = MLogP,y = OB_MW, color = cluster))+geom_point() 

fig <- plot_ly(df_clusters, x = ~MLogP, y = ~OB_MW, z = ~OB_logP, color = ~cluster)
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'MLogP'),yaxis = list(title = 'OB_MW'),zaxis = list(title = 'OB_logP')))
fig


