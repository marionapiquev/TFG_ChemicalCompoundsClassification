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
df_original <- read.csv("df_molecular_descriptors.csv")
smiles <- read.csv("df_smiles.csv")
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
grp <- cutree(hc, k = 20)
plot(hc,main = paste('Dendrogram'),xlab = 'Molecules', ylab = 'Euclidean distances') 
rect.hclust(hc, k = 20, border = 2:5)

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
df <- DF[,cnames[!cnames %in% c("X","id","X.1")]]

#Perform the clustering
#BO: eps = 2, minPts = 10
dbscan_result <- dbscan(df,eps = 2, minPts = 10)

#Add clustering results to a dataframe
df_dbscan <- df_original
df_dbscan['cluster'] <- dbscan_result$cluster

#Visualize the clustering
df_clusters <- filter(df_dbscan, cluster > 0)
ggplot(data = df_clusters, mapping = aes(x = MLogP,y = OB_MW, color = cluster))+geom_point() 

#Function to visualize molecules for each cluster
visualize_cluster_mols <- function(df_dbscan, cnames, nclusters){
  for (i in 1:nclusters){
    cluster1 <- filter(df_dbscan, cluster == i)
    c1_id <- cluster1[cnames[cnames %in% c("id")]]
    c1_smiles <- filter(smiles, id %in% as.list(c1_id)$id)
    c1_sdf <- c()
    for (i in as.list(c1_smiles['smile'])$smile){
      c1_sdf <- c(c1_sdf, smiles2sdf(i))
    }
    par(mfrow=c(4,4))
    for (i in 1:length(c1_sdf)){
      openBabelPlot(c1_sdf[[i]])
    }
  }
}
visualize_cluster_mols(df_dbscan, cnames,5)

#Function to get the similarity comparission of molecules on the same cluster
mols_similarity <- function(df_dbscan, cnames, nclusters){
  cluster1 <- filter(df_dbscan, cluster == 1)
  c1_id <- cluster1[cnames[cnames %in% c("id")]]
  c1_smiles <- filter(smiles, id %in% as.list(c1_id)$id)
  mol_similarity <- data.frame(mol1 = c1_smiles['smile']$smile[1], mol2 = c1_smiles['smile']$smile[1], Tanimoto_Coef = calcDrugMCSSim(c1_smiles['smile']$smile[1], c1_smiles['smile']$smile[1], type = 'smile')[[2]], Overlap_Coef = calcDrugMCSSim(c1_smiles['smile']$smile[1], c1_smiles['smile']$smile[1], type = 'smile')[[3]], cluster = 1)
  for (k in 1:nclusters){  
    cluster1 <- filter(df_dbscan, cluster == k)
    c1_id <- cluster1[cnames[cnames %in% c("id")]]
    c1_smiles <- filter(smiles, id %in% as.list(c1_id)$id)
    for (i in c1_smiles['smile']$smile){
      c_filter <- filter(c1_smiles, smile != i)
      for (j in c_filter['smile']$smile){
        sim <- calcDrugMCSSim(i, j, type = 'smile')
        mol_similarity[nrow(mol_similarity) + 1,] <- data.frame(mol1 = i, mol2 = j, Tanimoto_Coef = sim[[2]], Overlap_Coef = sim[[3]], cluster = k)
      }
    }
  }
  mol_similarity <- mol_similarity[-1,]
  return(mol_similarity)
}

df_mols_similarity <- mols_similarity(df_dbscan, cnames, 3)
