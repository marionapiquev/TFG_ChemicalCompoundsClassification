########################################
options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8192m"))
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
library(dplyr)
library(RxnSim)
library(hrbrthemes)
########################################
#LOAD THE DATAFRAME
#Set the directory 
setwd("/Users/MarionaP/Desktop/TFG/DrugBank")
Df_scaled <- read.csv("df_molecular_descriptors_scaled.csv")
df_original <- read.csv("df_molecular_descriptors.csv")
smiles <- read.csv("df_smiles.csv")
df_filter <- read.csv("df_molecular_descriptors_scaled_filtered.csv")
########################################
#Function that computes Tanimoto coeficient using a similarity matrix.
mols_similarity_matrix <- function(df_dbscan, cnames, nclusters){
  for (k in 1:nclusters){
    cluster <- filter(df_dbscan, cluster == k)
    c_id <- cluster[cnames[cnames %in% c("id")]]
    c_smiles <- filter(smiles, id %in% as.list(c_id)$id)
    m <- ms.compute.sim.matrix(c_smiles$smile, format = 'SMILES', fp.type='pubchem', sim.method='tanimoto')
    #print(m)
    m[lower.tri(m, diag=TRUE)] <- 0
    similarity <- sum(m)/sum(rowCounts(m > 0))
    cat('\nTanimoto coef mean cluster',k,':',similarity)
  }
}

#Function that calculates the performance of the HDBSCAN algorithm for different parameters values based on the Tanimoto Similarity
parameters_hc <- function(dist_matrix, numclusters){
  mol_similarity_hc <- data.frame(k = integer(), nclusters = integer(), non_class = integer(), cluster_max_size = integer(), Similarity = double())
  hc <- hclust(dist_matrix, method = 'ward.D')
  for (j in numclusters){ 
    print(j)
    unclass <- 0
    total <- 0
    grp <- cutree(hc, k = j)
    df_hc <- df_original
    df_hc['cluster'] <- grp
    max_size <- df_hc %>% count(cluster)
    c_max_size <- max(max_size['n'])
    print(c_max_size)
    nclusters <- length(unique(df_hc$cluster))
    if (nclusters > 0){
      for (k in 1:nclusters){
        cluster <- filter(df_hc, cluster == k)
        if (nrow(cluster) < 3000 && nrow(cluster) > 1){
          c_id <- cluster[cnames[cnames %in% c("id")]]
          c_smiles <- filter(smiles, id %in% as.list(c_id)$id)
          m <- ms.compute.sim.matrix(c_smiles$smile, format = 'SMILES', fp.type='pubchem', sim.method='tanimoto')
          m[lower.tri(m, diag=TRUE)] <- 0
          similarity <- sum(m)/sum(rowCounts(m > 0))
          if(similarity == 'NaN'){total <- total + 0 } #if similarity is 0 the result is NaN
          else{total <- total + similarity}
        }
        else if (nrow(cluster) == 1){uncalss <- uncalss + 1;print(uncalss)}#total <- total + 0}
        else{
          total <- total + 0
        }
      }
      mol_similarity_hc[nrow(mol_similarity_hc) + 1,] <- data.frame(k = j, nclusters = length(unique(df_hc$cluster)), non_class = unclass, cluster_max_size = c_max_size,  Similarity = total/nclusters)
    }
    else{
      mol_similarity_hc[nrow(mol_similarity_hc) + 1,] <- data.frame(k = j, nclusters = length(unique(df_hc$cluster)), non_class = unclass, cluster_max_size = c_max_size, Similarity = 0)
    }
  }  
  return(mol_similarity_hc)
}
########################################
#HIERARCHICAL CLUSTERING
#Prepare the data for HC
DF <- Df_scaled #Choose the dataframe you want to use
cnames <- colnames(DF)
df <- DF[,cnames[!cnames %in% c("X","id")]]
dist_matrix <- daisy(df, metric = "euclidean")

#Perform the clustering
hc <- hclust(dist_matrix, method = 'ward.D')
plot(hc,main = paste('Dendrogram'), xlab = 'Molecules', ylab = 'Euclidean distances')
library(dendextend)
dend1 <- as.dendrogram(hc)
max_depth(dend1)
#Cut tree in n clusters
grp <- cutree(hc, k = 600)
plot(hc,main = paste('Dendrogram'),xlab = 'Molecules', ylab = 'Euclidean distances') 
rect.hclust(hc, k = 50, border = 2:5)

#Add clustering results to a dataframe
df_hc <- df_original
df_hc['cluster'] <- grp

#Number of mols for cluster
cluster_dist <- df_hc %>% group_by(cluster) %>% summarise(total_count=n(), .groups = 'drop')
cluster_dist

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

#Tanimoto coeficient
nclusters <- length(unique(df_hc$cluster))
mols_similarity_matrix(df_hc, cnames, nclusters)

numclusters <- seq(100, 1500, by=100)
param_similarity <- parameters_hc(dist_matrix, numclusters)
write.csv(param_similarity, "param.csv")

plot_similarity <- ggplot(param_similarity, aes(k, Similarity)) +
  geom_point(color="steelblue")+geom_line(color="steelblue")+theme_ipsum()+ ggtitle("Mean Tanimoto Similarity") 
plot_similarity
