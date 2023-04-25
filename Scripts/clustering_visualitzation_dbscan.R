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

#Function to get the similarity comparison of molecules on the same cluster
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

df_mols_similarity <- mols_similarity(df_dbscan, cnames, 1)

for (i in 1:5){
  cluster <- filter(df_mols_similarity, cluster == i)
  cat('\nTanimoto coef mean cluster',i,':',sum(cluster['Tanimoto_Coef'])/nrow(cluster))
  cat('\nOverlap coef mean cluster',i,':',sum(cluster['Overlap_Coef'])/nrow(cluster))
}

