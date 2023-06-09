########################################
#LOADING R PACKAGES
options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8192m"))
library(rcdk)
library(ChemmineR) 
library(ChemmineOB) 
library(rpubchem)
library(chemblr)
library(Rcpi)
library(cluster)
library("factoextra")
library(fpc)
library(ggcorrplot)
library(plotly)
library(ggplot2)
library(dbscan)
library(matrixStats)
library(RxnSim)
library(hrbrthemes)
library(RColorBrewer)
########################################
#LOAD THE DATAFRAME
#Set the directory 
setwd("/Users/MarionaP/Desktop/TFG/DrugBank")
Df_scaled <- read.csv("df_molecular_descriptors_scaled.csv")
df_original <- read.csv("df_molecular_descriptors.csv")
smiles <- read.csv("df_smiles.csv")
df_filter <- read.csv("df_molecular_descriptors_scaled_filtered.csv")
df_info <- read.csv("df_info_mols.csv")
########################################
#FUNCTIONS
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

#Function to get the similarity comparison of molecules on the same cluster.
#It computes the Tanimoto coeficient and the Overlap Coeficient.
mols_similarity <- function(df_dbscan, cnames, nclusters){
  mol_similarity <- data.frame(mol1 =  character(), mol2 =  character(), Tanimoto_Coef = double(), Overlap_Coef = double(), cluster = integer())
  for (k in 1:nclusters){  
    cluster1 <- filter(df_dbscan, cluster == k)
    c1_id <- cluster1[cnames[cnames %in% c("id")]]
    c1_smiles <- filter(smiles, id %in% as.list(c1_id)$id)
    c_filter <- c1_smiles
    for (i in c1_smiles['smile']$smile){
      c_filter <- filter(c_filter, smile != i)
      for (j in c_filter['smile']$smile){
        sim <- calcDrugMCSSim(i, j, type = 'smile')
        mol_similarity[nrow(mol_similarity) + 1,] <- data.frame(mol1 = i, mol2 = j, Tanimoto_Coef = sim[[2]], Overlap_Coef = sim[[3]], cluster = k)
      }
    }
  }
  for (i in 1:nclusters){
    cluster <- filter(mol_similarity, cluster == i)
    cat('\nTanimoto coef mean cluster',i,':',sum(cluster['Tanimoto_Coef'])/nrow(cluster))
    cat('\nOverlap coef mean cluster',i,':',sum(cluster['Overlap_Coef'])/nrow(cluster))
  }
  return(mol_similarity)
}

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

#Function that calculates the performance of the DBSCAN algorithm for different parameters values based on the Tanimoto Similarity
parameters_dbscan <- function(d, eps, minPts, eps_cl){
  mol_similarity <- data.frame(eps = integer(), minPts = integer(), eps_cl = integer(), nclusters = integer(), non_class = integer(), cluster_max_size = integer(), Similarity = double())
  results <- c()
  for (i in eps){
    for (j in minPts){ 
      optics_res <- optics(d, eps = i, minPts = j)
      for (l in eps_cl){
        total <- 0
        extract_res <- extractDBSCAN(optics_res, eps_cl = l)
        print(extract_res)
        df_hdbscan <- df_original
        df_hdbscan['cluster'] <- extract_res$cluster
        no_class <- nrow(filter(df_hdbscan, cluster == 0))
        nclusters <- length(unique(extract_res$cluster))-1
        if (nclusters > 0){
          max_size <- filter(df_hdbscan, cluster > 0) %>% count(cluster)
          c_max_size <- max(max_size['n'])
          for (k in 1:nclusters){
            print(k)
            cluster <- filter(df_hdbscan, cluster == k)
            if (nrow(cluster) < 1000 && nrow(cluster) > 1){
              c_id <- cluster[cnames[cnames %in% c("id")]]
              c_smiles <- filter(smiles, id %in% as.list(c_id)$id)
              m <- ms.compute.sim.matrix(c_smiles$smile, format = 'SMILES', fp.type='pubchem', sim.method='tanimoto')
              m[lower.tri(m, diag=TRUE)] <- 0
              similarity <- sum(m)/sum(rowCounts(m > 0))
              if(similarity == 'NaN'){total <- total + 0 } #if similarity is 0 the result is NaN
              else{total <- total + similarity}
            }
            else{
              total <- total + 0
            }
          }
          mol_similarity[nrow(mol_similarity) + 1,] <- data.frame(eps = i, minPts = j, minPts = l, nclusters = (length(unique(extract_res$cluster))-1), non_class = no_class, cluster_max_size = c_max_size, Similarity = total/(length(unique(extract_res$cluster))-1))
        }
        else{
          mol_similarity[nrow(mol_similarity) + 1,] <- data.frame(eps = i, minPts = j, minPts = l, nclusters = (length(unique(extract_res$cluster))-1), non_class = no_class, cluster_max_size = 0, Similarity = 0)
        }
      }
    }
  }
  return(mol_similarity)
}

########################################
#Optics
#Prepare the data for clustering
DF <- Df_scaled
cnames <- colnames(DF)
df <- DF[,cnames[!cnames %in% c("X","id")]]

#Perform the clustering
dist_matrix <- daisy(df, metric = "euclidean")
d <- dist(dist_matrix)
optics_res <- optics(d, eps = 2, minPts = 2)
optics_res <- extractDBSCAN(optics_res, eps_cl = 2)

#Add clustering results to a dataframe
df_optics <- df_original
df_optics['cluster'] <- optics_res$cluster

#Get number of clusters
nclusters <- length(unique(df_optics$cluster))-1

#Visualize the clustering
df_clusters <- filter(df_optics, cluster > 0)
ggplot(data = df_clusters, mapping = aes(x = MLogP,y = OB_MW, color = cluster))+geom_point() 

#fig <- plot_ly(df_clusters, x = ~MLogP, y = ~OB_MW, z = ~OB_logP, color = ~cluster)
#fig <- fig %>% add_markers()
#fig <- fig %>% layout(scene = list(xaxis = list(title = 'MLogP'),yaxis = list(title = 'OB_MW'),zaxis = list(title = 'OB_logP')))
#fig

#Visulaize the molecules of each cluster
visualize_cluster_mols(df_optics, cnames, nclusters)

#Compute the similarity between molecules
#Tanimoto and Overlap coeficients (If the size of the clusters is big use the other function mols_similarity_matrix())
df_mols_similarity <- mols_similarity(df_optics, cnames, nclusters)

#Tanimoto coeficient
mols_similarity_matrix(df_optics, cnames, nclusters)

#Get info of the clusters
#SMILES
df_smiles <- smiles
df_smiles['cluster'] <- optics_res$cluster
df_smiles <- filter(df_smiles, cluster > 0)

#OTHER INFO
df_info['cluster'] <- optics_res$cluster
df_info <- filter(df_info, cluster > 0)
for (i in 1:nclusters){
  info <- filter(df_info, cluster == i)
  print(info)
}

#Get clusters that contain approved, vet approved or nutraceutical drugs
interest_cluster <- c()
for (i in 1:nclusters){
  info <- filter(df_info, cluster == i)
  if ('approved' %in% info$charMA.DRUG_GROUPS | 'vet_approved' %in% info$charMA.DRUG_GROUPS | 'nutraceutical' %in% info$charMA.DRUG_GROUPS){
    interest_cluster <- c(interest_cluster, i)
  }
}
df_clust_of_interest <- filter(df_info, cluster %in% interest_cluster)
colnames(df_clust_of_interest)[2] <- "id"
names <- colnames(df_clust_of_interest)
for (i in interest_cluster){
  cluster1 <- filter(df_clust_of_interest, cluster == i)
  c1_id <- cluster1[names[names %in% c("charMA.DRUGBANK_ID")]]
  c1_smiles <- filter(df_smiles, id %in% as.list(c1_id)$charMA.DRUGBANK_ID)
  print(c1_smiles)
  c1_sdf <- c()
  for (i in as.list(c1_smiles['smile'])$smile){
    c1_sdf <- c(c1_sdf, smiles2sdf(i))
  }
  par(mfrow=c(4,4))
  for (i in 1:length(c1_sdf)){
    openBabelPlot(c1_sdf[[i]])
  }
}

df_interest_similarity <- mols_similarity(df_clust_of_interest, names, length(interest_cluster))




#Calculate the performance of the algorithm with different eps and minPts based on the similarity between molecules on the same cluster
dist_matrix <- daisy(df, metric = "euclidean")
d <- dist(dist_matrix)
eps <- seq(2, 5, by=1)
minPts <- seq(2, 10, by=1)
eps_cl <- seq(1, 4, by=0.5)
param_similarity <- parameters_dbscan(d, eps, minPts, eps_cl)
#write.csv(param_similarity, "param.csv")

mycolors <- colorRampPalette(brewer.pal(9, "Blues"))(15)
plot_dbscan <- ggplot(param_similarity, aes(eps, Similarity, colour = factor(minPts))) +
  geom_point()+ scale_fill_manual(mycolors)+geom_line()+theme_ipsum() + ggtitle("Mean Tanimoto Similarity") 
plot_dbscan
plot_dbscan_clusters <- ggplot(param_similarity, aes(eps, nclusters, colour = factor(minPts))) +
  geom_point()+ scale_fill_manual(mycolors)+geom_line()+theme_ipsum() + ggtitle("Number of clusters") 
plot_dbscan_clusters

