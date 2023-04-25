########################################
#LOAD THE DATAFRAME
#Set the directory 
setwd("/Users/MarionaP/Desktop/TFG/DrugBank")
df <- read.csv("df_molecular_descriptors.csv")
########################################
#FILTERING FUNCTION
filter_func <- function(df, MW, Aromatic, MLogP) {
  if (MW == 'S'){
    df <- filter(df, OB_MW < 900)
  }
  else if (MW == 'B'){
    df <- filter(df, OB_MW > 900)
  }
  
  if (Aromatic == TRUE){
    df <- filter(df, naAromAtom > 0)
  }
  else if (Aromatic == FALSE){
    df <- filter(df, naAromAtom == 0)
  }
  
  if (MLogP == 'P'){
    df <- filter(df, MLogP > 0)
  }
  else if (MLogP == 'N'){
    df <- filter(df, MLogP < 0)
  }
  
  write.csv(df, "df_md_filtered.csv")
  return(df)
}

df_f <- filter_func(df, MW = 'S', Aromatic = TRUE, MLogP = FALSE)

filter_func <- function(df, filter) {
  if (filter == 'Rof5'){
    df <- filter(df, OB_HBD <=5, OB_HBA1 <= 10,  OB_MW < 500, OB_logP <= 5)
  }
  write.csv(df, "df_md_filter.csv")
  return(df)
}

df_f <- filter_func(df, 'Rof5')

#Filter function
df_f1 <- filter(prop_df, OB_MW < 500  & Aromatic > 0 & MLogP > 0)


