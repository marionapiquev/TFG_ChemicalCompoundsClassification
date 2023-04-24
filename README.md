# TFG_ChemicalCompoundsClassification

## Project Files Description

This project includes two directories, one for Scripts and another for Data.

#### Scripts:
- `get_molecular_descriptors.R`: Includes all functions required to extract molecular descriptors.
- `data_preprocessing.R`: Source file with data transformation process to prepare data for clustering.
- `clustering_molecular_descriptors.R`: 

#### Data:
- `df_molecular_descriptors.csv`: Dataframe obtained using *extract_molecular_descriptors* function from `get_molecular_descriptors.R` file. It contains the molecular descriptors of 9213 chemical compounds.
- `df_molecular_descriptors_scaled.csv`: Dataframe obtained executing the `data_preprocessing.R` file. It contains the molecular descriptors of 9213 chemical compounds after data preprocessing.
- `df_smiles.csv`: Dataframe obtained using *extract_smiles* function from `get_molecular_descriptors.R` file. It contains the smiles of 9213 chemical compounds.
