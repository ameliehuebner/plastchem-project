#IMPORT ALL REQUIRED LIBRARIES
from rdkit import Chem
from sklearn.manifold import TSNE
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit import DataStructs
from sklearn.cluster import KMeans
import pubchempy as pcp
from sklearn.model_selection import train_test_split
from sklearn.naive_bayes import GaussianNB
from scipy.stats import pearsonr
from sklearn.metrics import accuracy_score, confusion_matrix
from rdkit.Chem import rdDepictor
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D
import urllib.request
import random
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn import svm
from sklearn.metrics import accuracy_score
import os

##STAGE 1: Generate ACTIVE SMILES FILE

# load SMILES strings from a text file
cwd = os.getcwd()
parent_dir = os.path.abspath(os.path.join(cwd, os.pardir))
chemical_file=parent_dir+'\\0_skdata\chemical_names.txt'
with open(chemical_file, 'r') as f:
    names_list = f.readlines()

#names_list=['Citronellal','Geranial'] 
# search PubChem for each SMILES string and retrieve the corresponding PubChem ID
results = []
for name in names_list:
    name = name.strip()
    search_results = pcp.get_cids(name, 'name')
    if search_results:
        pubchem_id = search_results[0]
        compound = pcp.Compound.from_cid(pubchem_id)
        smiles = compound.canonical_smiles
        row = {'Chemical_Name': name, 'SMILES': smiles, 'PubChem_ID': pubchem_id,'activity': 1,'state':'ACTIVE'}
        results.append(row)

# save results to a CSV file
active_smiles_df = pd.DataFrame(results)
active_smiles_df.to_csv(parent_dir+'\\0_skdata\chemical_names_pubchem_active.csv', index=False)


#data=active_smiles_df+inactive_smiles_df
cwd = os.getcwd()
parent_dir = os.path.abspath(os.path.join(cwd, os.pardir))
active_data_df=pd.read_csv(parent_dir+'\\0_skdata\chemical_names_pubchem_active.csv')
active_data_df['PubChem_ID'] = active_data_df['PubChem_ID'].astype(str)
inactive_data_df=pd.read_csv(parent_dir+'\\0_skdata\inactive_dude_db_top4.csv')
data = pd.concat([active_data_df, inactive_data_df])
data=data.reset_index(drop=True)
data['data_index'] = data.index




#selected_cols=['activity', 'SMILES','Chemical_Name', 'PubChem_ID']
#data = pd.concat([active_smiles_df[selected_cols], inactive_smiles_df[selected_cols]])
#data
nBitsK=4096

data_fp=[]
active_molecules = [AllChem.MolFromSmiles(smiles) for smiles in data[data['activity'] == 1]['SMILES']]
active_fingerprints = [AllChem.GetMorganFingerprintAsBitVect(mol, 6, nBits=nBitsK) for mol in active_molecules]
active_fingerprints_array = []
for fp in active_fingerprints:  
    arr = np.zeros((1,))
    AllChem.DataStructs.ConvertToNumpyArray(fp, arr)
    active_fingerprints_array.append(arr)

inactive_molecules = [AllChem.MolFromSmiles(smiles) for smiles in data[data['activity'] == 0]['SMILES']]
inactive_fingerprints = [AllChem.GetMorganFingerprintAsBitVect(mol, 4, nBits=nBitsK) for mol in inactive_molecules]
inactive_fingerprints_array = []
for fp in inactive_fingerprints:
    arr = np.zeros((1,))
    AllChem.DataStructs.ConvertToNumpyArray(fp, arr)
    inactive_fingerprints_array.append(arr)

# Combine active and inactive fingerprints
fingerprints_array = active_fingerprints_array + inactive_fingerprints_array
data_fp = pd.DataFrame(fingerprints_array)

# Perform t-SNE to reduce dimensionality
tsne = TSNE(n_components=2, random_state=1)
reduced_data = tsne.fit_transform(data_fp)

# Perform K-means clustering
kmeans = KMeans(n_clusters=2)
kmeans.fit(reduced_data)

# Plot scatter plot of reduced data with cluster labels
plt.scatter(reduced_data[:len(active_fingerprints_array), 0], reduced_data[:len(active_fingerprints_array), 1],color='purple', label='Active')
plt.scatter(reduced_data[len(active_fingerprints_array):, 0], reduced_data[len(active_fingerprints_array):, 1],color='pink', label='Inactive')
plt.scatter(kmeans.cluster_centers_[:, 0], kmeans.cluster_centers_[:, 1], color='black', marker='x', s=200)
plt.xlabel('t-SNE1')
plt.ylabel('t-SNE2')
plt.title('t-SNE plot of Morgan fingerprints for Active and Inactive compounds')
plt.legend()
plt.show()




#data=active_smiles_df+inactive_smiles_df
cwd = os.getcwd()
parent_dir = os.path.abspath(os.path.join(cwd, os.pardir))
active_data_df=pd.read_csv(parent_dir+'\\0_skdata\chemical_names_pubchem_active.csv')
active_data_df['PubChem_ID'] = active_data_df['PubChem_ID'].astype(str)
inactive_data_df=pd.read_csv(parent_dir+'\\0_skdata\inactive_dude_db_top4.csv')
data = pd.concat([active_data_df, inactive_data_df])
data=data.reset_index(drop=True)
data['data_index'] = data.index




##STAGE 5.1: GENERATE DESCRIPTORS, FIND DESCRIPTORS WITH CORRELATION COFFIENT MORE THAN 0.5 AND APPLY ML MODEL

# Convert SMILES strings to molecules
mols = [Chem.MolFromSmiles(smi) for smi in data['SMILES']]
data['Molecule'] = mols

# Calculate molecular descriptors
#X = pd.DataFrame()
data['MW'] = [Descriptors.MolWt(mol) for mol in mols]
data['LogP'] = [Descriptors.MolLogP(mol) for mol in mols]
data['NumHAcceptors'] = [Descriptors.NumHAcceptors(mol) for mol in mols]
data['NumHDonors'] = [Descriptors.NumHDonors(mol) for mol in mols]
data['NumRotatableBonds'] = [Descriptors.NumRotatableBonds(mol) for mol in mols]

data['TPSA'] = [Descriptors.TPSA(mol) for mol in mols]
data['NumAromaticRings'] = [Descriptors.NumAromaticRings(mol) for mol in mols]
data['NumAliphaticRings'] = [Descriptors.NumAliphaticRings(mol) for mol in mols]
data['NumSaturatedRings'] = [Descriptors.NumSaturatedRings(mol) for mol in mols]
data['NumHeteroatoms'] = [Descriptors.NumHeteroatoms(mol) for mol in mols]
data['NumHeavyAtoms'] = [Descriptors.HeavyAtomCount(mol) for mol in mols]
data['NumRings'] = [Descriptors.RingCount(mol) for mol in mols]
data['NumValenceElectrons'] = [Descriptors.NumValenceElectrons(mol) for mol in mols]
data['NumRadicalElectrons'] = [Descriptors.NumRadicalElectrons(mol) for mol in mols]
data['ExactMolWt'] = [Descriptors.ExactMolWt(mol) for mol in mols]
data['FractionCSP3'] = [Descriptors.FractionCSP3(mol) for mol in mols]
X_train=[]
X_test=[]
y_test=[]
y_train=[]
X=[]
y=[]


# Set activity as target variable
y = data['activity']
#X=data
X =  data[['MW', 'LogP', 'NumHDonors', 'NumHAcceptors', 'NumRotatableBonds','TPSA','NumAromaticRings','NumAliphaticRings','NumSaturatedRings','NumHeteroatoms','NumHeavyAtoms','NumRings','NumValenceElectrons','NumRadicalElectrons','ExactMolWt','FractionCSP3', 'data_index']]
# Split data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=testsize, stratify=y, random_state=randomstate)
#print(X.shape, X_train.shape, X_test.shape)
selected_features=['MW', 'LogP', 'NumHDonors', 'NumHAcceptors', 'NumRotatableBonds','TPSA','NumAromaticRings','NumAliphaticRings','NumSaturatedRings','NumHeteroatoms','NumHeavyAtoms','NumRings','NumValenceElectrons','NumRadicalElectrons','ExactMolWt','FractionCSP3']


# Calculate Pearson correlation coefficient for each feature
corr_dict=dict()
for column in X_train[selected_features].columns:
    corr, _ = pearsonr(X_train[column], y_train)
    corr_dict[column] = corr

# Select highly correlated features
highly_corr_features = [feat for feat, corr in corr_dict.items() if abs(corr) > corr_factor]



##STAGE 5.2: Apply Machine Learning Model and Predict Probability
nb = GaussianNB()
nb.fit(X_train[highly_corr_features], y_train)
probabilities =       nb.predict_proba(X_test[highly_corr_features])[:, 1]
probabilities_train = nb.predict_proba(X_train[highly_corr_features])[:, 1]

y_pred = (probabilities > probability_factor).astype(int)

accuracy = accuracy_score( y_pred,y_test)

#create an actual column from pseduo column of index which can be used for joining two dataframe

X_test['ProbActivity'] = probabilities
X_test['test_type'] = 'Test Set'

X_train['ProbActivity'] = probabilities_train
X_train['test_type'] = 'Training Set'
print(X.shape, X_train.shape, X_test.shape)

custom_features=['data_index','ProbActivity','test_type']

#custom_features.extend(highly_corr_features)
processed_data=pd.concat([X_test[custom_features], X_train[custom_features]])




# Paper title:Projecting phytochemical bacoside A anti-mucorale agent: an in-silico and in-vitro assessment
# Authors: Komal Tilwani1,  Drashti Patel1, Prachi Soni1, Dr Gayatri Dave 1*
# Author Affiliations: 
# 1.P D Patel Institute of Applied Sciences, CHARUSAT, Changa- 388421, Anand Gujarat, India
# Corresponding Author:
# Dr Gayatri Dave, PhD.
# Associate Professor,
# Department of Biological Science,
# P D Patel Institute of Applied Sciences, CHARUSAT
# E-mail: gayatridave.bt@charusat.ac.in
# Phone: (+91)-7600414303

# ORCID ID:
# Dr Gayatri Dave- https://orcid.org/0000-0002-8510-0179
# Komal Tilwani -  https://orcid.org/0000-0001-7433-2831

