### Helper functions for the jupyter notebooks
from pathlib import Path
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA


# defining a dictionary for reindexing/later use
col_dict = {
    "Identifiers" : ["plastchem_ID", "cas", "cas_fixed", "identified_by", "pubchem_cid", "pubchem_name", "iupac_name"],

    "Properties" :	["molecular_formula", "molecular_weight", "canonical_smiles", "isomeric_smiles", "inchi", "inchikey", "xlogp", "exact_mass", "monoisotopic_mass", "tpsa", "complexity", "charge"],

    "PlastChem_lists" :	["PlastChem_lists",	"MEA_list",	"MEA_names", "Precedent_list", "Precedent_names", "Red_list", "Orange_list", "Watch_list", "White_list", "Grey_list"],

    "Hazard_information" : ["Hazard_score", "Hazard_score_sum", "Evidence_score", "Persistence_score", "Bioaccumulation_score", "Mobility_score", "Toxicity_score", "Carcinogenic", "Mutagenic", 
                            "Reproduction", "CMR", "STOT", "EDC", "Aquatic_toxicity", "PBT", "vPvB", "PMT", "vPvM"],

    "Groups": ["grouped", "overlap_groups", "priority_groups", "aromatic_amines", "aralkyladehydes", "alkylphenols", "salicylate_esters", "aromatic_ethers", "bisphenols", "orthophthalates", 
            "benzothiazole", "benzotriazoles", "organometallics", "parabens", "azodyes", "acetophenones_benzophenones", "chlorinated_paraffins", "PFASs"],
    
        "Other_groups": ["UVCBs", "polymers", "mixtures", "inorganic_compounds", "DDT_DDE_DDD", "dioxins", "PBDEs", "PCBs", "PBDD_PBDF_PCDD_PCDF",  "aldehydes_simple", "alkanes", "alkenes", "alkynes", 
                        "alkane_ethers", "aliphatic_ketones", "aliphatic_primary_amides", "alkyl_nitrates", "aromatic_hydrocarbons", "carboxylic_acids_salts", "cyclic_acetals", "cyclic_ethers", 
                        "dialiphatic_ethers_excluding_unsatured", "diazo_amino_hydroxyl_naphthalenedisulfonic_acid_dyes", "dibenzoyl_peroxide_derivatives", "dihydropurinediones", "ethanediols", 
                        "isophthalates_terephthalates_trimellitates", "ketones_simple", "organophosphates", "phenolic_antioxidants", "polychlorinated_naphthalenes", "pyrazoles", "salicyclic_acid",
                        "silanes_siloxanes_silicones"],

        "Homologous_series" : ["homo_CCO", "homo_CF2", "homo_CF2CF2O", "homo_CF2O", "homo_CH2"],

        "Elemental_composition" : ["arsenic", "cadmium", "mercury", "chromium", "antimony", "tin", "bromo", "chloro", "fluoro", "iodo", "manganese", "magnesium", "barium", "nickel", "lead", "vanadium",
                                    "tellurium", "thallium", "beryllium", "selenium", "ECHA_grouping"], 
    
    "Global_regulation" : ["MEA_list", "MEA_names", "Minamata", "Stockholm", "Montreal", "Basel"], 

    "Regional_national_regulation" : ["Precedent_list", "Precedent_names", "US_summary", "EU_summary", "Rotterdam_summary", "Japan_summary", "Korea_summary", "California_P65", "EU_toy_restricted",
                                    "EU_toy_allergenic_fragrances", "REACH_authorisation", "REACH_restriction", "REACH_SVHC", "ROHS_directive", "Rotterdam", "Rotterdam_PIC", "CSCL_Class1", 
                                    "CSCL_Class2", "ISHA_permission_required", "ISHA_prohibited", "Korea_accidents", "Korea_CMR", "Korea_hazardous_chemicals", "Korea_intensive_control"],

    "Use, presence, release (UPR)": ["UPR_summary", "Use", "Presence", "Release"],

        "Commodity_plastics": ["Commodity_plastics_summary", "PE_combined", "PE", "HDPE", "LDPE", "PP", "PS_combined", "PS", "EPS", "HIPS", "PVC", "PET", "PA", "PUR"],

        "Specialty_plastic": ["Specialty_plastics_summary", "ABS", "BPC", "EVA", "EVOH", "Melamine", "PAA", "PAN", "PC", "PCT", "PES", "PMMA", "PVA", "SAN", "SAP", "Tritan"],

        "Elastomers": ["Elastomers_summary", "Rubber", "Silicone"],

        "Bioplastics": ["Bioplastics_summary", "Cellulose-based", "PBAT", "PBS", "PHA", "PLA", "Starch", "Bioplastics_unspecified_1", "Bioplastics_unspecified_2"],

        "Unspecified_plastics": ["Unspecified_plastics_summary", "Multilayer", "Unspecified_plastics"],

    "Production Volume":[ "total_production_volume_tons", "EU_production_volume_tons", "nordic_countries_spin_production_volume_tons", "usa_production_volume_tons", "oecd_production_volume_tons"],

    "Functions":["Harmonized_functions", "original_function_plasticmap", "original_function_cpp", "original_primary_function_aurisano", "original_other_function_aurisano", "industrial_sector_plasticmap",
                "Industrial_sector_food_contact_plasticmap"],

    "Information_sources": ["in_aurisano", "in_cpp", "in_echa", "in_fcc", "in_fccmigex", "in_litchem", "in_plasticmap", "original_name_aurisano", "original_name_cpp", "original_name_echa", 
                            "original_name_fcc", "original_name_fccmigex", "original_name_litchem", "original_name_plasticmap"]
}

# import data
def load_plastchem(path: Path) -> pd.DataFrame:
    """
    Reads in the PlastChem database. 

    Parameters
    ----------
    path: Path
        System path to the database.

    Returns
    -------
    pd.Dataframe
        The PlastChem Database, multindexed using the cols_dict.
    """
    plast_chem = pd.read_csv(path, sep="\t", decimal=',', low_memory=False)
    # remove empty columns and rows
    plast_chem.dropna(how='all', axis=1, inplace=True)
    plast_chem.dropna(how='all', axis=0, inplace=True)

    # reindexing and multiindexing, then combining multiindex into one → easier to work with (naming convention is higherlevel_lower)
    tuples = [(k, i) for k, v in col_dict.items() for i in v]
    index = pd.MultiIndex.from_tuples(tuples)
    df_plast_chem = plast_chem[1:]
    df_plast_chem.columns = ["_".join(col) for col in index.to_flat_index()]
    
    # drop first row → now obsolete
    df_plast_chem.dropna(axis=0, subset=['Identifiers_plastchem_ID'])
    # replace all commas with dots
    pcc = df_plast_chem.columns
    df_plast_chem = df_plast_chem[2:].stack().str.replace(',','.').unstack()
    return df_plast_chem, pcc


def get_MI_tupels(data, vars: list):
    cols_to_use = []
    for top_level in vars:
        # Get all columns under the given top-level category from multiindex
        cols_to_use.extend([col for col in data.columns if col[0] == top_level])
    return cols_to_use

def extract_feature_vecs(data: pd.DataFrame, vars: list, n):
    """
    Extract feature vectors from data using scikits implementation of PCA. 

    Parameters
    ----------
    data: pd.DataFrame
        PlastChem database. Multiindexed. Shape (n_chemicals, n_features).
    vars: list
        List of tuples of variable names for which to extract features.
    """
    # turn pandas df into numpy for PCA feature vector extraction
    assert all(var in data.columns for var in vars), "Target variables not found in dataset"
    subset = data[vars]
    cleaned_subset = subset.dropna(axis=0)
    X = cleaned_subset.to_numpy()
    feature_vectors = PCA(n_components=n).fit_transform(X)
    return feature_vectors


