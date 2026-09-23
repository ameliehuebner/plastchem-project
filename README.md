# Toxicity Prediction of Plastic Chemicals

Project thesis predicting the hazard information of plastic-associated chemicals from molecular structure, using the [PlastChem database](https://doi.org/10.5281/zenodo.10701706) and a XGBoost classifier.

Full methodology, results, and discussion are in [`report/main.pdf`](report/main.pdf).

## Repository structure

```
analysis/
├── data/                    
│   ├── best_model.joblib     # Retrained model
│   ├── plastchem_db_v1.1.tsv # PlastChem database (not included)
│   └── ...
├── binary.ipynb              # Pipeline for model training
├── environment.yml           # CONDA environment file
├── utils.py                  # Shared helper functions (incl. focal loss)
├── visual_dataset.ipynb      # Exploratory analysis and figures
└── visual_model.ipynb        # Result figures
report/                       # Figures and PDF of final report
```

## Data

Requires the PlastChem database (version 1.1) TSV file, placed at `data/plastchem_db_v1.1.tsv`.
Not included in this repository, needs to be downloaded from the [PlastChem Zenodo record](https://doi.org/10.5281/zenodo.10701706).

## Environment

Key dependencies:
```
python==3.14.6
pandas==3.0.3
scikit-learn==1.9.0
xgboost==3.3.0
imbalanced-learn==0.14.2
rdkit=2026.03.5
```
Build from environment.yml:
```
conda env create -f environment.yml
```

## Reproducing results

1. Place the PlastChem TSV under `data/` as described.
2. Setup the environment using the `environment.yml`.
3. Run `binary.ipynb` for data preprocessing, the grid search, and final test-set evaluation.
4. Run `visual.ipynb` for exploratory figures used in the report.

The final grid search, including the retrained model, can also be loaded directly from `best_model.joblib` / `grid_search.joblib`.
