# Toxicity Prediction of Plastic Chemicals

Project thesis predicting the hazard information of plastic-associated chemicals from molecular structure, using the [PlastChem database](https://doi.org/10.5281/zenodo.10701706) and an XGBoost classifier.

Full methodology, results, and discussion are in [`report/main.pdf`](report/main.pdf).

## Repository structure

```
analysis/
├── utils.py             # Shared helper functions (incl. focal loss)
├── data/                # PlastChem database (not included)
├── binary.ipynb         # pipeline and test figures
├── predict_visual.ipynb # result figures (separate fron gridsearch)
└── visual.ipynb         # Exploratory analysis and figures 
report/                  # figures and pdf of final report
```

## Data

Requires the PlastChem database (version 1.1) TSV file, placed at `data/plastchem_db_v1.1.tsv`.
Not included in this repository, needs to be downloaded from the [PlastChem Zenodo record](https://doi.org/10.5281/zenodo.10701706).

## Environment

Key dependencies:
```
python==3.14.6
numpy==2.5.0
pandas==3.0.3
scipy==1.18.0
scikit-learn==1.9.0
xgboost==3.3.0
imbalanced-learn==0.14.2
rdkit
```
Note: training the grid takes around 2 hours on 32 CPU cores

## Reproducing results

1. Place the PlastChem TSV under `data/` as described.
2. Run `binary.ipynb` for data preprocessing, the grid search, and final test-set evaluation.
3. Run `visual.ipynb` for exploratory figures used in the report.

The final grid search, including the retrained model, can also be loaded directly from `best_model.joblib` / `grid_search.joblib`.
