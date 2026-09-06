#!/bin/bash
#SBATCH --job-name=jupyterlab
#SBATCH --partition=vis,standard,interactive
#SBATCH --cpus-per-task=64
#SBATCH --mem=0
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=amelie.huebner@uni-jena.de
#-- Activate your conda environment which contains the jupyter package
source /home/$USER/.bashrc
conda activate plastchem 
jupyter nbconvert   --to notebook \
                    --execute /work/ki64cah/plast-chem/analysis/binary.ipynb \
                    --output /work/ki64cah/plast-chem/analysis/binary.ipynb \
                    --ExecutePreprocessor.timeout=None
