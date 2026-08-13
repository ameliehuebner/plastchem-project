#!/bin/bash
#SBATCH --job-name=jupyterlab
#SBATCH --partition=vis,standard,interactive
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#-- Activate your conda environment which contains the jupyter package
source /home/$USER/.bashrc
conda activate plastchem        # NOTE: you have to create an environment with the `jupyter` package
#-- Jupyter is served at a dedicated port. To not conflict with other users we choose a random free port
read LOWERPORT UPPERPORT < /proc/sys/net/ipv4/ip_local_port_range
while :
do
   PORT="`shuf -i $LOWERPORT-$UPPERPORT -n 1`"
   ss -lpn | grep -q ":$PORT " || break
done
#-- Print what we need to copy/past to connect
echo "INFO: Create ssh tunnel from your pc using: ssh -L $PORT:$HOSTNAME:$PORT  $USER@login1.draco.uni-jena.de"
echo "INFO: You can also use draco instead of ara"
#-- Launch jupyterlab ------------------------
jupyter lab --no-browser --ip=0.0.0.0 --port=${PORT}
