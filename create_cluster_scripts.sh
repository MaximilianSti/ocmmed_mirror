#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH -J create_cluster_scripts
#SBATCH -o create_cluster_scripts_out.out
#SBATCH -e create_cluster_scripts_err.out

cd $SLURM_SUBMIT_DIR
module purge
module load system/Python-3.7.4

source env/bin/activate

python create_cluster_scripts.py
