#!/bin/bash
#SBATCH --mail-type=NONE
#SBATCH -J imat_DMSO
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -t 00:10:00
#SBATCH -o pilot_data_cluster_files/imat_DMSO_out.out
#SBATCH -e pilot_data_cluster_files/imat_DMSO_err.out
cd $SLURM_SUBMIT_DIR
module purge
module load system/Python-3.7.4
source env/bin/activate
export PYTHONPATH=${PYTHONPATH}:"/home/mstingl/save/CPLEX_Studio1210/cplex/python/3.7/x86-64_linux"
python utilities_cluster/cluster_weights_imat.py -c DMSO