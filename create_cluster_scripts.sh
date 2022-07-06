#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH -J create_cluster_scripts
#SBATCH -o create_cluster_scripts_out.out
#SBATCH -e create_cluster_scripts_err.out

module purge
module load system/Python-3.7.4

source env/bin/activate

pip install dexom-python
pip install miom[full]

cd $SLURM_SUBMIT_DIR
python create_cluster_scripts.py
