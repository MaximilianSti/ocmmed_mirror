#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH -J cluster_install
#SBATCH -o cluster_install_out.out
#SBATCH -e cluster_install_err.out
#SBATCH -c 8

cd $SLURM_SUBMIT_DIR
git clone https://forge.inrae.fr/metexplore/cbm/ocmmed.git ocmmed
cd ocmmed

echo "creating python environment"

module purge
module load devel/python/Python-3.7.9

python -m venv env
source env/bin/activate

pip install --upgrade pip

echo "installing packages"

pip install dexom-python
pip install snakemake
pip install pulp==2.7.0
pip  uninstall --yes gurobipy

echo "downloading test model"

python utilities_cluster/cluster_installation_helper.py

echo "installation complete"
