#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH -J cluster_install
#SBATCH -o cluster_install_out.out
#SBATCH -e cluster_install_err.out
#SBATCH -c 8

cd $SLURM_SUBMIT_DIR
git clone https://forgemia.inra.fr/metexplore/cbm/ocmmed.git ocmmed_mouse
cd ocmmed_mouse

echo "creating python environment'"

module purge
module load system/Python-3.7.4

pip install --upgrade pip

python -m venv env
source env/bin/activate

echo "installing packages"

pip install dexom-python
pip install miom[all]
pip install snakemake

echo "installation complete"
