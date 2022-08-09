#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH -J cluster_install
#SBATCH -o cluster_install_out.out
#SBATCH -e cluster_install_err.out

cd $SLURM_SUBMIT_DIR
git clone https://forgemia.inra.fr/metexplore/cbm/ocmmed.git ocmmed
cd ocmmed

module purge
module load system/Python-3.7.4

python -m venv env 
source env/bin/activate

pip install --upgrade pip
pip install dexom-python
pip install miom[all]

echo "installation complete"
