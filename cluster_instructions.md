
# OCMMED cluster pipeline

| Table of contents                                                                          |
|--------------------------------------------------------------------------------------------|
| [Installation instructions](cluster_instructions.md#Installation-instructions)             |
| [Installing cplex on the cluster](cluster_instructions.md#Installing-cplex-on-the-cluster) |
| [Setting parameters](cluster_instructions.md#Setting-parameters)                           |
| [Running the pipeline](cluster_instructions.md#Running-the-pipeline)                       | 

## Installation instructions

The `cluster_install.sh` file contains the commands for downloading and installing OCMMED on a slurm cluster.  
Copy this file to your working directory on the slurm cluster, then run the command `sbatch cluster_install.sh`.    
Below is an explanation of the commands that are used in this script.

```
#!/bin/bash  
#SBATCH --mail-type=ALL  
#SBATCH -J cluster_install 
#SBATCH -o cluster_install_out.out  
#SBATCH -e cluster_install_err.out  
```
These commands indicate that this is a bash script, and when it is sent to the slurm cluster via the sbatch command, it will create a job with the name "cluster_install" which will write its output to the "cluster_install_out.out" file and any potential errors or warnings into the "cluster_install_err.out" file. The cluster will also send you emails anytime the run starts, ends, or fails.

```
cd $SLURM_SUBMIT_DIR  
git clone https://forgemia.inra.fr/metexplore/cbm/ocmmed ocmmed 
cd ocmmed
```
With these commands, the script will copy the ocmmed repository into the current directory, in a folder called "ocmmed".

```
module purge  
module load system/Python-3.7.4
```
`module purge` deactivates any currently loaded modules.  
Here, we load a module called `Python-3.7.4` which is installed in the `system` folder of the slurm cluster.   
If your python distribution is installed elsewhere, or if you are using a different version of python, you must modify this line. (In this case, you must also modify the same line in the `create_cluster_scripts.sh` file, as well as the `pythonpath` parameter in the `cluster_params.yaml` file)

```
python -m venv env  
source env/bin/activate
```
Here, we create a python virtual environment with the name "env", and we activate this environment (assuming that the cluster is running linux).

```
pip install --update pip  
pip install dexom-python  
pip install miom[all]
```
These commands upgrade pip and install the necessary packages for the OCMMED pipeline.

## Installing cplex on the cluster
For installing cplex, you will need to create an IBM account with an academic email address.  
https://myibm.ibm.com/dashboard/

Once logged in, you can search for "cplex Optimization Studio"  
https://www.ibm.com/products/ilog-cplex-optimization-studio

Click on “Get no-cost academic edition”, and give your academic e-mail address for confirmation.  
https://www.ibm.com/academic/technology/data-science

Then, in the header “Software”, click on the download button for “ILOG CPLEX Optimization Studio”. This will take you to the IBM software download page.  

The website will automatically direct you to the latest version of the program. If you want to install an older version for backwards compatibility, click on “Search for Software” and enter “cplex optimization studio” to choose the version you want.

You will be taken to a list of downloadable files for different OS. For installing on the cluster, you will need the Linux X86-64 distribution.

Check the box(es) next to the download(s) you want, scroll down to the bottom of the page, check "I agree" and click on "Download Now". Your download will open with IBM Download Director (it may ask for permission: select "Yes").

After the download, you will have a .bin file, which you can transfer to the cluster with WinSCP

To run the file, execute the following commands:  
```
chmod +x cplex_studio1210.linux-x86-64.bin

./cplex_studio1210.linux-x86-64.bin
```
During the installation, you will be asked to specify a path. You must select a path that you have permission for (for example in my case `/home/username/save/CPLEX`)


## Setting parameters
- Modify the necessary parameters in `parameters.yaml`.

In particular, make sure that the `output_path` parameter is an existing folder in the ocmmed diretory (Or leave empty if you want to save results directly in the ocmmed folder).  
If the folder doesn't exist, you will get an error when the program tries to save its results.  

- Modifying the parameters in `additional_params.yaml` is optional. 

If your iterations take very long, you may want to set a higher `timelimit` or increase one of the tolerance parameters.

- Modify the parameters in `cluster_params.yaml`

If, during the installation, you modified the path to the python module `system/Python-3.7.4`, then you must also modify the `pythonpath` parameter  
The `cluster_files` parameter must also be an existing folder in the ocmmed directory.  
It's recommended to use a separate folder for the cluster files, as the program will create many files which can be deleted after the run is finished.

Note that on the cluster, the `rxn_enum_iterations` and `div_enum_iterations` parameter in the `parameters.yaml` file are ignored.  
Instead, you must use the `batch_num`, `batch_rxn_sols` and `batch_div_sols` parameters in the `cluster_params.yaml` file

The `separate` approach is divided into 4 different scripts: 
- one for computing reaction-weights and a starting iMAT solution
- one for launching batches of reaction-enumeration
- one for launching batches of diversity-enumeration
- one for producing the final output

The `grouped` approach is condensed into 3 scripts (reaction-enumeration and diversity-enumeration are joined). There is a slightly higher probability of finding duplicate solutions, but for large solution numbers this effect is negligible.

## Running the pipeline

Once you have changed all the necessary parameters in the different .yaml files, you can run create_cluster_scripts:  
```
sbatch create_cluster_scripts.sh
```
Reminder: If, during the installation, you modified the path to the python module `system/Python-3.7.4`, then you must also do it here

If you are using the `grouped` approach, this will create 3 files:  
```
cluster_script_1.sh
cluster_script_2.sh
cluster_script_3.sh
```
If you are using the `separate` approach, it will also create `cluster_script_4.sh`

You can then launch the first script: `sbatch cluster_script_1.sh`

You must then wait until all the computation is finished before launching the next script.  
The script is finished once your output-path contains reaction-weights and imat-solutions for every experimental condition.

You can then launch the second script: `sbatch cluster_script_2.sh`

This script will take longer before it is finished.  
In the separate approach, you must wait until your cluster-files path contains a file called `rxn_enum_solutions_conditionname_batchnumber.csv` for every condition & every batch.  
In the grouped appraoch, you must additionally wait until you have every `div_enum_solutions_conditionname_batchnumber.csv`

You can then launch the second script: `sbatch cluster_script_3.sh`
In the separate approach, this will now create the `div_enum_solutions_conditionname_batchnumber.csv` files  
In the grouped approach, this will create the `cellspecific_model.xml` file in your output-path

In the separate approach, `sbatch cluster_script_4.sh` will create the `cellspecific_model.xml` file in your output-path.