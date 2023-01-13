
# OCMMED cluster pipeline

| Table of contents                                                                                                    |
|----------------------------------------------------------------------------------------------------------------------|
| [Installation instructions](cluster_instructions.md#Installation-instructions)                                       |
| [Installing cplex on the cluster](cluster_instructions.md#Installing-cplex-on-the-cluster)                           |
| [Setting parameters](cluster_instructions.md#Setting-parameters)                                                     |
| [Running the pipeline with snakemake](cluster_instructions.md#Running-the-pipeline-with-snakemake)                   |
| [Running the pipeline with slurm (deprecated)](cluster_instructions.md#Running-the-pipeline-with-slurm-(deprecated)) | 
| [Files created by the cluster scripts](cluster_instructions.md#Files-created-by-the-cluster-scripts)                 | 

## Installation instructions

The `cluster_install.sh` file contains the commands for downloading and installing OCMMED on a slurm cluster.  
Copy this file to your working directory on the slurm cluster, then run the command `sbatch cluster_install.sh`   
If you get the following error: `sbatch: error: Batch script contains DOS line breaks (\r\n) instead of expected UNIX line breaks (\n).`, you must use the command `dos2unix cluster_install.sh` before running sbatch.

Below is an explanation of the commands that are used in the `cluster_install.sh` script.

```
#!/bin/bash  
#SBATCH --mail-type=ALL  
#SBATCH -J cluster_install 
#SBATCH -o cluster_install_out.out  
#SBATCH -e cluster_install_err.out  
```
These commands indicate that this is a bash script, and when it is sent to the slurm cluster via the sbatch command, it will create a job with the name "cluster_install" which will write its output to the "cluster_install_out.out" file and any potential errors or warnings into the "cluster_install_err.out" file. The cluster will also send you an email anytime the run starts, ends, or fails.

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
If your python distribution is installed elsewhere, or if you are using a different version of python, you must modify this line. (In this case, you must also modify the same line in the `create_cluster_scripts.sh` file, as well as the `pythonmodule` parameter in the `cluster_params.yaml` file)

```
python -m venv env  
source env/bin/activate
```
Here, we create a python virtual environment with the name "env", and we activate this environment (assuming that the cluster is running linux).

```
pip install --update pip  
pip install dexom-python  
pip install miom[full]
pip install snakemake
```
These commands update the package installer for python and then install the necessary packages for the OCMMED pipeline.

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

After the download, you will have a .bin file, which you can transfer to the cluster.

To run the file, execute the following commands (this example uses CPLEX version 12.10):  
```
chmod +x cplex_studio1210.linux-x86-64.bin
./cplex_studio1210.linux-x86-64.bin
```
During the installation, you will be asked to specify a path. You must select a path in which you have permission and enough space to install the software (for example in my case `/home/username/save/CPLEX`)


## Setting parameters
- Modify the necessary parameters in `parameters.yaml`.

- Modifying the parameters in `additional_params.yaml` is optional. 

If your iterations take very long and often hit the timelimit, you may want to set a higher `timelimit` or increase one of the tolerance parameters.

- Modify the parameters in `cluster_params.yaml`

If, during the installation, you modified the path to the python module `system/Python-3.7.4`, then you must also modify the `pythonmodule` parameter  
For the `cluster_files` parameter, I recommended using a different folder than `output_path`, as the program will create many files which are not necessary to keep after the run is finished.

Note that on the cluster, the `rxn_enum_iterations` and `div_enum_iterations` parameters in the `parameters.yaml` file are ignored.  
Instead, you must use the `batch_num`, `batch_rxn_sols` and `batch_div_sols` parameters in the `cluster_params.yaml` file to set the number of iterations.

## Running the pipeline with snakemake

For better reproducibility and scalability, we now use the [snakemake](https://snakemake.github.io/) workflow managment tool to execute the OCMMED pipeline on computer clusters.  
Once you have installed the necessary packages and defined the parameters in the .yaml files, you can launch the snakemake workflow which is detailed in the `Snakefile`.  

If your cluster uses **slurm**, you can launch the entire pipeline simply via `sbatch submit_slurm.sh`. Once again, make sure that the file contains the correct paths to your cluster's python module, environment, and CPLEX installation.  
The first snakemake command in `submit_slurm.sh` will create `dag.pdf`, which contains a graph of all the snakemake rules that will be executed during the pipeline.  You can use this graph to check that the pipeline is being executed on the correct conditions and with the correct number of parallel batches.  
The second snakemake command uses the `submit_slurm.py` script as a wrapper for running the workflow. Note that the parameter  `-j 500` signifies that snakemake will submit up to 500 jobs at once; if you require more parallel batches, you may need increase this parameter.

For running the pipeline on other types of clusters, please refer to the [snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cluster.html).  
Configurable profiles for various clusters can be found here: [https://github.com/Snakemake-Profiles](https://github.com/Snakemake-Profiles).

## Running the pipeline with slurm (deprecated)

**This way of running the OCMMED pipeline on slurm cluster is deprecated. We recommend using snakemake instead.**

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

The `separate` approach is divided into 4 different scripts: 
- one for computing reaction-weights and a starting iMAT solution
- one for launching batches of reaction-enumeration
- one for launching batches of diversity-enumeration
- one for producing the final output

The `grouped` approach is condensed into 3 scripts (reaction-enumeration and diversity-enumeration are joined). There is a slightly higher probability of finding duplicate solutions, but for large solution numbers this effect is negligible.

Launch the first script: `sbatch cluster_script_1.sh`  
You must then wait until all the computation is finished before launching the next script.  
On a slurm cluster, you can use the command `squeue -u username` to check if any of your scripts are still running.  
The script is finished once your output-path contains reaction-weights and imat-solutions for every experimental condition.

You can then launch the second script: `sbatch cluster_script_2.sh`  
This script will take longer before it is finished.  
In the _separate_ approach, you must wait until your cluster-files path contains a file called `rxn_enum_solutions_conditionname_batchnumber.csv` for every condition & every batch.  
In the _grouped_ appraoch, you must additionally wait until every `div_enum_solutions_conditionname_batchnumber.csv` has been created.

You can then launch the next script: `sbatch cluster_script_3.sh`  
In the _separate_ approach, this will now create the `div_enum_solutions_conditionname_batchnumber.csv` files  
In the _grouped_ approach, this will create the `cellspecific_model.xml` file in your output-path

In the _separate_ approach, `sbatch cluster_script_4.sh` will create the `cellspecific_model.xml` file in your output-path.

## Files created by the cluster scripts

Below is an overview of all the files that are produced by the cluster scripts.   
This can be useful to check if a script has correctly finished running before launching the next script.

**Note: when using snakemake, the pipeline is no longer divided into several scripts, but it still produces the same output files in the same order**

**cluster_script_1**  
If the parameter `qualitative` is set to `true`, this will produce the file geneweights_qualitative.csv in the output_path.  
Then, for each condition, this script produces the following files:  
- **in output_path:**  
  - reaction_weights_condition.csv 
  - imat_solution_condition.csv
- **in cluster_files:**  
  - imat_condition_err.out  
  - imat_condition_out.out  

Therefore, the total number of newly created files is: num_conditions * 4  (+ 1 if geneweights_qualtitative is present)  

### approach: separate

**cluster_script_2**  
This script only produces files in the cluster_files folder.  
- **for every condition:**  
    - rxnenum_conditionnumber_err.out  
    - rxnenum_conditionnumber_out.out  
- **for every condition and every batch:**  
    - rxnenum_conditionnumber_batchnumber_err.out  
    - rxnenum_conditionnumber_batchnumber_out.out  
    - rxn_enum_solutions_condition_batchnumber.csv  

The total number of newly created files is: num_conditions * (num_batches * 3 + 2 )

**cluster_script_3**  
This script only produces files in the cluster_files folder.  
- **for every condition:**  
    - full_rxn_enum_solutions_condition.csv  
    - divenum_conditionnumber_err.out  
    - divenum_conditionnumber_out.out  
- **for every condition and every batch:**  
    - divenum_conditionnumber_batchnumber_err.out  
    - divenum_conditionnumber_batchnumber_out.out  
    - div_enum_condition_batchnumber_solutions.csv  
    - div_enum_condition_batchnumber_stats.csv  

The total number of newly created files is: num_conditions * (num_batches * 4 + 3)  

**cluster_script_4**  
- **in cluster_files, for every condition:**  
    - full_div_enum_solutions_condition.csv  
- **in output_path, for every condition:**  
    - all_DEXOM_solutions_condition.csv  
- **in output_path:**  
    - all_DEXOM_solutions.csv  
    - activation_frequency_reactions.csv  
    - cellspecific_model.xml
    - inactive_pathways.csv
    - inactive_pathways_relative.csv

The total number of newly created files is: num_conditions * 2 + 5  

### approach: grouped
**cluster_script_2**  
This script only produces files in the cluster_files folder.  
- **for every condition:**  
	- dexom_conditionnumber_err.out  
	- dexom_conditionnumber_out.out  
- **for every condition and every batch:**  
    - dex_conditionnumber_batchnumber_err.out  
    - dex_conditionnumber_batchnumber_out.out  
    - rxn_enum_solutions_condition_batchnumber.csv  
    - div_enum_condition_batchnumber_solutions.csv  
    - div_enum_condition_batchnumber_stats.csv  

The total number of newly created files is: num_conditions * (num_batches * 5 + 2)  

**cluster_script_3**  
- **in cluster_files, for every condition:**  
    - full_rxn_enum_solutions_condition.csv  
    - full_div_enum_solutions_condition.csv  
- **in output_path, for every condition:**  
    - all_dexom_solutions_condition.csv  
- **in output_path:**  
    - all_DEXOM_solutions.csv  
    - activation_frequency_reactions.csv  
    - cellspecific_model.xml  
    - inactive_pathways.csv
    - inactive_pathways_relative.csv

The total number of newly created files is: num_conditions * 3 + 5  