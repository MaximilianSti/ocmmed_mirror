# OCMMED
Obtaining cell-specific metabolic models through enumeration with DEXOM

## Requirements
- Python 3.7+
- dexom-python
- cplex 12.10+

Install dexom-python with:

`pip install dexom-python`

## Inputs
The `config.yaml` file contains the main parameters, including the path to the metabolic model and the path to the gene expression file  
`additional_params.yaml` contains more optional parameters

## Main pipeline
The gene expression data is converted into one set of reaction weights for each experimental condition.  
These weights are then used to find a context-specific network which fits the transcriptomics.  
Afterwards, DEXOM is used to enumerate multiple solutions (reaction-enumeration & diversity-enumeration)
The DEXOM solutions are merged and used to construct the new metabolic model which is saved as an SBML file.

## Cluster pipeline
Because the enumeration of multiple solutions with DEXOM can be very slow when using large metabolic networks and/or large transcriptomic datasets, the pipeline has been adapted for use on a slurm computation cluster.  
Parameters for the cluster pipeline can be found in `cluster_params.yaml`  

## Todo
write documentation for cluster pipeline
