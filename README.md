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

## Pipeline
The gene expression data is converted into one set of reaction weights for each experimental condition.  
These weights are then used to find a context-specific network which fits the transcriptomics.  
Afterwards, DEXOM is used to enumerate multiple solutions (reaction-enumeration & diversity-enumeration)
The DEXOM solutions are merged and used to construct the new metabolic model which is saved as an SBML file.

## Todo
integrate slurm cluster pipeline?
