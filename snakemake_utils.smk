import ruamel.yaml as yaml
import pandas as pd
import os
import dexom_python as dp

# read configuration from YAML files
yaml_reader = yaml.YAML(typ='safe')
with open('parameters.yaml', 'r') as file:
    a = file.read()
doc = yaml_reader.load(a)

with open('params_additional.yaml', 'r') as file:
    b = file.read()
pad = yaml_reader.load(b)

with open('params_cluster.yaml', 'r') as file:
    c = file.read()
clus = yaml_reader.load(c)

if doc['output_path']:
    outpath = doc['output_path']
    os.makedirs(outpath, exist_ok=True)
    if outpath[-1] not in ['/', '\\']:
        outpath += '/'
else:
    outpath = ''

if clus['cluster_files']:
    cluspath = clus['cluster_files']
    os.makedirs(cluspath, exist_ok=True)
    if cluspath[-1] not in ['/', '\\']:
        cluspath += '/'
else:
    cluspath = outpath

    print('writing Snakefile')

if doc['gene_expression_columns']:
    gene_conditions = [x.strip() for x in doc['gene_expression_columns'].split(',')]
else:
    genes = pd.read_csv(doc['expressionfile'], sep=';|,|\t', engine='python').set_index(doc['gene_ID_column'])
    gene_conditions = genes.columns.to_list()

def get_conditions():
    return gene_conditions

def get_parallel():
    return list(range(clus['batch_num']))

