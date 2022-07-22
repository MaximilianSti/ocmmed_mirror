import ruamel.yaml as yaml
import pandas as pd
import argparse
from pathlib import Path

yaml_reader = yaml.YAML(typ='safe')
with open('parameters.yaml', 'r') as file:
    a = file.read()
doc = yaml_reader.load(a)

with open('cluster_params.yaml', 'r') as file:
    c = file.read()
clus = yaml_reader.load(c)

if doc['output_path']:
    outpath = doc['output_path']
else:
    outpath = ''

if clus['cluster_files']:
    cluspath = clus['cluster_files']
else:
    cluspath = outpath

expressionpath = doc['expressionpath']


if __name__ == '__main__':
    description = 'Concatenates all reaction-enumeration solutions'
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    args = parser.parse_args()

    genes = pd.read_csv(expressionpath).set_index(doc['gene_ID_column'])
    if doc['gene_expression_columns']:
        gene_conditions = [x.strip() for x in doc['gene_expression_columns'].split(',')]
    else:
        gene_conditions = genes.columns.to_list()
    for condition in gene_conditions:
        solutions = []
        solfiles = Path(cluspath).glob('rxn_enum_solutions_%s_*.csv' % condition)
        for f in solfiles:
            sol = pd.read_csv(f, index_col=0)
            solutions.append(sol)
        rxn_sols = pd.concat(solutions).drop_duplicates(ignore_index=True)
        rxn_sols.to_csv(cluspath + 'full_rxn_enum_solutions_%s.csv' % condition)
