import ruamel.yaml as yaml
import pandas as pd
import argparse
from pathlib import Path


# read configuration from YAML file
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


if __name__ == '__main__':
    description = 'Concatenates all reaction-enumeration solutions'

    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    args = parser.parse_args()

    gene_conditions = [c.strip() for c in doc['gene_expression_columns'].split(',')]
    for condition in gene_conditions:
        solutions = []
        solfiles = Path(cluspath).glob('rxn_enum_solutions_%s_*.csv' % condition)
        print("rxnenumfiles", solfiles)
        for f in solfiles:
            print(f)
            sol = pd.read_csv(f, index_col=0)
            print(type(sol))
            solutions.append(sol)
            print(len(solutions))
        rxn_sols = pd.concat(solutions).drop_duplicates(ignore_index=True)
        rxn_sols.to_csv(cluspath + 'full_rxn_enum_solutions_%s.csv' % condition)
