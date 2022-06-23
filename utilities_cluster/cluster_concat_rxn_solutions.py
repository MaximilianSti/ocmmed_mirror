import ruamel.yaml as yaml
import pandas as pd
import argparse
from pathlib import Path


# read configuration from YAML file
yaml_reader = yaml.YAML(typ='safe')
with open('config.yaml', 'r') as file:
    a = file.read()
doc = yaml_reader.load(a)


if doc['output_path']:
    outpath = doc['output_path']
else:
    outpath = ''


if __name__ == '__main__':
    description = 'Concatenates all reaction-enumeration solutions'

    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    args = parser.parse_args()

    gene_conditions = doc['gene_expression_columns'].split(',')
    for condition in gene_conditions:
        solutions = []
        solfiles = Path(outpath).glob('rxn_enum_solutions_%s_*.csv' % condition)
        for f in solfiles:
            sol = pd.read_csv(f, index_col=0)
            solutions.append(sol)
        solutions = pd.concat(solutions).drop_duplicates(ignore_index=True)
        solutions.to_csv(outpath+'full_rxn_enum_solutions_%s.csv' % condition)
