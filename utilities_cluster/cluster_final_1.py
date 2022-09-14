import dexom_python
import ruamel.yaml as yaml
import pandas as pd
from pathlib import Path
from cobra import Configuration
import os

yaml_reader = yaml.YAML(typ='safe')
with open('parameters.yaml', 'r') as file:
    a = file.read()
doc = yaml_reader.load(a)

with open('cluster_params.yaml', 'r') as file:
    c = file.read()
clus = yaml_reader.load(c)

if doc['output_path']:
    outpath = doc['output_path']
    if outpath[-1] not in ['/', '\\']:
        outpath += os.sep
else:
    outpath = ''

if clus['cluster_files']:
    cluspath = clus['cluster_files']
    if cluspath[-1] not in ['/', '\\']:
        cluspath += os.sep
else:
    cluspath = outpath

expressionfile = doc['expressionfile']
cobra_config = Configuration()
cobra_config.solver = 'cplex'


if __name__ == '__main__':
    genes = pd.read_csv(expressionfile).set_index(doc['gene_ID_column'])
    if doc['gene_expression_columns']:
        gene_conditions = [x.strip() for x in doc['gene_expression_columns'].split(',')]
    else:
        gene_conditions = genes.columns.to_list()
    all_sols = []
    for condition in gene_conditions:
        solutions = []
        solfiles = Path(cluspath).glob('div_enum_solutions_%s_*.csv' % condition)
        for f in solfiles:
            sol = pd.read_csv(f, index_col=0)
            solutions.append(sol)
        divsols = pd.concat(solutions).drop_duplicates(ignore_index=True)
        divsols.to_csv(cluspath + 'full_div_enum_solutions_%s.csv' % condition)
        rxnsols = pd.read_csv(cluspath + 'full_rxn_enum_solutions_%s.csv' % condition, index_col=0)
        dexomsols = pd.concat([divsols, rxnsols]).drop_duplicates(ignore_index=True)
        dexomsols.to_csv(outpath + 'all_DEXOM_solutions_%s.csv' % condition)
        all_sols.append(dexomsols)
    dex = pd.concat(all_sols).drop_duplicates(ignore_index=True)
    dex.to_csv(outpath + 'all_DEXOM_solutions.csv')
    print("concatenated all DEXOM solutions")
    model = dexom_python.read_model(doc['modelpath'], solver='cplex')
    dex.columns = [r.id for r in model.reactions]
    frequencies = dex.sum()
    frequencies.columns = ['frequency']
    frequencies.to_csv(outpath + 'activation_frequency_reactions.csv')
