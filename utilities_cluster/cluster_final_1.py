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

with open('params_cluster.yaml', 'r') as file:
    c = file.read()
clus = yaml_reader.load(c)

if doc['output_path']:
    outpath = doc['output_path']
    if outpath[-1] not in ['/', '\\']:
        outpath += '/'
else:
    outpath = ''

if clus['cluster_files']:
    cluspath = clus['cluster_files']
    if cluspath[-1] not in ['/', '\\']:
        cluspath += '/'
else:
    cluspath = outpath

expressionfile = doc['expressionfile']
cobra_config = Configuration()
cobra_config.solver = 'cplex'


if __name__ == '__main__':
    genes = pd.read_csv(expressionfile, sep=';|,|\t', engine='python').set_index(doc['gene_ID_column'])
    genes = genes.loc[genes.index.dropna()]
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
        div_sols = pd.concat(solutions, ignore_index=True).drop_duplicates()
        
        fluxes = []
        fluxfiles = Path(cluspath).glob('div_enum_fluxes_%s_*.csv' % condition)
        for f in fluxfiles:
            fl = pd.read_csv(f, index_col=0)
            fluxes.append(fl)
        div_fluxes = pd.concat(fluxes, ignore_index=True).loc[div_sols.index]

        div_sols.reset_index(inplace=True, drop=True)
        div_fluxes.reset_index(inplace=True, drop=True)
        
        div_sols.to_csv(cluspath + 'full_div_enum_solutions_%s.csv' % condition)
        div_fluxes.to_csv(cluspath + 'full_div_enum_fluxes_%s.csv' % condition)
        rxn_sols = pd.read_csv(cluspath + 'full_rxn_enum_solutions_%s.csv' % condition, index_col=0)
        dexomsols = pd.concat([div_sols, rxn_sols]).drop_duplicates(ignore_index=True)
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
