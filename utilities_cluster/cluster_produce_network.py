import dexom_python
import ruamel.yaml as yaml
import pandas as pd
import argparse
from pathlib import Path
from warnings import warn
from utilities.minimal import mba
from cobra.io import write_sbml_model

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
    description = 'Concatenates all dexom solutions and saves the network'

    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    args = parser.parse_args()

    genes = pd.read_csv(expressionpath).set_index(doc['gene_ID_column'])
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

    model = dexom_python.read_model(doc['modelpath'])
    dex.columns = [r.id for r in model.reactions]
    frequencies = dex.sum()
    frequencies.columns = ['frequency']
    frequencies.to_csv(outpath + 'activation_frequency_reactions.csv')

    if doc['final_network'] == 'union':
        rem_rxns = dex.columns[frequencies == 0].to_list()  # remove reactions which are active in zero solutions
        model.remove_reactions(rem_rxns, remove_orphans=True)
    elif doc['final_network'] == 'minimal':
        model = mba(model_keep=model, enum_solutions=dex, essential_reactions=doc['active_reactions'])
    else:
        warn('Invalid value for "final_network" in parameters.yaml, returning original network.')
    model.id += '_cellspecific'
    write_sbml_model(model, outpath+'cellspecific_model.xml')
