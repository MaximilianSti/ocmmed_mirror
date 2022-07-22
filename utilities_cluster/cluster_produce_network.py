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

    dex = pd.read_csv(outpath + 'all_DEXOM_solutions.csv')
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
