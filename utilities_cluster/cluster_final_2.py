import dexom_python
import ruamel.yaml as yaml
import pandas as pd
from pathlib import Path
from warnings import warn
from utilities.minimal import mba
from utilities.inactive_pathways import compute_inactive_pathways
from cobra.io import write_sbml_model
from cobra import Configuration

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
    model = dexom_python.read_model(doc['modelpath'], solver='cplex')
    frequencies = pd.read_csv(outpath + 'activation_frequency_reactions.csv', index_col=0)
    freq = frequencies[frequencies.columns[0]]
    if doc['final_network'] == 'union':
        rem_rxns = freq[freq == 0].index.to_list()  # remove reactions which are active in zero solutions
        for rxn in rem_rxns:
            model.remove_reactions([rxn], remove_orphans=True)
    elif doc['final_network'] == 'minimal':
        model = mba(model_keep=model, frequency_table=frequencies, essential_reactions=doc['active_reactions'])
    else:
        warn('Invalid value for "final_network" in parameters.yaml, returning original network.')
    model.id += '_cellspecific'
    print('number of reactions in the model:', len(model.reactions))
    write_sbml_model(model, outpath+'cellspecific_model.xml')
    compute_inactive_pathways(model)
