import dexom_python
import ruamel.yaml as yaml
import pandas as pd
from utilities.minimal import mba
from utilities.force import force_active_rxns, force_reaction_bounds
from utilities.inactive_pathways import compute_inactive_pathways
from cobra.io import write_sbml_model
from cobra import Configuration
from warnings import warn
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
    model = dexom_python.read_model(doc['modelpath'], solver='cplex')
    if doc['force_flux_bounds']:
        force_reaction_bounds(model, doc['force_flux_bounds'])
    if doc['force_active_reactions']:
        force_active_rxns(model, doc['force_active_reactions'], doc['fluxvalue'])
    frequencies = pd.read_csv(outpath + 'activation_frequency_reactions.csv', index_col=0)
    freq = frequencies[frequencies.columns[0]]
    if doc['final_network'] == 'union':
        cutoff = doc['union_cutoff']
        if isinstance(cutoff, str):
            if cutoff[-1] == '%':
                cutoff = len(freq) * float(cutoff[:-1]) / 100
            else:
                warn('Unrecognized character in union_cutoff parameter, default to 0.')
                cutoff = 0
        rem_rxns = freq[freq <= cutoff].index.to_list()  # remove reactions which are active in less than union_cutoff solutions
        for rxn in rem_rxns:
            model.remove_reactions([rxn], remove_orphans=True)
    elif doc['final_network'] == 'minimal':
        model = mba(model_keep=model, frequency_table=frequencies, essential_reactions=doc['force_active_reactions'])
    else:
        raise ValueError('Invalid value for "final_network" in parameters.yaml.')
    model.id += '_cellspecific'
    write_sbml_model(model, outpath+'cellspecific_model.xml')
    compute_inactive_pathways(model)
