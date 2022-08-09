import dexom_python
import ruamel.yaml as yaml
import pandas as pd
import numpy as np
from cobra.flux_analysis import find_blocked_reactions
from cobra import Configuration
import os

yaml_reader = yaml.YAML(typ='safe')
with open('parameters.yaml', 'r') as file:
    a = file.read()
doc = yaml_reader.load(a)
with open('additional_params.yaml', 'r') as file:
    b = file.read()
params = yaml_reader.load(b)
if doc['output_path']:
    outpath = doc['output_path']
    if outpath[-1] not in ['/', '\\']:
        outpath += os.sep
else:
    outpath = ''
modelpath = doc['modelpath']

cobra_config = Configuration()
cobra_config.solver = 'gurobi'


def compute_inactive_pathways(model):
    fullmodel = dexom_python.read_model(modelpath)
    rxns_cell = set([r.id for r in model.reactions])
    rxns_full = set([r.id for r in fullmodel.reactions])
    if params['blocked_rxns']:
        with open(params['blocked_rxns']) as file:
            reader = file.read()
        if '\n' in reader:
            rxns_inactive = set(reader.split('\n'))
        elif ';' in reader:
            rxns_inactive = set(reader.split(';'))
        else:
            rxns_inactive = set(reader.split(','))
    else:
        potential_inactive = [r.id for r in fullmodel if r.id not in model.reactions]
        if hasattr(fullmodel, "_sbml"):
            fullmodel._sbml['created'] = None
        rxns_inactive = set(find_blocked_reactions(fullmodel, reaction_list=potential_inactive))
    rxns_flux = rxns_full - rxns_inactive
    paths = pd.Series(dtype=int)
    for g in fullmodel.groups:
        paths[g.name] = len([m for m in g.members if m.id in (rxns_flux - rxns_cell)])
    paths.sort_values(ascending=False, inplace=True)
    paths.to_csv(outpath + 'inactive_pathways.csv')
    paths_rel = pd.Series(dtype=float)
    for g in fullmodel.groups:
        paths_rel[g.name] = np.around(100 * paths[g.name] / len(g.members), 1)
    paths_rel.sort_values(ascending=False, inplace=True)
    paths_rel.to_csv(outpath + 'inactive_pathways_relative.csv')
    return paths, paths_rel
