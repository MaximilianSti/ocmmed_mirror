import pandas as pd
import numpy as np
from cobra.flux_analysis import find_blocked_reactions


def compute_inactive_pathways(model, fullmodel, outpath, blocked_rxns=None):

    rxns_cell = set([r.id for r in model.reactions])
    rxns_full = set([r.id for r in fullmodel.reactions])
    if blocked_rxns is not None:
        with open(blocked_rxns) as file:
            reader = file.read()
        if '\n' in reader:
            rxns_inactive = set(reader.split('\n'))
        elif ';' in reader:
            rxns_inactive = set(reader.split(';'))
        else:
            rxns_inactive = set(reader.split(','))
    else:
        with open(outpath + 'blocked_reactions.txt') as file:
            rxns_inactive = set(file.read().split('\n'))
    rxns_flux = rxns_full - rxns_inactive
    print('The new model has %i active reactions, out of a maximum of %i reactions' % (len(rxns_cell), len(rxns_flux)))
    print('Percentage of active reactions: {i}%'.format(i=int(100*(len(rxns_cell)/len(rxns_flux)))))
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
