import dexom_python
import ruamel.yaml as yaml
import pandas as pd
from utilities.force import force_active_rxns, force_reaction_bounds
import argparse
from warnings import warn
import os

yaml_reader = yaml.YAML(typ='safe')
with open('parameters.yaml', 'r') as file:
    a = file.read()
doc = yaml_reader.load(a)

with open('additional_params.yaml', 'r') as file:
    b = file.read()
params = yaml_reader.load(b)

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

mp = params['model_params']
ip = params['imat_params']
ep = params['enum_params']
rp = params['rxn_enum_params']

modelpath = doc['modelpath']

if __name__ == '__main__':
    description = 'For a given condition, performs reaction-enumeration'
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-c', '--condition',
                        help='column of the gene expression file containing the data for one condition')
    parser.add_argument('-r', '--rxn_range', default='_', help='range of reactions')
    parser.add_argument('-p', '--parallel_id', type=str, default='', help='id of the parallel thread')
    args = parser.parse_args()
    # read model
    model = dexom_python.read_model(modelpath, solver=mp['solver'])
    model = dexom_python.check_model_options(model, timelimit=mp['timelimit'], feasibility=mp['feasibility'],
                                             mipgaptol=mp['mipgaptol'], verbosity=mp['verbosity'])
    condition = args.condition

    if doc['force_active_reactions']:
        force_active_rxns(model, doc['active_reactions'], doc['fluxvalue'])
    if doc['force_flux_bounds']:
        force_reaction_bounds(model, doc['force_flux_bounds'])

    rw = dexom_python.load_reaction_weights(filename=outpath+'reaction_weights_%s.csv' % condition)
    imatsol, imatbin = dexom_python.read_solution(filename=outpath+'imat_solution_%s.csv' % condition)

    if rp['reaction_list']:
        df = pd.read_csv(rp['reaction_list'], header=None)
        reactions = [x for x in df.unstack().values]
        if False in [rid in model.reactions for rid in reactions]:
            warn('One or more of the reactions in "reaction_list" are not present in the model. '
                 'The provided reaction list will not be used.')
            reactions = [r.id for r in model.reactions]
    else:
        reactions = [r.id for r in model.reactions]
    if params['blocked_rxns']:
        with open(params['blocked_rxns']) as file:
            reader = file.read()
        if '\n' in reader:
            rxns_inactive = set(reader.split('\n'))
        elif ';' in reader:
            rxns_inactive = set(reader.split(';'))
        else:
            rxns_inactive = set(reader.split(','))
        reactions = list(set(reactions) - rxns_inactive)
    rxn_range = args.rxn_range.split('_')
    if rxn_range[0] == '':
        start = 0
    else:
        start = int(rxn_range[0])
    if rxn_range[1] == '':
        rxn_list = reactions[start:]
    elif int(rxn_range[1]) >= len(reactions):
        rxn_list = reactions[start:]
    else:
        rxn_list = reactions[start:int(rxn_range[1])]
    rxnsol = dexom_python.enum_functions.rxn_enum(model=model, reaction_weights=rw, prev_sol=imatsol, rxn_list=rxn_list,
                                                  obj_tol=ep['objective_tolerance'], eps=ip['epsilon'], thr=ip['threshold'])
    uniques = pd.DataFrame(rxnsol.unique_binary)
    uniques.to_csv(cluspath + 'rxn_enum_solutions_%s_%s.csv' % (condition, args.parallel_id))

    # if clus['approach'] == 'grouped':
    #     for i in range(1, len(rxnsol.unique_solutions)):
    #         dexom_python.write_solution(model, rxnsol.unique_solutions[i], ip['threshold'], clus['cluster_files']
    #                                     +'rxn_enum_solution_%s_%s_%i.csv' % (condition, args.parallel_id, i))
