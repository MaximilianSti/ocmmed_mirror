import dexom_python
import ruamel.yaml as yaml
import pandas as pd
from utilities.force import force_active_rxns
import argparse


# read configuration from YAML file
yaml_reader = yaml.YAML(typ='safe')
with open('config.yaml', 'r') as file:
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
else:
    outpath = ''

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
    parser.add_argument('-a', '--approach', type=int, default=2, help='DEXOM approach')
    args = parser.parse_args()
    # read model
    model = dexom_python.read_model(modelpath, solver=mp['solver'])
    model = dexom_python.check_model_options(model, timelimit=mp['timelimit'], feasibility=mp['feasibility'],
                                             mipgaptol=mp['mipgaptol'], verbosity=mp['verbosity'])
    condition = args.condition

    if doc['force_active_reactions']:
        force_active_rxns(model, doc['active_reactions'], doc['fluxvalue'])

    rw = dexom_python.load_reaction_weights(filename=outpath+'reaction_weights_%s' % condition)
    imatsol, imatbin = dexom_python.read_solution(filename=outpath+'imat_solution_%s.csv' % condition)

    if rp['reaction_list']:
        df = pd.read_csv(rp['reaction_list'], header=None)
        reactions = [x for x in df.unstack().values]
    else:
        reactions = [r.id for r in model.reactions]

    rxn_range = args.rxn_range.split('_')
    if rxn_range[0] == '':
        start = 0
    else:
        start = int(rxn_range[0])
    if rxn_range[1] == '':
        rxn_list = reactions[start:]
    elif int(rxn_range[1]) > len(reactions):
        rxn_list = reactions[start:]
    else:
        rxn_list = reactions[start:int(rxn_range[1])]
    rxnsol = dexom_python.enum_functions.rxn_enum(model=model, reaction_weights=rw, prev_sol=imatsol, rxn_list=rxn_list,
                                                  obj_tol=ep['objective_tolerance'], thr=ip['threshold'])
    uniques = pd.DataFrame(rxnsol.unique_binary)
    uniques.to_csv(outpath + 'rxn_enum_solutions_%s_%s.csv' % (condition, args.parallel_id))

    if args.approach == 2:
        for i in range(1, len(rxnsol.unique_solutions)):
            dexom_python.write_solution(model, rxnsol.unique_solutions[i], ip['threshold'], clus['cluster_files']
                                        +'rxn_enum_solution_%s_%s_%i.csv' % (condition, args.parallel_id, i))
