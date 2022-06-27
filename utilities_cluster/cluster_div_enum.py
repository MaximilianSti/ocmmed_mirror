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
dp = params['div_enum_params']

modelpath = doc['modelpath']


if __name__ == '__main__':

    description = 'For a given condition calculates reaction weights and computes iMAT solution'

    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-c', '--condition', help='column of the gene expression file containing the data for one condition')
    parser.add_argument('-d', '--dist_anneal', type=float, default=-1., help='diversity-enum dist_anneal parameter')
    parser.add_argument('-i', '--maxiter', type=int, default=10, help='number of iterations')
    parser.add_argument('-p', '--parallel_id', type=str, default='', help='id of the parallel thread')
    args = parser.parse_args()
    # read model
    model = dexom_python.read_model(modelpath, solver=mp['solver'])
    model = dexom_python.check_model_options(model, timelimit=mp['timelimit'], feasibility=mp['feasibility'],
                                             mipgaptol=mp['mipgaptol'], verbosity=mp['verbosity'])

    condition = args.condition

    if doc['force_active_reactions']:
        force_active_rxns(model, doc['active_reactions'], doc['fluxvalue'])
    rw = dexom_python.load_reaction_weights(filename=outpath + 'reaction_weights_%s' % condition)

    prevsol, prevbin = dexom_python.read_solution(outpath + 'imat_solution_%s.csv' % condition)
    if clus['approach'] == 'grouped':
        # prevsol, prevbin = dexom_python.read_solution(clus['cluster_files']+'rxn_enum_solution_%s_%s_1.csv' % (condition, args.parallel_id))
        rxnenum_sols = pd.read_csv(outpath + 'rxn_enum_solutions_%s_%s.csv' % (condition, args.parallel_id))
        rxnenum_sols.columns = prevsol.fluxes.index
        prevsol.fluxes = rxnenum_sols.iloc[args.parallel_id]
    elif clus['approach'] == 'separate':
        rxnenum_sols = pd.read_csv(outpath + 'all_rxn_enum_solutions_%s.csv' % condition, index_col=0)
        rxnenum_sols.columns = prevsol.fluxes.index
        prevsol.fluxes = rxnenum_sols.iloc[args.parallel_id]
    else:
        print("could not recognise approach, using iMAT solution as previous solution")

    if args.dist_anneal >= 0:
        distanneal = args.dist_anneal
    else:
        distanneal = dp['dist_anneal']

    divsol, divres = dexom_python.enum_functions.diversity_enum(model=model, reaction_weights=rw, prev_sol=prevsol,
                                dist_anneal=distanneal, eps=ip['epsilon'], thr=ip['threshold'], icut=dp['icut'],
                                maxiter=args.maxiter, obj_tol=ep['objective_tolerance'], full=dp['full'])

    solutions = pd.DataFrame(divsol.binary)
    solutions.to_csv(outpath + 'div_enum_solutions_%s_%s.csv' % (condition, args.parallel_id))
    divres.to_csv(outpath + 'div_enum_stats_%s_%s.csv' % (condition, args.parallel_id))
