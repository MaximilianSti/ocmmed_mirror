import dexom_python
import ruamel.yaml as yaml
import pandas as pd
import numpy as np
import optlang
from utilities.force import force_active_rxns
import argparse


# read configuration from YAML file
yaml_reader = yaml.YAML(typ='safe')
with open('../config.yaml', 'r') as file:
    a = file.read()
doc = yaml_reader.load(a)

with open('../additional_params.yaml', 'r') as file:
    b = file.read()
params = yaml_reader.load(b)

if doc['output_path']:
    outpath = doc['output_path']
else:
    outpath = ''

mp = params['model_params']
ip = params['imat_params']
ep = params['enum_params']

modelpath = doc['modelpath']
expressionpath = doc['expressionpath']


if __name__ == '__main__':

    description = 'For a given condition calculates reaction weights and computes iMAT solution'

    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-c', '--condition', help='column of the gene expression file containing the data for one condition')
    args = parser.parse_args()
    # read model
    model = dexom_python.read_model(modelpath, solver=mp['solver'])
    model = dexom_python.check_model_options(model, timelimit=mp['timelimit'], feasibility=mp['feasibility'],
                                             mipgaptol=mp['mipgaptol'], verbosity=mp['verbosity'])

    condition = args.condition
    # read and process gene expression file
    genes = pd.read_csv(expressionpath).set_index(doc['gene_ID_column'])
    if doc['gpr_parameters']['qualitative']:
        genes = dexom_python.expression2qualitative(genes=genes, column_list=[condition],
                                                    proportion=doc['gpr_parameters']['percentile'],
                                                    outpath=outpath+'geneweights_qualitative')

    # create reaction weights from gene expression
    print('computing reaction weights for condition '+condition)
    gene_weights = pd.Series(genes[condition].values, index=genes.index, dtype=float)
    for x in set(gene_weights.index):
        if type(gene_weights[x]) != np.float64:
            if len(gene_weights[x].value_counts()) > 1:
                gene_weights.pop(x)
    gene_weights = gene_weights.to_dict()
    rw = dexom_python.apply_gpr(model=model, gene_weights=gene_weights, modelname=doc['modelname'], save=True,
                                filename=outpath+'reaction_weights_%s' % condition)

    # compute imat solution from reaction weights
    print('performing iMAT for condition ' + condition)

    if doc['force_active_reactions']:
        force_active_rxns(model, doc['active_reactions'], doc['fluxvalue'])
    try:
        imatsol = dexom_python.imat(model=model, reaction_weights=rw, epsilon=ip['epsilon'], threshold=ip['threshold'])
    except optlang.exceptions.SolverError:
        print('Solver could not find a solution with forced active reactions. '
              'Attempting to find a solution without forced flux.')
        imatsol = dexom_python.imat(model=model, reaction_weights=rw, epsilon=ip['epsilon'],
                                    threshold=ip['threshold'])
        print('The solver could not find a solution with forced active reactions, '
              'but found one without forcing flux through the selected reactions. '
              'Check if any of the reactions in the "active_reactions" list are blocked ')
    dexom_python.write_solution(model=model, solution=imatsol, threshold=ip['threshold'],
                                filename=outpath+'imat_solution_%s.csv' % condition)
