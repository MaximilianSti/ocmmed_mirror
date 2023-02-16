import dexom_python
import ruamel.yaml as yaml
import pandas as pd
import optlang
from utilities.force import force_active_rxns, force_reaction_bounds
import argparse
from warnings import warn, filterwarnings, resetwarnings, catch_warnings

yaml_reader = yaml.YAML(typ='safe')
with open('parameters.yaml', 'r') as file:
    a = file.read()
doc = yaml_reader.load(a)

with open('params_additional.yaml', 'r') as file:
    b = file.read()
params = yaml_reader.load(b)

with open('params_cluster.yaml', 'r') as file:
    c = file.read()
clus = yaml_reader.load(c)

if doc['output_path']:
    outpath = doc['output_path']
    if outpath[-1] not in ['/', '\\']:
        outpath += '/'
else:
    outpath = ''

mp = params['model_params']
ip = params['imat_params']
ep = params['enum_params']

modelpath = doc['modelpath']
expressionfile = doc['expressionfile']


if __name__ == '__main__':
    description = 'For a given condition calculates reaction weights and computes iMAT solution'
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-c', '--condition', help='column of the gene expression file containing the data for one condition')
    args = parser.parse_args()
    # read model
    model_keep = dexom_python.read_model(modelpath, solver=mp['solver'])
    model_keep = dexom_python.check_model_options(model_keep, timelimit=mp['timelimit'], feasibility=mp['feasibility'],
                                                  mipgaptol=mp['mipgaptol'], verbosity=mp['verbosity'])
    condition = args.condition
    # read and process gene expression file
    genes = pd.read_csv(expressionfile, sep=';|,|\t', engine='python').set_index(doc['gene_ID_column'])
    genes = genes.loc[genes.index.dropna()]
    if doc['gpr_parameters']['qualitative'] and not doc['reaction_scores']:
        genes = dexom_python.expression2qualitative(genes=genes, column_list=[condition],
                                                    proportion=doc['gpr_parameters']['percentile'],
                                                    outpath=outpath+'geneweights_qualitative_%s' % condition)
    # create reaction weights from gene expression
    model = model_keep.copy()
    print('computing reaction weights for condition '+condition)
    gene_weights = pd.Series(genes[condition].values, index=genes.index, dtype=float)
    if doc['reaction_scores']:
        rw = {}
        for rxn in model.reactions:
            rw[rxn.id] = float(gene_weights.to_dict().get(rxn.id, 0.))
        dexom_python.save_reaction_weights(rw, outpath + 'reaction_weights_%s.csv' % condition)
    else:
        rw = dexom_python.apply_gpr(model=model, gene_weights=gene_weights, duplicates=doc['duplicates'], save=True,
                                    filename=outpath+'reaction_weights_%s' % condition)

    # compute imat solution from reaction weights

    print('performing iMAT for condition ' + condition)
    if doc['force_flux_bounds']:
        force_reaction_bounds(model, doc['force_flux_bounds'])
    if doc['force_active_reactions']:
        force_active_rxns(model, doc['force_active_reactions'], doc['fluxvalue'])

    solver_ready = True
    if clus['force_cplex'] and not hasattr(optlang, 'cplex_interface'):
        solver_ready = False

    if solver_ready:
        with catch_warnings():
            filterwarnings('error')
            try:
                imatsol = dexom_python.imat(model=model, reaction_weights=rw, epsilon=ip['epsilon'], threshold=ip['threshold'])
            except UserWarning:
                resetwarnings()
                warn('Solver could not find a solution with forced active reactions. '
                     'Attempting to find a solution without forced flux.')
                model = model_keep.copy()
                imatsol = dexom_python.imat(model=model, reaction_weights=rw, epsilon=ip['epsilon'],
                                            threshold=ip['threshold'])
                print('The solver could not find a solution with forced active reactions, '
                      'but found one without forcing flux through the selected reactions. '
                      'Check if any of the reactions in the "active_reactions" list are blocked ')
            dexom_python.write_solution(model=model, solution=imatsol, threshold=ip['threshold'],
                                        filename=outpath+'imat_solution_%s.csv' % condition)
    else:
        warn('cplex is not properly installed and the force_cplex parameter is set to True')
        with open(outpath + 'imat_CPLEXERROR_%s.txt' % condition, 'w+') as file:
            file.write('cplex is not properly installed and the force_cplex parameter is set to True')
