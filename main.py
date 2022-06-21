import dexom_python
import ruamel.yaml as yaml
import pandas as pd
import numpy as np
from cobra.io import write_sbml_model
import optlang
from utilities.force import force_active_rxns
from utilities.minimal import mba


# read configuration from YAML file
yaml_reader = yaml.YAML(typ="safe")
with open("config.yaml", "r") as file:
    a = file.read()
doc = yaml_reader.load(a)

with open("additional_params.yaml", "r") as file:
    b = file.read()
params = yaml_reader.load(b)

if doc['output_path']:
    outpath = doc['output_path']
else:
    outpath = ''

mp = params['model_params']
ip = params['imat_params']
ep = params['enum_params']
rp = params['rxn_enum_params']
dp = params['div_enum_params']

modelpath = doc['modelpath']
expressionpath = doc['expressionpath']


if __name__ == "__main__":
    # read model
    model_keep = dexom_python.read_model(modelpath, solver=mp['solver'])
    model_keep = dexom_python.check_model_options(model_keep, timelimit=mp['timelimit'], feasibility=mp['feasibility'],
                                                  mipgaptol=mp['mipgaptol'], verbosity=mp['verbosity'])
    new_model = model_keep.copy()

    # read and process gene expression file
    gene_conditions = doc['gene_expression_columns'].split(',')
    genes = pd.read_csv(expressionpath).set_index(doc['gene_ID_column'])
    if doc['gpr_parameters']['qualitative']:
        genes = dexom_python.expression2qualitative(genes=genes, column_list=gene_conditions,
                                                    proportion=doc['gpr_parameters']['percentile'],
                                                    outpath=outpath+'geneweights_qualitative')

    dexom_sols = []
    for condition in gene_conditions:
        # create reaction weights from gene expression
        print("computing reaction weights for condition "+condition)
        gene_weights = pd.Series(genes[condition].values, index=genes.index, dtype=float)
        for x in set(gene_weights.index):
            if type(gene_weights[x]) != np.float64:
                if len(gene_weights[x].value_counts()) > 1:
                    gene_weights.pop(x)
        gene_weights = gene_weights.to_dict()
        rw = dexom_python.apply_gpr(model=model_keep, gene_weights=gene_weights, modelname=doc['modelname'], save=True,
                                    filename=outpath+'reaction_weights_%s.csv' % condition)

        # compute imat solution from reaction weights
        print("performing iMAT for condition " + condition)
        model = model_keep.copy()
        if doc['force_active_reactions']:
            force_active_rxns(model, doc['active_reactions'], doc['fluxvalue'])
        try:
            imatsol = dexom_python.imat(model=model, reaction_weights=rw, epsilon=ip['epsilon'], threshold=ip['threshold'])
        except optlang.exceptions.SolverError:
            print('Solver could not find a solution with forced active reactions. '
                  'Attempting to find a solution without forced flux.')
            model = model_keep.copy()
            imatsol = dexom_python.imat(model=model, reaction_weights=rw, epsilon=ip['epsilon'],
                                        threshold=ip['threshold'])
            print('The solver could not find a solution with forced active reactions, '
                  'but found one without forcing flux through the selected reactions. '
                  'Check if any of the reactions in the "active_reactions" list are blocked ')
        dexom_python.write_solution(model=model, solution=imatsol, threshold=ip['threshold'],
                                    filename=outpath+'imat_solution_%s.csv' % condition)

        # enumerate solutions with the DEXOM method
        print("performing reaction-enumeration for condition " + condition)
        df = pd.read_csv(rp['reaction_list'], header=None)
        reactions = [x for x in df.unstack().values]
        rxnlist = reactions[:rp['num_rxns']]
        rxnsol = dexom_python.enum_functions.rxn_enum(model=model, reaction_weights=rw, prev_sol=imatsol, rxn_list=rxnlist,
                                                      obj_tol=ep['obj_tol'], thr=ip['threshold'])
        uniques = pd.DataFrame(rxnsol.unique_binary)
        uniques.to_csv(outpath+'rxn_enum_solutions_%s.csv' % condition)
        print("performing diversity-enumeration for condition " + condition)
        divsol, divres = dexom_python.enum_functions.diversity_enum(model=model, reaction_weights=rw, prev_sol=imatsol,
                                                                    eps=ip['epsilon'], thr=ip['threshold'],
                                                                    obj_tol=ep['obj_tol'], maxiter=dp['iterations'],
                                                                    dist_anneal=dp['dist_anneal'])
        divs = pd.DataFrame(divsol.binary)
        divs.to_csv(outpath+'div_enum_solutions_%s.csv' % condition)
        dexomsol = pd.concat([uniques, divs])
        dexom_sols.append(dexomsol)

    dexom_sols = pd.concat(dexom_sols).drop_duplicates(ignore_index=True)
    dexom_sols.to_csv(outpath + 'dexom_solutions.csv')
    dexom_sols.columns = [r.id for r in model_keep.reactions]

    print("producing final network")
    if doc['final_network'] == 'union':
        rem_rxns = dexom_sols.columns[dexom_sols.sum() == 0].to_list()  # remove reactions which are active in zero solutions
        new_model.remove_reactions(rem_rxns, remove_orphans=True)
    elif doc['final_network'] == 'minimal':
        new_model = mba(model_keep=new_model, enum_solutions=dexom_sols, essential_reactions=doc['active_reactions'])
    else:
        print('Invalid value for "final_network" in config.yaml, returning original network.')
    new_model.id += '_cellspecific'
    write_sbml_model(new_model, outpath+'cellspecific_model.csv')
