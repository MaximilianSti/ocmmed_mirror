import dexom_python
import ruamel.yaml as yaml
import pandas as pd
import numpy as np

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

modelpath = doc['modelpath']
expressionpath = doc['expressionpath']

# read model
model_keep = dexom_python.read_model(modelpath, solver=mp['solver'])
model_keep = dexom_python.check_model_options(model_keep, timelimit=mp['timelimit'], feasibility=mp['feasibility'],
                                              mipgaptol=mp['mipgaptol'], verbosity=mp['verbosity'])

# read and process gene expression file
gene_conditions = doc['gene_expression_columns'].split(',')
genes = pd.read_csv(expressionpath).set_index(doc['gene_ID_column'])
if doc['gpr_parameters']['qualitative']:
    genes = dexom_python.expression2qualitative(genes=genes, column_list=gene_conditions,
                                                proportion=doc['gpr_parameters']['percentile'],
                                                outpath=outpath+'geneweights_qualitative')

# create reaction weights from gene expression
reaction_weight_dict = {}
imat_solution_dict = {}
for condition in gene_conditions:
    print("computing reaction weights for condition "+condition)
    gene_weights = pd.Series(genes[condition].values, index=genes.index, dtype=float)
    for x in set(gene_weights.index):
        if type(gene_weights[x]) != np.float64:
            if len(gene_weights[x].value_counts()) > 1:
                gene_weights.pop(x)
    gene_weights = gene_weights.to_dict()
    rw = dexom_python.apply_gpr(model=model_keep, gene_weights=gene_weights, modelname=doc['modelname'], save=True,
                                filename=outpath+'reaction_weights_%s.csv' % condition)
    reaction_weight_dict[condition] = rw

    # compute imat solution from reaction weights
    print("performing iMAT for condition " + condition)
    model = model_keep.copy()
    imatsol = dexom_python.imat(model=model, reaction_weights=rw, epsilon=ip['epsilon'], threshold=ip['threshold'])
    imat_solution_dict[condition] = imatsol
    dexom_python.write_solution(model=model, solution=imatsol, threshold=ip['threshold'],
                                filename=outpath+'imat_solution_%s.csv' % condition)

    print("performing reaction-enumeration for condition " + condition)
    df = pd.read_csv(rp['reaction_list'], header=None)
    reactions = [x for x in df.unstack().values]
    rxnlist = reactions[:rp['num_rxns']]
    rxnsol = dexom_python.enum_functions.rxn_enum(model=model, reaction_weights=rw, prev_sol=imatsol, rxn_list=rxnlist,
                                                  obj_tol=ep['obj_tol'])
    uniques = pd.DataFrame(rxnsol.unique_binary)
    uniques.to_csv(outpath+'rxn_enum_solutions_%s.csv' % condition)
