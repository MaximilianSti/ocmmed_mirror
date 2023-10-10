import dexom_python
import ruamel.yaml as yaml
import pandas as pd
from pathlib import Path
import numpy as np
from cobra.flux_analysis import find_blocked_reactions
from cobra import Configuration
from utilities.force import force_active_rxns, force_reaction_bounds

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

modelpath = doc['modelpath']
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

cobra_config = Configuration()
cobra_config.solver = 'cplex'


def compute_differentially_activated_reactions(control):
    if doc['gene_expression_columns']:
        conditions = [x.strip() for x in doc['gene_expression_columns'].split(',')]
    else:
        genes = pd.read_csv(doc['expressionfile'], sep=';|,|\t', engine='python').set_index(doc['gene_ID_column'])
        conditions = genes.columns.to_list()
    if control not in conditions:
        raise ValueError('No control condition named %s' % control)

    paths = Path(cluspath).glob('all_DEXOM_solutions_*.csv')
    dataframes = {}
    for p in paths:
        if control in str(p):
            ctrl = pd.read_csv(p, index_col=0)
        else:
            df = pd.read_csv(p, index_col=0)
            dataframes[str(p)[20:-4]] = df
    fctrl = ctrl.sum() / len(ctrl)
    darnumbers = pd.Series(index=conditions)
    for c in conditions:
        fc = dataframes[c].sum() / len(c)
        R2c = (fctrl - fc) ** 2
        DARc = fc[R2c > 0.2] - fctrl[R2c > 0.2]
        DARc.to_csv(outpath + 'DAR_significant_' + c + '.csv')
        R2c.to_csv(outpath + 'DAR_R2statistic_' + c + '.csv')
        (fc - fctrl).to_csv(outpath + 'DAR_frequency_diff_' + c + '.csv')
        darnumbers[c] = len(DARc)
    darnumbers.to_csv(outpath + 'DAR_allconditions.csv')
