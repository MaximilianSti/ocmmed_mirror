import ruamel.yaml as yaml
import pandas as pd
import os
import dexom_python as dp
from pathlib import Path

# read configuration from YAML files
yaml_reader = yaml.YAML(typ='safe')
with open('parameters.yaml', 'r') as file:
    a = file.read()
params = yaml_reader.load(a)

if params['output_path']:
    outpath = params['output_path']
    if outpath[-1] not in ['/', '\\']:
        outpath += '/'
else:
    outpath = ''

cluspath = outpath[:-1] + 'clusterfiles/'

# check if modelpath exists
if not Path(params['modelpath']).exists():
    raise FileNotFoundError('Model file not found, check if you provided the correct path: %s' % params['modelpath'])

# check if expressionfile exists
if not Path(params['expressionfile']).exists():
    raise FileNotFoundError('Gene expression file not found, check if you provided the correct path: %s' % params['expressionfile'])

# helper function for getting conditions
if params['gene_expression_columns']:
    gene_conditions = [x.strip() for x in params['gene_expression_columns'].split(',')]
else:
    genes = pd.read_csv(params['expressionfile'], sep=';|,|\t', engine='python').set_index(params['gene_ID_column'])
    gene_conditions = genes.columns.to_list()

def get_conditions():
    return gene_conditions

# helper function for getting batch numbers
def get_parallel():
    return list(range(params['batch_num']))

# check if blocked_rxns exists
if params['blocked_rxns']:
    if not Path(params['blocked_rxns']).exists():
        raise FileNotFoundError('Blocked reaction file not found, check if you provided the correct path: %s' % params['blocked_rxns'])


# helper function for getting full_rxn_eum batch numbers, check if reaction_list exists
if params['rxn_enum_params']['reaction_list']:
    if not Path(params['rxn_enum_params']['reaction_list']).exists():
        raise FileNotFoundError('Reaction-list file not found, check if you provided the correct path: %s' % params['rxn_enum_params']['reaction_list'])
    df = pd.read_csv(params['rxn_enum_params']['reaction_list'], header=None)
    reactionlen = len([x for x in df.unstack().values])
else:
    if Path(outpath + 'reactions_shuffled.txt').exists():
        with open(outpath+'reactions_shuffled.txt', 'r') as f:
            reactionlen = len(f.read().split('\n'))
    else:
        model = dp.read_model(params['modelpath'])
        reactionlen = len(model.reactions)

def get_batchnum():
    batchnum = (reactionlen // params['batch_rxn_sols']) + 1
    return list(range(batchnum))

final_output_full_rxn_enum = ''
rxn_enum_prefix = 'all_rxn_enum_'
if params['full_rxn_enum']:
    rxn_enum_prefix += 'full_'
    final_output_full_rxn_enum = expand(cluspath + 'fullrxnenumdone_{condition}.txt', condition=get_conditions())

# check if force_flux_bounds exists
if isinstance(params['force_flux_bounds'], dict):
    pass
elif params['force_flux_bounds'] and not Path(params['force_flux_bounds']).exists():
    raise FileNotFoundError('Flux bounds file not found, check if you provided the correct path: %s' % ['force_flux_bounds'])

# check if force_active_reactions exists
if isinstance(params['force_active_reactions'], dict):
    pass
elif params['force_active_reactions'] and not Path(params['force_active_reactions']).exists():
    raise FileNotFoundError('Flux bounds file not found, check if you provided the correct path: %s' % ['force_active_reactions'])
