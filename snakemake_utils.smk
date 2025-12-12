import ruamel.yaml as yaml
import pandas as pd
import os
import dexom_python as dp
import random
import time
from pathlib import Path

# read configuration from YAML files
yaml_reader = yaml.YAML(typ='safe')
with open('parameters.yaml', 'r') as file:
    a = file.read()
params = yaml_reader.load(a)

if params['output_path']:
    outpath = params['output_path']
    os.makedirs(outpath, exist_ok=True)
    if outpath[-1] not in ['/', '\\']:
        outpath += '/'
else:
    outpath = ''

if params['cluster_files']:
    cluspath = params['cluster_files']
    os.makedirs(cluspath, exist_ok=True)
    if cluspath[-1] not in ['/', '\\']:
        cluspath += '/'
else:
    cluspath = outpath


if not Path(params['modelpath']).exists():
    raise FileNotFoundError('Model file not found, check if you provided the correct path: %s' % params['modelpath'])

if not Path(params['expressionfile']).exists():
    raise FileNotFoundError('Gene expression file not found, check if you provided the correct path: %s' % params['expressionfile'])


if params['gene_expression_columns']:
    gene_conditions = [x.strip() for x in params['gene_expression_columns'].split(',')]
else:
    genes = pd.read_csv(params['expressionfile'], sep=';|,|\t', engine='python').set_index(params['gene_ID_column'])
    gene_conditions = genes.columns.to_list()

def get_conditions():
    return gene_conditions

def get_parallel():
    return list(range(params['batch_num']))

if params['rxn_enum_params']['reaction_list']:
    if not Path(params['rxn_enum_params']['reaction_list']).exists():
        raise FileNotFoundError('Reaction-list file not found, check if you provided the correct path: %s' % params['rxn_enum_params']['reaction_list'])
    df = pd.read_csv(params['rxn_enum_params']['reaction_list'], header=None)
    reactions = [x for x in df.unstack().values]
else:
    model = dp.read_model(params['modelpath'])
    reactions = [r.id for r in model.reactions]
    random.shuffle(reactions)
    with open(outpath + 'reactions_shuffled.txt', 'w+') as file:
        file.write('\n'.join(reactions))

if params['blocked_rxns']:
    if not Path(params['blocked_rxns']).exists():
        raise FileNotFoundError('Blocked reaction file not found, check if you provided the correct path: %s' % params['blocked_rxns'])

def get_batchnum():
    batchnum = (len(reactions) // params['batch_rxn_sols']) + 1
    return list(range(batchnum))

final_output_full_rxn_enum = ''
rxn_enum_prefix = 'all_rxn_enum_'
if params['rxn_enum_params']['full_rxn_enum']:
    rxn_enum_prefix += 'full_'
    final_output_full_rxn_enum = expand(cluspath + 'fullrxnenumdone_{condition}.txt', condition=get_conditions())

if isinstance(params['force_flux_bounds'], dict):
    pass
elif params['force_flux_bounds'] and not Path(params['force_flux_bounds']).exists():
    raise FileNotFoundError('Flux bounds file not found, check if you provided the correct path: %s' % ['force_flux_bounds'])

if isinstance(params['force_active_reactions'], dict):
    pass
elif params['force_active_reactions'] and not Path(params['force_active_reactions']).exists():
    raise FileNotFoundError('Flux bounds file not found, check if you provided the correct path: %s' % ['force_active_reactions'])

yaml_writer = yaml.YAML()
with open(outpath + 'parameters_used_for_run_%.0f.yaml' % time.time(), 'w+') as file:
    yaml_writer.dump(params, file)
