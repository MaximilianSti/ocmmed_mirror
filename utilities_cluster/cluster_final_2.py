import dexom_python
import ruamel.yaml as yaml
import pandas as pd
from utilities.minimal import maximal_frequency
from utilities.force import force_active_rxns, force_reaction_bounds
from utilities.inactive_pathways import compute_inactive_pathways
from utilities.differentially_activated_reactions import compute_differentially_activated_reactions
from cobra.io import write_sbml_model
from cobra.flux_analysis import find_blocked_reactions
from cobra import Configuration
from warnings import warn

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

cluspath = outpath[:-1] + '_clusterfiles/'

expressionfile = params['expressionfile']
cobra_config = Configuration()
cobra_config.solver = 'cplex'


if __name__ == '__main__':
    model = dexom_python.read_model(params['modelpath'], solver='cplex')
    if isinstance(params['force_flux_bounds'], dict):
        force_reaction_bounds(model, params['force_flux_bounds'])
    if isinstance(params['force_active_reactions'], list):
        force_active_rxns(model, params['force_active_reactions'], params['fluxvalue'])
    frequencies = pd.read_csv(outpath + 'activation_frequency_reactions.csv', index_col=0)
    freq = frequencies[frequencies.columns[0]]
    if params['final_network'] == 'union':
        cutoff = params['union_cutoff']
        if isinstance(cutoff, str):
            if cutoff[-1] == '%':
                cutoff = freq.max() * float(cutoff[:-1]) / 100
            else:
                warn('Unrecognized character in union_cutoff parameter, default to 0.')
                cutoff = 0
        rem_rxns = freq[freq <= cutoff].index.to_list()  # remove reactions which are active in less than union_cutoff solutions
        model.remove_reactions(rem_rxns, remove_orphans=True)
        if cutoff > 0:
            blocked = find_blocked_reactions(model)
            model.remove_reactions(blocked, remove_orphans=True)
    elif params['final_network'] == 'minimal':
        model = maximal_frequency(model_keep=model, frequency_table=frequencies, essential_reactions=params['force_active_reactions'])
    elif params['final_network'] == 'none':
        pass
    else:
        raise ValueError('Invalid value for "final_network" in parameters.yaml.')
    model.id += '_cellspecific'
    write_sbml_model(model, outpath+'cellspecific_model.xml')

    fullmodel = dexom_python.read_model(params['modelpath'], solver='cplex')
    compute_inactive_pathways(model=model, fullmodel=fullmodel, outpath=outpath, blocked_rxns=params['blocked_rxns'])

    if params['gene_expression_columns']:
        conditions = [x.strip() for x in params['gene_expression_columns'].split(',')]
    else:
        genes = pd.read_csv(params['expressionfile'], sep=';|,|\t', engine='python').set_index(params['gene_ID_column'])
        conditions = genes.columns.to_list()
    compute_differentially_activated_reactions(inputpath=outpath, conditions=conditions)
