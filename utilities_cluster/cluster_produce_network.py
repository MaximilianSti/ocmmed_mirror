import dexom_python
import ruamel.yaml as yaml
import pandas as pd
from warnings import warn
from utilities.minimal import mba
from utilities.inactive_pathways import compute_inactive_pathways
from cobra.io import write_sbml_model
from cobra import Configuration

yaml_reader = yaml.YAML(typ='safe')
with open('parameters.yaml', 'r') as file:
    a = file.read()
doc = yaml_reader.load(a)

if doc['output_path']:
    outpath = doc['output_path']
else:
    outpath = ''

cobra_config = Configuration()
cobra_config.solver = 'gurobi'


if __name__ == '__main__':
    description = 'Concatenates all dexom solutions and saves the network'
    dex = pd.read_csv(outpath + 'all_DEXOM_solutions.csv', index_col=0)
    print("loaded all DEXOM solutions")
    model = dexom_python.read_model(doc['modelpath'], solver='gurobi')
    dex.columns = [r.id for r in model.reactions]
    frequencies = dex.sum()
    frequencies.columns = ['frequency']
    frequencies.to_csv(outpath + 'activation_frequency_reactions.csv')
    print("saved frequencies")
    if doc['final_network'] == 'union':
        rem_rxns = dex.columns[frequencies == 0].to_list()  # remove reactions which are active in zero solutions
        model.remove_reactions(rem_rxns, remove_orphans=True)
    elif doc['final_network'] == 'minimal':
        model = mba(model_keep=model, enum_solutions=dex, essential_reactions=doc['active_reactions'])
    else:
        warn('Invalid value for "final_network" in parameters.yaml, returning original network.')
    model.id += '_cellspecific'
    print('number of reactions in the model:', len(model.reactions))
    write_sbml_model(model, outpath+'cellspecific_model.xml')
    compute_inactive_pathways(model)
