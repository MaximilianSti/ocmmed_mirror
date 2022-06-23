import dexom_python
import ruamel.yaml as yaml
import pandas as pd
import argparse
from pathlib import Path
from warnings import warn
from utilities.minimal import mba
from cobra.io import write_sbml_model


# read configuration from YAML file
yaml_reader = yaml.YAML(typ='safe')
with open('config.yaml', 'r') as file:
    a = file.read()
doc = yaml_reader.load(a)


if doc['output_path']:
    outpath = doc['output_path']
else:
    outpath = ''


if __name__ == '__main__':
    description = 'Concatenates all dexom solutions and saves the network'

    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    args = parser.parse_args()

    gene_conditions = doc['gene_expression_columns'].split(',')
    all_sols = []
    for condition in gene_conditions:
        solutions = []
        solfiles = Path(outpath).glob('div_enum_solutions_%s_*.csv' % condition)
        for f in solfiles:
            sol = pd.read_csv(f, index_col=0)
            solutions.append(sol)
        divsols = pd.concat(solutions).drop_duplicates(ignore_index=True)
        divsols.to_csv(outpath + 'full_div_enum_solutions_%s.csv' % condition)
        rxnsols = pd.read_csv(outpath + 'full_rxn_enum_solutions_%s.csv' % condition, index_col=0)
        dexomsols = pd.concat([divsols, rxnsols]).drop_duplicates(ignore_index=True)
        dexomsols.to_csv(outpath + 'all_DEXOM_solutions_%s.csv' % condition)
        all_sols.append(dexomsols)
    dex = pd.concat(all_sols).drop_duplicates(ignore_index=True)

    model = dexom_python.read_model(doc['modelpath'])
    dex.columns = [r.id for r in model.reactions]

    if doc['final_network'] == 'union':
        rem_rxns = dex.columns[dex.sum() == 0].to_list()  # remove reactions which are active in zero solutions
        model.remove_reactions(rem_rxns, remove_orphans=True)
    elif doc['final_network'] == 'minimal':
        model = mba(model_keep=model, enum_solutions=dex, essential_reactions=doc['active_reactions'])
    else:
        warn('Invalid value for "final_network" in config.yaml, returning original network.')
    model.id += '_cellspecific'
    write_sbml_model(model, outpath+'cellspecific_model.xml')
