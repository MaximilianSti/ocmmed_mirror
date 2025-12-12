import ruamel.yaml as yaml
import pandas as pd
import argparse
from pathlib import Path
from sklearn.cluster import KMeans


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

if params['cluster_files']:
    cluspath = params['cluster_files']
    if cluspath[-1] not in ['/', '\\']:
        cluspath += '/'
else:
    cluspath = outpath

expressionfile = params['expressionfile']


if __name__ == '__main__':
    description = 'Concatenates all reaction-enumeration solutions'
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    args = parser.parse_args()

    prefix = 'rxn_enum_'
    if params['rxn_enum_params']['full_rxn_enum']:
        prefix += 'full_'


    genes = pd.read_csv(expressionfile, sep=';|,|\t', engine='python').set_index(params['gene_ID_column'])
    genes = genes.loc[genes.index.dropna()]
    if params['gene_expression_columns']:
        gene_conditions = [x.strip() for x in params['gene_expression_columns'].split(',')]
    else:
        gene_conditions = genes.columns.to_list()
    for condition in gene_conditions:
        solutions = []
        solfiles = Path(cluspath).glob(prefix + 'solutions_%s_*.csv' % condition)
        for f in solfiles:
            sol = pd.read_csv(f, index_col=0)
            solutions.append(sol)
        rxn_sols = pd.concat(solutions, ignore_index=True).drop_duplicates()

        fluxes = []
        fluxfiles = Path(cluspath).glob(prefix + 'fluxes_%s_*.csv' % condition)
        for f in fluxfiles:
            fl = pd.read_csv(f, index_col=0)
            fluxes.append(fl)
        rxn_fluxes = pd.concat(fluxes, ignore_index=True).loc[rxn_sols.index]

        rxn_sols.reset_index(inplace=True, drop=True)
        rxn_fluxes.reset_index(inplace=True, drop=True)

        clustering = KMeans(n_clusters=params['batch_num']).fit(rxn_sols)  # form batch_num kmeans clusters
        clusterdf = pd.DataFrame(clustering.transform(rxn_sols))
        sol_index = clusterdf.idxmin().values.tolist()  # we take the solution closest to each cluster center
        first_pos = list(set(range(10)) - set(sol_index))
        sol_index = list(set(sol_index) - set(range(10)))
        rename_dic = {}
        for i, j in zip(first_pos, sol_index):
            rename_dic[i] = j
            rename_dic[j] = i
        new_sols = rxn_sols.rename(rename_dic).sort_index()  # exchange the first solutions with the center solutions
        new_fluxes = rxn_fluxes.rename(rename_dic).sort_index()
        new_sols.to_csv(outpath + 'all_' + prefix + 'solutions_%s.csv' % condition)
        new_fluxes.to_csv(outpath + 'all_' + prefix + 'fluxes_%s.csv' % condition)
        if params['rxn_enum_params']['full_rxn_enum']:
            with open(cluspath + 'fullrxnenumdone_%s.txt' %condition, 'w+') as file:
                file.write(condition+' done')
