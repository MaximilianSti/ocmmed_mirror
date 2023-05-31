import ruamel.yaml as yaml
import pandas as pd
import argparse
from pathlib import Path
from sklearn.cluster import KMeans


yaml_reader = yaml.YAML(typ='safe')
with open('parameters.yaml', 'r') as file:
    a = file.read()
doc = yaml_reader.load(a)

with open('params_cluster.yaml', 'r') as file:
    c = file.read()
clus = yaml_reader.load(c)

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

expressionfile = doc['expressionfile']


if __name__ == '__main__':
    description = 'Concatenates all reaction-enumeration solutions'
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    args = parser.parse_args()

    genes = pd.read_csv(expressionfile, sep=';|,|\t', engine='python').set_index(doc['gene_ID_column'])
    genes = genes.loc[genes.index.dropna()]
    if doc['gene_expression_columns']:
        gene_conditions = [x.strip() for x in doc['gene_expression_columns'].split(',')]
    else:
        gene_conditions = genes.columns.to_list()
    for condition in gene_conditions:
        solutions = []
        solfiles = Path(cluspath).glob('rxn_enum_solutions_%s_*.csv' % condition)
        for f in solfiles:
            sol = pd.read_csv(f, index_col=0)
            solutions.append(sol)
        rxn_sols = pd.concat(solutions, ignore_index=True).drop_duplicates()

        fluxes = []
        fluxfiles = Path(cluspath).glob('rxn_enum_fluxes_%s_*.csv' % condition)
        for f in fluxfiles:
            fl = pd.read_csv(f, index_col=0)
            fluxes.append(fl)
        rxn_fluxes = pd.concat(fluxes, ignore_index=True).loc[rxn_sols.index]

        rxn_sols.reset_index(inplace=True, drop=True)
        rxn_fluxes.reset_index(inplace=True, drop=True)

        clustering = KMeans(n_clusters=clus['batch_num']).fit(rxn_sols)  # form batch_num kmeans clusters
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
        new_sols.to_csv(cluspath + 'full_rxn_enum_solutions_%s.csv' % condition)
        new_fluxes.to_csv(cluspath + 'full_rxn_enum_fluxes_%s.csv' % condition)
