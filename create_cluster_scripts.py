import dexom_python
import ruamel.yaml as yaml
import pandas as pd
from utilities.force import force_active_rxns
import argparse
from warnings import warn


# read configuration from YAML file
yaml_reader = yaml.YAML(typ='safe')
with open('config.yaml', 'r') as file:
    a = file.read()
doc = yaml_reader.load(a)

with open('additional_params.yaml', 'r') as file:
    b = file.read()
params = yaml_reader.load(b)

with open('cluster_params.yaml', 'r') as file:
    c = file.read()
clus = yaml_reader.load(c)

if doc['output_path']:
    outpath = doc['output_path']
else:
    outpath = ''

if clus['cluster_files']:
    cluspath = clus['cluster_files']
else:
    cluspath = ''

expressionpath = doc['expressionpath']

if __name__ == '__main__':
    print('writing cluster scripts')

    if doc['gene_expression_columns']:
        gene_conditions = [x.strip() for x in doc['gene_expression_columns'].split(',')]
    else:
        genes = pd.read_csv(expressionpath).set_index(doc['gene_ID_column'])
        gene_conditions = genes.columns.to_list()

    for i, condition in enumerate(gene_conditions):
        with open(cluspath+'imat_%i.sh' % i, 'w+') as f:
            f.write('#!/bin/bash\n#SBATCH --mail-type=ALL\n#SBATCH -J imat_{s}\n#SBATCH -c {c}\n#SBATCH --mem={m}G\n'
                    '#SBATCH -t 1:00:00\n#SBATCH -o imat_{s}_out.out\n#SBATCH -e imat_{s}_err.out\ncd $SLURM_SUBMIT_DIR'
                    '\ncd ..\nmodule purge\nmodule load {py}\npython utilities_cluster/cluster_weights_imat.py -c {s}'
                    ''.format(s=condition, c=clus['cores'], m=clus['memory'], py=clus['pythonpath']))
    with open('cluster_script_1.sh', 'w+') as f:
        f.write('#!/bin/bash\n#SBATCH --mail-type=ALL\n#SBATCH -J script_1\n#SBATCH -o script_1_out.out\n#SBATCH '
                '-e script_1_err.out\ncd $SLURM_SUBMIT_DIR\ncd ..\nfor i in {0..{l}}\ndo\n    dos2unix {p}/imat_"$i".sh'
                '\n    sbatch {p}/imat_"$i".sh\ndone'
                ''.format(p=cluspath, l=len(gene_conditions)-1))

    if clus['approach'] == 'grouped':
        for i, condition in enumerate(gene_conditions):
            for j in range(clus['batch_num']):
                rxn_range = str(clus['batch_rxn_sols']*j) + '_' + str(clus['batch_rxn_sols']*(j+1))
                dist_a = (1 - 1 / (clus['batch_num'] * 2 * (clus['batch_div_sols'] / 10))) ** j
                with open(cluspath + "batch_%i_%i.sh" % (i, j), "w+") as f:
                    f.write('#!/bin/bash\n#SBATCH -p workq\n#SBATCH --mail-type=ALL\n#SBATCH --mem={c}G\n#SBATCH -c {c}'
                            '\n#SBATCH -t {t}\n#SBATCH -J dexom_{i}_{j}\n#SBATCH -o dexout_{i}_{j}.out\n#SBATCH -e '
                            'dexerr_{i}_{j}.out\n'
                            ''.format(s=condition, c=clus['cores'], m=clus['memory'], t=clus['tme'],i=i, j=j))
                    f.write('cd $SLURM_SUBMIT_DIR\ncd ..\nmodule purge\nmodule load %s\nsource %s/bin/activate\nexport '
                            'PYTHONPATH=${PYTHONPATH}:"%s"\n' %
                            (clus['pythonpath'], clus['envpath'], clus['cplexpath']))
                    f.write('python utilities_cluster/cluster_rxn_enum.py -c {c} -r {r} -p {p}\n'
                            ''.format(c=condition, r=rxn_range, p=j))
                    f.write('python utilities_cluster/cluster_div_enum.py -c {c} -d {d} -i {j}'
                            ''.format(c=condition, d=dist_a, j=j))
        #to do: write script_2 for launching all of this
        # write script 3
        pass
    elif clus['approach'] == 'separate':
        pass
    else:
        warn("Could not identify approach")