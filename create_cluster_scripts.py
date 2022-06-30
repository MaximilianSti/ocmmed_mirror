import ruamel.yaml as yaml
import pandas as pd
from warnings import warn


# read configuration from YAML file
yaml_reader = yaml.YAML(typ='safe')
with open('parameters.yaml', 'r') as file:
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
                    '#SBATCH -t {t}\n#SBATCH -o imat_{s}_out.out\n#SBATCH -e imat_{s}_err.out\n'
                    ''.format(s=condition, c=clus['cores'], m=clus['memory'], py=clus['pythonpath'], t=clus['time']))
            f.write('cd $SLURM_SUBMIT_DIR\ncd ..\nmodule purge\nmodule load %s\nsource %s/bin/activate\nexport '
                    'PYTHONPATH=${PYTHONPATH}:"%s"\n' %
                    (clus['pythonpath'], clus['envpath'], clus['cplexpath']))
            f.write('python utilities_cluster/cluster_weights_imat.py -c {s}'
                    ''.format(s=condition))
    with open('cluster_script_1.sh', 'w+') as f:
        f.write('#!/bin/bash\n#SBATCH --mail-type=ALL\n#SBATCH -J script_1\n#SBATCH -o script_1_out.out\n#SBATCH '
                '-e script_1_err.out\ncd $SLURM_SUBMIT_DIR\nfor i in {l}\ndo\n    dos2unix {p}imat_"$i".sh'
                '\n    sbatch {p}imat_"$i".sh\ndone'
                ''.format(l='{0..%i}' % (len(gene_conditions)-1), p=cluspath))
    if clus['approach'] == 'grouped':
        for i, condition in enumerate(gene_conditions):
            for j in range(clus['batch_num']):
                rxn_range = str(clus['batch_rxn_sols']*j) + '_' + str(clus['batch_rxn_sols']*(j+1))
                dist_a = (1 - 1 / (clus['batch_num'] * 2 * (clus['batch_div_sols'] / 10))) ** j
                with open(cluspath + "batch_%i_%i.sh" % (i, j), "w+") as f:
                    f.write('#!/bin/bash\n#SBATCH -p workq\n#SBATCH --mail-type=ALL\n#SBATCH --mem={m}G\n#SBATCH -c {c}'
                            '\n#SBATCH -t {t}\n#SBATCH -J dexom_{i}_{j}\n#SBATCH -o dexout_{i}_{j}.out\n#SBATCH -e '
                            'dexerr_{i}_{j}.out\n'
                            ''.format(s=condition, c=clus['cores'], m=clus['memory'], t=clus['time'], i=i, j=j))
                    f.write('cd $SLURM_SUBMIT_DIR\ncd ..\nmodule purge\nmodule load %s\nsource %s/bin/activate\nexport '
                            'PYTHONPATH=${PYTHONPATH}:"%s"\n' %
                            (clus['pythonpath'], clus['envpath'], clus['cplexpath']))
                    f.write('python utilities_cluster/cluster_rxn_enum.py -c {c} -r {r} -p {p}\n'
                            ''.format(c=condition, r=rxn_range, p=j))
                    f.write('python utilities_cluster/cluster_div_enum.py -c {c} -d {d} -i {maxiter} -p {p}'
                            ''.format(c=condition, d=dist_a, p=j, maxiter=clus['batch_div_sols']))
            with open(cluspath + 'dexom_%i.sh' % i, 'w+') as f:
                f.write('#!/bin/bash\n#SBATCH --mail-type=ALL\n#SBATCH -J dexom_{i}\n#SBATCH -o dexom_{i}_out.out\n'
                        '#SBATCH -e dexom_{i}_err.out\ncd $SLURM_SUBMIT_DIR\nfor j in {l}\ndo'
                        '\n    dos2unix batch_{i}_"$j".sh\n    sbatch batch_{i}_"$j".sh\ndone'
                        ''.format(i=i, l='{0..%i}' % (clus['batch_num']-1)))
        with open('cluster_script_2.sh', 'w+') as f:
            f.write('#!/bin/bash\n#SBATCH --mail-type=ALL\n#SBATCH -J script_2\n#SBATCH -o script_2_out.out\n#SBATCH '
                    '-e script_2_err.out\ncd $SLURM_SUBMIT_DIR\nfor i in {l}\ndo'
                    '\n    dos2unix {p}dexom_"$i".sh\n    sbatch {p}dexom_"$i".sh\ndone'
                    ''.format(l='{0..%i}' % (len(gene_conditions)-1), p=cluspath))
        # write script 3
        with open('cluster_script_3.sh', 'w+') as f:
            f.write('#!/bin/bash\n#SBATCH --mail-type=ALL\n#SBATCH -J script_3\n#SBATCH -o script_3_out.out\n#SBATCH '
                    '-e script_3_err.out\n#SBATCH --mem={m}G\n#SBATCH -c {c}\n#SBATCH -t {t}\n'
                    ''.format(p=cluspath, c=clus['cores'], m=clus['memory'], t=clus['time']))
            f.write('cd $SLURM_SUBMIT_DIR\nmodule purge\nmodule load %s\nsource %s/bin/activate\nexport '
                    'PYTHONPATH=${PYTHONPATH}:"%s"' %
                    (clus['pythonpath'], clus['envpath'], clus['cplexpath']))
            f.write('\npython utilities_cluster/cluster_concat_rxn_solutions.py\n'
                    'python utilities_cluster/cluster_produce_network.py'
                    ''.format(p=cluspath, c=clus['cores'], m=clus['memory'], t=clus['time']))

    elif clus['approach'] == 'separate':
        for i, condition in enumerate(gene_conditions):
            for j in range(clus['batch_num']):
                rxn_range = str(clus['batch_rxn_sols']*j) + '_' + str(clus['batch_rxn_sols']*(j+1))
                with open(cluspath + "batch_r_%i_%i.sh" % (i, j), "w+") as f:
                    f.write('#!/bin/bash\n#SBATCH -p workq\n#SBATCH --mail-type=ALL\n#SBATCH --mem={m}G\n#SBATCH -c {c}'
                            '\n#SBATCH -t {t}\n#SBATCH -J rxnenum_{i}_{j}\n#SBATCH -o rxnenumout_{i}_{j}.out\n#SBATCH '
                            '-e rxnenumerr_{i}_{j}.out\n'
                            ''.format(s=condition, c=clus['cores'], m=clus['memory'], t=clus['time'], i=i, j=j))
                    f.write('cd $SLURM_SUBMIT_DIR\ncd ..\nmodule purge\nmodule load %s\nsource %s/bin/activate\nexport '
                            'PYTHONPATH=${PYTHONPATH}:"%s"\n' %
                            (clus['pythonpath'], clus['envpath'], clus['cplexpath']))
                    f.write('python utilities_cluster/cluster_rxn_enum.py -c {c} -r {r} -p {p}\n'
                            ''.format(c=condition, r=rxn_range, p=j))
                dist_a = (1 - 1 / (clus['batch_num'] * 2 * (clus['batch_div_sols'] / 10))) ** j
                with open(cluspath + "batch_d_%i_%i.sh" % (i, j), "w+") as f:
                    f.write('#!/bin/bash\n#SBATCH -p workq\n#SBATCH --mail-type=ALL\n#SBATCH --mem={m}G\n#SBATCH -c {c}'
                            '\n#SBATCH -t {t}\n#SBATCH -J divenum_{i}_{j}\n#SBATCH -o divenum_{i}_{j}.out\n#SBATCH -e '
                            'divenum_{i}_{j}.out\n'
                            ''.format(s=condition, c=clus['cores'], m=clus['memory'], t=clus['time'], i=i, j=j))
                    f.write('cd $SLURM_SUBMIT_DIR\ncd ..\nmodule purge\nmodule load %s\nsource %s/bin/activate\nexport '
                            'PYTHONPATH=${PYTHONPATH}:"%s"\n' %
                            (clus['pythonpath'], clus['envpath'], clus['cplexpath']))
                    f.write('python utilities_cluster/cluster_div_enum.py -c {c} -d {d} -i {maxiter} -p {p}'
                            ''.format(c=condition, d=dist_a, p=j, maxiter=clus['batch_div_sols']))
            with open(cluspath + 'rxnenum_%i.sh' % i, 'w+') as f:
                f.write('#!/bin/bash\n#SBATCH --mail-type=ALL\n#SBATCH -J rxnenum_{i}\n#SBATCH -o rxnenum_{i}_out.out'
                        '\n#SBATCH -e rxnenum_{i}_err.out\ncd $SLURM_SUBMIT_DIR\nfor j in {l}\ndo'
                        '\n    dos2unix batch_r_{i}_"$j".sh\n    sbatch batch_r_{i}_"$j".sh\ndone'
                        ''.format(i=i, l='{0..%i}' % (clus['batch_num']-1)))
            with open(cluspath + 'divenum_%i.sh' % i, 'w+') as f:
                f.write('#!/bin/bash\n#SBATCH --mail-type=ALL\n#SBATCH -J divenum_{i}\n#SBATCH -o divenum_{i}_out.out'
                        '\n#SBATCH -e divenum_{i}_err.out\ncd $SLURM_SUBMIT_DIR\nfor j in {l}\ndo'
                        '\n    dos2unix batch_d_{i}_"$j".sh\n    sbatch batch_d_{i}_"$j".sh\ndone'
                        ''.format(i=i, l='{0..%i}' % (clus['batch_num']-1)))
        with open('cluster_script_2.sh', 'w+') as f:
            f.write('#!/bin/bash\n#SBATCH --mail-type=ALL\n#SBATCH -J script_2\n#SBATCH -o script_2_out.out\n#SBATCH '
                    '-e script_2_err.out\ncd $SLURM_SUBMIT_DIR\nfor i in {l}\ndo'
                    '\n    dos2unix {p}rxnenum_"$i".sh\n    sbatch {p}rxnenum_"$i".sh\ndone'
                    ''.format(l='{0..%i}' % (len(gene_conditions)-1), p=cluspath))
        with open('cluster_script_3.sh', 'w+') as f:
            f.write('#!/bin/bash\n#SBATCH --mail-type=ALL\n#SBATCH -J script_3\n#SBATCH -o script_3_out.out\n#SBATCH '
                    '-e script_3_err.out\ncd $SLURM_SUBMIT_DIR\n'
                    'python utilities_cluster/cluster_concat_rxn_solutions.py\nfor i in {l}\ndo'
                    '\n    dos2unix {p}divenum_"$i".sh\n    sbatch {p}divenum_"$i".sh\ndone'
                    ''.format(l='{0..%i}' % (len(gene_conditions)-1), p=cluspath))
        # write script 4
        with open('cluster_script_4.sh', 'w+') as f:
            f.write('#!/bin/bash\n#SBATCH --mail-type=ALL\n#SBATCH -J script_4\n#SBATCH -o script_4_out.out\n#SBATCH '
                    '-e script_4_err.out\n#SBATCH --mem={m}G\n#SBATCH -c {c}\n#SBATCH -t {t}\n'
                    ''.format(p=cluspath, c=clus['cores'], m=clus['memory'], t=clus['time']))
            f.write('cd $SLURM_SUBMIT_DIR\nmodule purge\nmodule load %s\nsource %s/bin/activate\nexport '
                    'PYTHONPATH=${PYTHONPATH}:"%s"' %
                    (clus['pythonpath'], clus['envpath'], clus['cplexpath']))
            f.write('\npython utilities_cluster/cluster_produce_network.py'
                    ''.format(p=cluspath, c=clus['cores'], m=clus['memory'], t=clus['time']))
    else:
        warn("Could not identify approach")
