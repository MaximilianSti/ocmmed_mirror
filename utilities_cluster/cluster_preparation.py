import dexom_python
import ruamel.yaml as yaml
import random
import os

yaml_reader = yaml.YAML(typ='safe')
with open('parameters.yaml', 'r') as file:
    a = file.read()
params = yaml_reader.load(a)


if __name__ == '__main__':

    if params['output_path']:
        outpath = params['output_path']
        os.makedirs(outpath, exist_ok=True)
        if outpath[-1] not in ['/', '\\']:
            outpath += '/'
    else:
        outpath = ''

    cluspath = outpath[:-1] + 'clusterfiles/'
    os.makedirs(cluspath, exist_ok=True)


    if not params['rxn_enum_params']['reaction_list']:
        model = dexom_python.read_model(params['modelpath'])
        reactions = [r.id for r in model.reactions]
        random.shuffle(reactions)
        with open(outpath + 'reactions_shuffled.txt', 'w+') as file:
            file.write('\n'.join(reactions))

    yaml_writer = yaml.YAML()
    with open(outpath + 'parameters_used_for_run.yaml', 'w+') as file:
        yaml_writer.dump(params, file)
