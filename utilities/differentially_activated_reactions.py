import ruamel.yaml as yaml
import pandas as pd
import os

yaml_reader = yaml.YAML(typ='safe')
with open('parameters.yaml', 'r') as file:
    a = file.read()
doc = yaml_reader.load(a)
if doc['output_path']:
    outpath = doc['output_path']
    if outpath[-1] not in ['/', '\\']:
        outpath += '/'
else:
    outpath = ''


def compute_differentially_activated_reactions(input_folder=None, control=None, r2_threshold=0.2):
    """
    Computes Differentially Activated Reactions, as defined in https://www.biorxiv.org/content/10.1101/2023.06.30.547200v1.full
    The control condition is used as a reference for computing DARs. Defaults to the first entry in the conditions list
    """
    if doc['gene_expression_columns']:
        conditions = [x.strip() for x in doc['gene_expression_columns'].split(',')]
    else:
        genes = pd.read_csv(doc['expressionfile'], sep=';|,|\t', engine='python').set_index(doc['gene_ID_column'])
        conditions = genes.columns.to_list()
    if control is None:
        control = conditions[0]
    elif control not in conditions:
        raise ValueError('No condition named %s in condition list' % control)
    if input_folder is None:
        inputpath = outpath
    dataframes = {}
    for c in conditions:
        path = inputpath + 'all_DEXOM_solutions_%s.csv' % c
        if c == control:
            ctrl = pd.read_csv(path, index_col=0)
        else:
            df = pd.read_csv(path, index_col=0)
            dataframes[c] = df
    fctrl = ctrl.sum() / len(ctrl)
    darnumbers = pd.Series(index=conditions, dtype=int)

    outpath_DAR = outpath + 'DAR_analysis/'
    os.makedirs(outpath_DAR, exist_ok=True)
    for c in conditions:
        if c != control:
            fc = dataframes[c].sum() / len(dataframes[c])
            R2c = (fctrl - fc) ** 2
            DARc = fc[R2c > r2_threshold] - fctrl[R2c > r2_threshold]
            DARc.to_csv(outpath_DAR + 'DAR_significant_' + c + '.csv')
            R2c.to_csv(outpath_DAR + 'DAR_R2statistic_' + c + '.csv')
            (fc - fctrl).to_csv(outpath_DAR + 'DAR_frequency_diff_' + c + '.csv')
            darnumbers[c] = len(DARc)
    darnumbers.to_csv(outpath_DAR + 'DAR_allconditions.csv')
