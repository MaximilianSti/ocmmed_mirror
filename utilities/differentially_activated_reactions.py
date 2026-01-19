import pandas as pd
import os


def compute_differentially_activated_reactions(inputpath, conditions, control=None, r2_threshold=0.2):
    """
    Computes Differentially Activated Reactions, as defined in https://doi.org/10.1186/s12859-024-05845-z
    A control condition can be used as a reference for computing DARs. Defaults to None, in which case all pairwise comparisons are performed.
    """
    if control is None:
        print('No control condition specified, all pairwise comparisons will be performed.')
    elif control not in conditions:
        raise ValueError('No condition named %s in condition list' % control)
    dataframes = {}
    for c in conditions:
        path = inputpath + 'all_DEXOM_solutions_%s.csv' % c
        df = pd.read_csv(path, index_col=0)
        dataframes[c] = df
    darnumbers = pd.DataFrame(index=conditions, columns=conditions, dtype=int)

    outpath_DAR = inputpath + 'DAR_analysis/'
    os.makedirs(outpath_DAR, exist_ok=True)
    if control is None:
        for c1 in conditions:
            f1 = dataframes[c1].sum() / len(dataframes[c1])
            for c2 in conditions:
                if c2 != c1:
                    fc2 = dataframes[c2].sum() / len(dataframes[c2])
                    R2c = (f1 - fc2) ** 2
                    DARc = fc2[R2c > r2_threshold] - f1[R2c > r2_threshold]
                    DARc.to_csv(outpath_DAR + 'DAR_significant_' + c2 + '_vs_' + c1 + '.csv')
                    R2c.to_csv(outpath_DAR + 'DAR_R2statistic_' + c2 + '_vs_' + c1 + '.csv')
                    (fc2 - f1).to_csv(outpath_DAR + 'DAR_frequency_diff_' + c2 + '_vs_' + c1 + '.csv')
                    darnumbers.loc[c2, c1] = len(DARc)
                else:
                    darnumbers.loc[c2, c1] = 0
    else:
        ctrl = dataframes[control]
        fctrl = ctrl.sum() / len(ctrl)
        for c in conditions:
            if c != control:
                fc = dataframes[c].sum() / len(dataframes[c])
                R2c = (fctrl - fc) ** 2
                DARc = fc[R2c > r2_threshold] - fctrl[R2c > r2_threshold]
                DARc.to_csv(outpath_DAR + 'DAR_significant_' + c + '.csv')
                R2c.to_csv(outpath_DAR + 'DAR_R2statistic_' + c + '.csv')
                (fc - fctrl).to_csv(outpath_DAR + 'DAR_frequency_diff_' + c + '.csv')
                darnumbers[ctrl, c] = len(DARc)
    darnumbers.astype(int).to_csv(outpath_DAR + 'DAR_allconditions.csv')
