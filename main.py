import dexom_python
import ruamel.yaml as yaml
import pandas as pd
from cobra.io import write_sbml_model
from cobra.flux_analysis import find_blocked_reactions
from utilities.force import force_active_rxns, force_reaction_bounds
from utilities.minimal import maximal_frequency
from utilities.inactive_pathways import compute_inactive_pathways
from utilities.differentially_activated_reactions import compute_differentially_activated_reactions
from warnings import warn, catch_warnings, filterwarnings, resetwarnings
import os


# read configuration from YAML files
yaml_reader = yaml.YAML(typ='safe')
with open('parameters.yaml', 'r') as file:
    a = file.read()
doc = yaml_reader.load(a)

with open('params_additional.yaml', 'r') as file:
    b = file.read()
params = yaml_reader.load(b)

if doc['output_path']:
    outpath = doc['output_path']
    os.makedirs(outpath, exist_ok=True)
    if outpath[-1] not in ['/', '\\']:
        outpath += '/'
else:
    outpath = ''

mp = params['model_params']
ip = params['imat_params']
ep = params['enum_params']
rp = params['rxn_enum_params']
dp = params['div_enum_params']

modelpath = doc['modelpath']
expressionfile = doc['expressionfile']


if __name__ == '__main__':
    # read model
    model_keep = dexom_python.read_model(modelpath, solver=mp['solver'])
    model_keep = dexom_python.check_model_options(model_keep, timelimit=mp['timelimit'], feasibility=mp['feasibility'],
                                                  mipgaptol=mp['mipgaptol'], verbosity=mp['verbosity'])
    # read and process gene expression file
    genes = pd.read_csv(expressionfile, sep=';|,|\t', engine='python').set_index(doc['gene_ID_column'])
    genes = genes.loc[genes.index.dropna()]
    if doc['gene_expression_columns']:
        gene_conditions = [x.strip() for x in doc['gene_expression_columns'].split(',')]
    else:
        gene_conditions = genes.columns.to_list()
    if doc['gpr_parameters']['qualitative'] and not doc['reaction_scores']:
        genes = dexom_python.expression2qualitative(genes=genes, column_list=gene_conditions,
                                                    proportion=doc['gpr_parameters']['percentile'],
                                                    outpath=outpath+'geneweights_qualitative')
    dexom_sols = []
    for condition in gene_conditions:
        new_model = model_keep.copy()
        if doc['force_flux_bounds']:
            force_reaction_bounds(new_model, doc['force_flux_bounds'], condition)
        if doc['force_active_reactions']:
            force_active_rxns(new_model, doc['force_active_reactions'], doc['fluxvalue'], condition)
        # create reaction weights from gene expression
        print('computing reaction weights for condition '+condition)
        gene_weights = pd.Series(genes[condition].values, index=genes.index, dtype=float)

        if doc['reaction_scores']:
            rw = {}
            for rxn in model_keep.reactions:
                rw[rxn.id] = float(gene_weights.to_dict().get(rxn.id, 0.))
            dexom_python.save_reaction_weights(rw, outpath+'reaction_weights_%s' % condition)
        else:
            rw = dexom_python.apply_gpr(model=model_keep, gene_weights=gene_weights, duplicates=doc['duplicates'],
                                        save=True, filename=outpath+'reaction_weights_%s' % condition)

        # compute imat solution from reaction weights
        print('performing iMAT for condition ' + condition)
        model = new_model.copy()
        model = dexom_python.check_model_options(model, timelimit=mp['timelimit'], feasibility=mp['feasibility'],
                                                 mipgaptol=mp['mipgaptol'], verbosity=mp['verbosity'])
        with catch_warnings():
            filterwarnings('error')
            try:
                imatsol = dexom_python.imat(model=model, reaction_weights=rw, epsilon=ip['epsilon'],
                                            threshold=ip['threshold'])
            except UserWarning:
                resetwarnings()
                warn('Solver could not find a solution with forced active reactions. '
                      'Attempting to find a solution without forced flux.')
                model = model_keep.copy()
                model = dexom_python.check_model_options(model, timelimit=mp['timelimit'],
                                                         feasibility=mp['feasibility'],
                                                         mipgaptol=mp['mipgaptol'], verbosity=mp['verbosity'])
                imatsol = dexom_python.imat(model=model, reaction_weights=rw, epsilon=ip['epsilon'],
                                            threshold=ip['threshold'])
                print('The solver could not find a solution with forced active reactions, '
                      'but found one without forcing flux through the selected reactions. '
                      'Check if any of the reactions in the "active_reactions" list are blocked ')
        dexom_python.write_solution(model=model, solution=imatsol, threshold=ip['threshold'],
                                    filename=outpath+'imat_solution_%s.csv' % condition)

        # enumerate solutions with the DEXOM method
        print('performing reaction-enumeration for condition ' + condition)
        if rp['reaction_list']:
            df = pd.read_csv(rp['reaction_list'], header=None)
            reactions = [x for x in df.unstack().values]
            wrongrids = [rid for rid in reactions if rid not in [r.id for r in model.reactions]]
            for rid in wrongrids:
                warn('reaction %s is not in the model, this reaction will be skipped' % rid)
            reactions = list(set(reactions) - set(wrongrids))
        else:
            reactions = [r.id for r in model.reactions]
        if params['blocked_rxns']:
            seps = ['\t', ';', ',', '\n']  # list of potential separators for the file
            with open(params['blocked_rxns']) as file:
                reader = file.read()
            for sep in seps:
                reader = reader.replace(sep, ' ')
            rxns_inactive = reader.split()
            reactions = list(set(reactions) - set(rxns_inactive))
        rxnlist = reactions[:doc['rxn_enum_iterations']]
        rxnsol = dexom_python.enum_functions.rxn_enum(model=model, reaction_weights=rw, prev_sol=imatsol,
                                                      rxn_list=rxnlist,obj_tol=ep['objective_tolerance'],
                                                      eps=ip['epsilon'], thr=ip['threshold'])
        uniques = pd.DataFrame(rxnsol.unique_binary)
        uniques.columns = [r.id for r in model.reactions]
        uniques.to_csv(outpath+'rxn_enum_solutions_%s.csv' % condition)
        fluxes = pd.concat([s.fluxes for s in rxnsol.unique_solutions], axis=1).T.reset_index().drop('index', axis=1)
        fluxes.to_csv(outpath+'rxn_enum_fluxes_%s.csv' % condition)
        print("length reaction-enumeration solution:", len(uniques))
        print('performing diversity-enumeration for condition ' + condition)
        divsol, divres = dexom_python.enum_functions.diversity_enum(model=model, reaction_weights=rw, prev_sol=imatsol,
                                dist_anneal=dp['dist_anneal'], eps=ip['epsilon'], thr=ip['threshold'], icut=dp['icut'],
                                maxiter=doc['div_enum_iterations'], obj_tol=ep['objective_tolerance'],  full=dp['full'])
        divs = pd.DataFrame(divsol.binary)
        divs.columns = [r.id for r in model.reactions]
        divs.to_csv(outpath+'div_enum_solutions_%s.csv' % condition)
        fluxd = pd.concat([s.fluxes for s in divsol.solutions], axis=1).T.reset_index().drop('index', axis=1)
        fluxd.to_csv(outpath+'div_enum_fluxes_%s.csv' % condition)
        dexomsol = pd.concat([uniques, divs]).drop_duplicates(ignore_index=True)
        dexomsol.to_csv(outpath+'all_DEXOM_solutions_%s.csv' % condition)
        dexom_sols.append(dexomsol)

    dexom_sols = pd.concat(dexom_sols).drop_duplicates(ignore_index=True)
    dexom_sols.to_csv(outpath + 'dexom_solutions.csv')
    dexom_sols.columns = [r.id for r in model_keep.reactions]
    frequencies = dexom_sols.sum()._set_name('frequency')
    frequencies.to_csv(outpath + 'activation_frequency_reactions.csv')

    print('producing final network')
    if doc['final_network'] == 'union':
        cutoff = doc['union_cutoff']
        if isinstance(cutoff, str):
            if cutoff[-1] == '%':
                cutoff = frequencies.max() * float(cutoff[:-1]) / 100
            else:
                warn('Unrecognized character in union_cutoff parameter, default to 0.')
                cutoff = 0
        rem_rxns = dexom_sols.columns[frequencies <= cutoff].to_list()  # remove reactions which are active in zero solutions
        new_model.remove_reactions(rem_rxns, remove_orphans=True)
        if cutoff > 0:
            blocked_reactions = find_blocked_reactions(new_model)
            new_model.remove_reactions(blocked_reactions, remove_orphans=True)
    elif doc['final_network'] == 'minimal':
        new_model = maximal_frequency(model_keep=new_model, frequency_table=frequencies, essential_reactions=doc['force_active_reactions'])
    else:
        warn('Invalid value for "final_network" in parameters.yaml, returning original network.')
    new_model.id += '_cellspecific'
    write_sbml_model(new_model, outpath+'cellspecific_model.xml')
    compute_inactive_pathways(new_model)
    compute_differentially_activated_reactions()
