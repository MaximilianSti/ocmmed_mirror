from warnings import warn
import pandas as pd


def force_active_rxns(model, rxns, fluxvalue=0.1, condition=None):
    if isinstance(rxns, list):
        rxnlist = rxns
    elif isinstance(rxns, str):
        file = pd.read_csv(rxns, sep=';|,|\t', engine='python')
        if condition.isnumeric():
            try:
                file = file[file.columns[int(condition)]]
            except IndexError as e:
                raise Exception(f'Condition number {condition} is not present in the {rxns} file') from e
        else:
            try:
                file = file[condition]
            except KeyError as e:
                raise Exception(f'Condition {condition} is not present in the {rxns} file') from e
        rxnlist = list(file[file==file].values)
    else:
        raise TypeError('Wrong type for parameter force_active_reactions. Expected either list or str, got %s instead.'
                        % type(rxns))
    for rid in rxnlist:
        try:
            reaction = model.reactions.get_by_id(rid)
        except KeyError:
            warn('KeyError in force_active_reactions: reaction %s was not found in the model, this reaction will be '
                 'skipped.' % rid)
        else:
            if reaction.lower_bound >= 0 and reaction.upper_bound > 0:
                if reaction.lower_bound > 0:
                    warn(f'The positive lower bound of reaction {rid} was overwritten: '
                         f'{reaction.lower_bound} -> {fluxvalue}')
                reaction.lower_bound = fluxvalue
            elif reaction.lower_bound < 0 and reaction.upper_bound <= 0:
                if reaction.upper_bound < 0:
                    warn(f'The negative upper bound of reaction {rid} was overwritten: '
                         f'{reaction.upper_bound} -> {-fluxvalue}')
                reaction.upper_bound = -fluxvalue
            elif reaction.lower_bound < 0 and reaction.upper_bound > 0:
                y_pos = model.solver.interface.Variable('force_%s_pos' % rid, type='binary')
                y_neg = model.solver.interface.Variable('force_%s_neg' % rid, type='binary')
                pos_constraint = model.solver.interface.Constraint(
                    reaction.flux_expression + y_pos * (reaction.lower_bound - fluxvalue),
                    lb=reaction.lower_bound, name='force_%s_lower' % rid)
                neg_constraint = model.solver.interface.Constraint(
                    reaction.flux_expression + y_neg * (reaction.upper_bound + fluxvalue),
                    ub=reaction.upper_bound, name='force_%s_upper' % rid)
                model.solver.add(y_pos)
                model.solver.add(y_neg)
                model.solver.add(pos_constraint)
                model.solver.add(neg_constraint)
                tot_constraint = model.solver.interface.Constraint(y_pos + y_neg, lb=1., name='force_active_reactions')
                model.solver.add(tot_constraint)
            else:
                raise Exception(f'reaction {rid} is blocked in this model, '
                                f'flux bounds must be changed before activity can be forced')


def force_reaction_bounds(model, rxnbounds, condition=None):
    if isinstance(rxnbounds, dict):
        rxndict = rxnbounds
    elif isinstance(rxnbounds, str):
        file = pd.read_csv(rxnbounds, sep=';|,|\t', engine='python', index_col=0)
        if condition.isnumeric():
            try:
                index = int(condition)
                file = file[[file.columns[index*2], file.columns[index*2+1]]]
            except IndexError as e:
                raise Exception(f'Condition number {condition} is not present in the {rxnbounds} file') from e
        else:
            try:
                file = file[[condition+'_lb', condition+'_ub']]
            except KeyError as e:
                raise Exception(f'Condition {condition} is not present in the {rxnbounds} file') from e
        rxndict = {r: str(v[0]) + ',' + str(v[1]) for r, v in file.iterrows()}
    else:
        raise TypeError('Wrong type of parameter for force_flux_bounds. Expected either dict or str, got %s instead.'
                        % type(rxnbounds))
    for rid in rxndict.keys():
        try:
            reaction = model.reactions.get_by_id(rid)
        except KeyError as e :
            raise Exception(f'KeyError in force_reaction_bounds: reaction {rid} was not found in the model') from e
        else:
            bounds = [float(x.strip()) for x in rxndict[rid].split(',')]
            reaction.bounds = bounds
