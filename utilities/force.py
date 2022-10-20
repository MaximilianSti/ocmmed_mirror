from warnings import warn
from pathlib import Path


def force_active_rxns(model, rxns, fluxvalue=0.1):
    if isinstance(rxns, str):
        rxnpath = Path(rxns)
        seps = ['\t', ';', ',', '\n']  # list of potential separators for the file
        if rxnpath.is_file():
            with open(rxnpath, 'r') as file:
                rxns = file.read()
        for sep in seps:
            rxns = rxns.replace(sep, ' ')
        rxnlist = rxns.split()
    elif isinstance(rxns, list):
        rxnlist = rxns
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
            if reaction.lower_bound == 0 and reaction.upper_bound > 0:
                reaction.lower_bound = fluxvalue
            elif reaction.lower_bound < 0 and reaction.upper_bound == 0:
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
                warn('reaction %s is blocked in this model, flux bounds must be changed before activity can be forced'
                     % rid)


def force_reaction_bounds(model, rxnbounds):
    if isinstance(rxnbounds, dict):
        rxndict = rxnbounds
    elif isinstance(rxnbounds, str):
        rxndict = {}
        with open(rxnbounds, 'r') as file:
            reader = file.read()
        rows = reader.split('\n')
        seps = ['\t', ';', ',']  # list of potential separators for the file
        for row in rows:
            for sep in seps:
                row = row.replace(sep, ' ')
            temp = row.split()
            rxndict[temp[0]] = ','.join(temp[1:])
    else:
        raise TypeError('Wrong type of parameter for force_flux_bounds. Expected either dict or str, got %s instead.'
                        % type(rxnbounds))
    for rid in rxndict.keys():
        try:
            reaction = model.reactions.get_by_id(rid)
        except KeyError:
            warn('KeyError in force_reaction_bounds: reaction %s was not found in the model, this reaction will be '
                 'skipped.' % rid)
        else:
            bounds = [float(x.strip()) for x in rxndict[rid].split(',')]
            reaction.bounds = bounds
