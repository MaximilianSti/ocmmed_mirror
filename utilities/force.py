

def force_active_rxns(model, rxnlist, fluxvalue=0.1):
    for rid in rxnlist:
        reaction = model.reactions.get_by_id(rid)
        if reaction.lower_bound == 0 and reaction.upper_bound > 0:
            reaction.lower_bound = fluxvalue
        elif reaction.lower_bound < 0 and reaction.upper_bound == 0:
            reaction.upper_bound = -fluxvalue
        elif reaction.lower_bound < 0 and reaction.upper_bound > 0:
            y_pos = model.solver.interface.Variable("force_%s_pos" % rid, type="binary")
            y_neg = model.solver.interface.Variable("force_%s_neg" % rid, type="binary")
            pos_constraint = model.solver.interface.Constraint(
                reaction.flux_expression + y_pos * (reaction.lower_bound - fluxvalue),
                lb=reaction.lower_bound, name="force_%s_lower" % rid)
            neg_constraint = model.solver.interface.Constraint(
                reaction.flux_expression + y_neg * (reaction.upper_bound + fluxvalue),
                ub=reaction.upper_bound, name="force_%s_upper" % rid)
            model.solver.add(y_pos)
            model.solver.add(y_neg)
            model.solver.add(pos_constraint)
            model.solver.add(neg_constraint)
            tot_constraint = model.solver.interface.Constraint(y_pos + y_neg, lb=1., name="force_active_reactions")
            model.solver.add(tot_constraint)
        else:
            reaction.upper_bound = 1000.
            reaction.lower_bound = fluxvalue
            print("reaction %s is blocked in this model, upper bound was set to 1000" % rid)


def force_reaction_bounds(model, rxndict):
    for rid in rxndict.keys():
        reaction = model.reactions.get_by_id(rid)
        bounds = [float(x.strip()) for x in rxndict[rid].split(',')]
        reaction.bounds = bounds
