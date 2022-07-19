import pandas as pd
import miom
from cobra.io import write_sbml_model
from dexom_python import read_model, check_model_options


def mba(model_keep, enum_solutions, essential_reactions):

    enum_solutions.columns = [r.id for r in model_keep.reactions]
    frequencies = enum_solutions.sum().sort_values(ascending=False)
    vals = list(set(frequencies.values))
    vals.sort(reverse=True)

    allrecs = [r.id for r in model_keep.reactions]
    current_rxns = essential_reactions.copy()
    indexes = []
    newvals = vals

    if len(vals) >= 20:
        for i in range(1, 11):
            print("i", i)
            temp_model = model_keep.copy()
            index_value = i*int(len(vals)/10) if i != 10 else -1
            rxn_ids = frequencies[frequencies >= vals[index_value]].index.to_list()
            rxn_ids = list(set(rxn_ids) - set(current_rxns))
            current_rxns.extend(rxn_ids)
            rem_rxns = [temp_model.reactions.get_by_id(r) for r in list(set(allrecs) - set(current_rxns))]
            temp_model.remove_reactions(rem_rxns)

            miom_model = miom.load(miom.mio.cobra_to_miom(temp_model), "cplex")
            miom_model.steady_state().subset_selection(1).solve()
            flux_consistent_model = miom_model.select_subnetwork()

            if not False in [r in flux_consistent_model.network.R["id"] for r in essential_reactions]:
                print("tenth-step ends at iteration:", i)
                indexes.append(i)
                for rid in rxn_ids:
                    current_rxns.remove(rid)
                break
            newvals = vals[index_value:]

    if len(vals) >= 200:
        for j in range(1, 11):
            print("j", j)
            temp_model = model_keep.copy()
            index_value = (i-1)*int(len(vals)/10) + j*int(len(vals)/100) if j != 10 else i*int(len(vals)/10)
            if index_value >= len(vals):
                index_value = -1
            rxn_ids = frequencies[frequencies >= vals[index_value]].index.to_list()
            rxn_ids = list(set(rxn_ids) - set(current_rxns))
            current_rxns.extend(rxn_ids)
            rem_rxns = [temp_model.reactions.get_by_id(r) for r in list(set(allrecs) - set(current_rxns))]
            temp_model.remove_reactions(rem_rxns)

            miom_model = miom.load(miom.mio.cobra_to_miom(temp_model), "cplex")
            miom_model.steady_state().subset_selection(1).solve()
            flux_consistent_model = miom_model.select_subnetwork()

            if not False in [r in flux_consistent_model.network.R["id"] for r in essential_reactions]:
                print("hundredth-step ends at iteration:", j)
                indexes.append(j)
                for rid in rxn_ids:
                    current_rxns.remove(rid)
                break

        index_value = (i-1)*int(len(vals)/10) + (j-1)*int(len(vals)/100)
        newvals = vals[index_value:]

    for idx, v in enumerate(newvals[:]):
        print("v", v)
        temp_model = model_keep.copy()
        rxn_ids = frequencies[frequencies >= v].index.to_list()
        rxn_ids = list(set(rxn_ids) - set(current_rxns))
        current_rxns.extend(rxn_ids)
        rem_rxns = [temp_model.reactions.get_by_id(r) for r in list(set(allrecs) - set(current_rxns))]
        temp_model.remove_reactions(rem_rxns)

        miom_model = miom.load(miom.mio.cobra_to_miom(temp_model), "cplex")
        miom_model.steady_state().subset_selection(1).solve()
        flux_consistent_model = miom_model.select_subnetwork()

        if not False in [r in flux_consistent_model.network.R["id"] for r in essential_reactions]:
            print("reached frequency value:", v)
            break

    keep_recs = list(flux_consistent_model.network.R["id"])
    new_model = model_keep.copy()
    rem_recs = list(set([r.id for r in model_keep.reactions]) - set(keep_recs))
    new_model.remove_reactions(rem_recs, remove_orphans=True)

    new_model.id += "_mba"

    return new_model


if __name__ == "__main__":
    m = read_model("pilot_data/Human-GEM11.xml")
    check_model_options(m, feasibility=1e-6)
    sols = pd.read_csv("../pilot_data/dexom_solutions.csv", index_col=0)

    rxns = ["MAR09931"]
    m.reactions.get_by_id("MAR09931").upper_bound = 1000.

    new_model = mba(model_keep=m, enum_solutions=sols, essential_reactions=rxns)
    write_sbml_model(new_model, "../new_mba_model.xml")
