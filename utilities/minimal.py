import pandas as pd
import numpy as np
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

    # # rw = {r: 1 for r in essential_reactions}
    # # imatsol = imat(newmodel, rw)
    # # binary = (np.abs(imatsol.fluxes) >= 1e-5 - 1e-6).astype(int)
    # # print(binary)
    # # print(binary[essential_reactions])
    # # if binary[essential_reactions].sum() == len(essential_reactions):
    # #     indexes.append(idx)
    # #     print("reached frequency value:", v)
    # #     break

    # keep_recs = [r.id for r in newmodel.reactions]
    #
    # flux_consistent = False
    #
    # while not flux_consistent:
    #     newmodel = model.copy()
    #     index_value = i * int(len(vals) / 10) if i != 10 else -1
    #     rxn_ids = frequencies[frequencies >= vals[index_value]].index.to_list()
    #     rxn_ids = list(set(rxn_ids) - set(current_rxns))
    #     current_rxns.extend(rxn_ids)
    #     rem_rxns = [newmodel.reactions.get_by_id(r) for r in list(set(allrecs) - set(current_rxns))]
    #     newmodel.remove_reactions(rem_rxns)
    #
    #     temp_model = miom.load(miom.mio.cobra_to_miom(newmodel), "cplex")
    #     temp_model.steady_state().subset_selection(1).solve()
    #     test_model = temp_model.select_subnetwork()
    #     binary_temp = (np.abs(temp_model.variables.flux_values) >= 1e-5 - 1e-6).astype(int)
    #
    #     if False in [r in test_model.network.R["id"] for r in keep_recs]:
    #         break
    #     else:
    #         for rid in rxn_ids:
    #             current_rxns.remove(rid)
    #         print("flux-consistent in this i")
    #
    #     newmodel = model.copy()
    #     index_value = (i-1)*int(len(vals)/10) + j*int(len(vals)/100) if j != 10 else i*int(len(vals)/10)
    #     if index_value >= len(vals):
    #         index_value = -1
    #     rxn_ids = frequencies[frequencies >= vals[index_value]].index.to_list()
    #     rxn_ids = list(set(rxn_ids) - set(current_rxns))
    #     current_rxns.extend(rxn_ids)
    #     rem_rxns = [newmodel.reactions.get_by_id(r) for r in list(set(allrecs) - set(current_rxns))]
    #     newmodel.remove_reactions(rem_rxns)
    #
    #     temp_model = miom.load(miom.mio.cobra_to_miom(newmodel), "cplex")
    #     temp_model.steady_state().subset_selection(1).solve()
    #     test_model = temp_model.select_subnetwork()
    #     binary_temp = (np.abs(temp_model.variables.flux_values) >= 1e-5 - 1e-6).astype(int)
    #
    #     if False in [r in test_model.network.R["id"] for r in keep_recs]:
    #         break
    #     else:
    #         for rid in rxn_ids:
    #             current_rxns.remove(rid)
    #         print("flux-consistent in this j")
    #
    #     for idx, v in enumerate(newvals[idx:]):
    #         newmodel = model.copy()
    #         rxn_ids = frequencies[frequencies >= v].index.to_list()
    #         rxn_ids = list(set(rxn_ids) - set(current_rxns))
    #         current_rxns.extend(rxn_ids)
    #         rxns = [model.reactions.get_by_id(r) for r in rxn_ids]
    #         # newmodel.add_reactions(rxns)
    #         rem_rxns = [newmodel.reactions.get_by_id(r) for r in list(set(allrecs) - set(current_rxns))]
    #         newmodel.remove_reactions(rem_rxns)
    #
    #         temp_model = miom.load(miom.mio.cobra_to_miom(newmodel), "cplex")
    #         temp_model.steady_state().subset_selection(1).solve()
    #         test = temp_model.select_subnetwork()
    #         binary_temp = (np.abs(temp_model.variables.flux_values) >= 1e-5 - 1e-6).astype(int)
    #
    #         if not False in [r in test_model.network.R["id"] for r in keep_recs]:
    #             flux_consistent = True
    #             for rid in rxn_ids:
    #                 current_rxns.remove(rid)
    #             break
    #     print("this statement should not be reached (end of the while loop)")

    # if not flux_consistent:
    #     #continue the previous loop, but this time th eend condition is keep_recs
    #     for i in range(indexes[0], 11):
    #         print("i", i)
    #         newmodel = model.copy()
    #         index_value = i * int(len(vals) / 10) if i != 10 else -1
    #         rxn_ids = frequencies[frequencies >= vals[index_value]].index.to_list()
    #         rxn_ids = list(set(rxn_ids) - set(current_rxns))
    #         current_rxns.extend(rxn_ids)
    #         rem_rxns = [newmodel.reactions.get_by_id(r) for r in list(set(allrecs) - set(current_rxns))]
    #         newmodel.remove_reactions(rem_rxns)
    #
    #         temp_model = miom.load(miom.mio.cobra_to_miom(newmodel), "cplex")
    #         temp_model.steady_state().subset_selection(1).solve()
    #         test = temp_model.select_subnetwork()
    #         binary_temp = (np.abs(temp_model.variables.flux_values) >= 1e-5 - 1e-6).astype(int)
    #
    #         if not False in [r in test_model.network.R["id"] for r in keep_recs]:
    #             flux_consistent = True
    #             print("tenth-step ends at iteration:", i)
    #             for rid in rxn_ids:
    #                 current_rxns.remove(rid)
    #             break
    #
    #     for j in range(1, 11):
    #         print("j", j)
    #         newmodel = model.copy()
    #         index_value = (i - 1) * int(len(vals) / 10) + j * int(len(vals) / 100) if j != 10 else i * int(
    #             len(vals) / 10)
    #         if index_value >= len(vals):
    #             index_value = -1
    #         rxn_ids = frequencies[frequencies >= vals[index_value]].index.to_list()
    #         rxn_ids = list(set(rxn_ids) - set(current_rxns))
    #         current_rxns.extend(rxn_ids)
    #         rem_rxns = [newmodel.reactions.get_by_id(r) for r in list(set(allrecs) - set(current_rxns))]
    #         newmodel.remove_reactions(rem_rxns)
    #
    #         temp_model = miom.load(miom.mio.cobra_to_miom(newmodel), "cplex")
    #         temp_model.steady_state().subset_selection(1).solve()
    #         test = temp_model.select_subnetwork()
    #         binary_temp = (np.abs(temp_model.variables.flux_values) >= 1e-5 - 1e-6).astype(int)
    #
    #         if not False in [r in test_model.network.R["id"] for r in keep_recs]:
    #             flux_consistent = True
    #             print("hundredth-step ends at iteration:", j)
    #             for rid in rxn_ids:
    #                 current_rxns.remove(rid)
    #             break
    #
    #     index_value = (i - 1) * int(len(vals) / 10) + (j - 1) * int(len(vals) / 100)
    #     newvals = vals[index_value:]
    #
    #     for idx, v in enumerate(newvals):
    #         print("v", v)
    #         newmodel = model.copy()
    #         rxn_ids = frequencies[frequencies >= v].index.to_list()
    #         rxn_ids = list(set(rxn_ids) - set(current_rxns))
    #         current_rxns.extend(rxn_ids)
    #         rxns = [model.reactions.get_by_id(r) for r in rxn_ids]
    #         # newmodel.add_reactions(rxns)
    #         rem_rxns = [newmodel.reactions.get_by_id(r) for r in list(set(allrecs) - set(current_rxns))]
    #         newmodel.remove_reactions(rem_rxns)
    #
    #         temp_model = miom.load(miom.mio.cobra_to_miom(newmodel), "cplex")
    #         temp_model.steady_state().subset_selection(1).solve()
    #         test = temp_model.select_subnetwork()
    #         binary_temp = (np.abs(temp_model.variables.flux_values) >= 1e-5 - 1e-6).astype(int)
    #
    #         if not False in [r in test_model.network.R["id"] for r in keep_recs]:
    #             flux_consistent = True
    #             print("reached frequency value:", v)
    #             break