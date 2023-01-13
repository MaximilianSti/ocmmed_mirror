
def get_conditions():
    return {CONDITIONS}

def get_parallel():
    return list(range({BATCH_NUM}))

rule end:
    input:
        "{OUT_FOLDER}cellspecific_model.xml"

rule produce_cellspecific_model:
    input:
        "{OUT_FOLDER}activation_frequency_reactions.csv"
    output:
        "{OUT_FOLDER}cellspecific_model.xml"
    log:
        "logs/cluster_final_2.log"
    shell:
        "python utilities_cluster/cluster_final_2.py"

rule concat_div_sols:
    input:
        expand("{CLUSTER_FILES}div_enum_solutions_{condition}_{parallel}.csv", condition=get_conditions(), parallel=get_parallel())
    output:
        "{OUT_FOLDER}all_DEXOM_solutions.csv",
        "{OUT_FOLDER}activation_frequency_reactions.csv"
    log:
        "logs/cluster_final_1.log"
    shell:
        "python utilities_cluster/cluster_final_1.py"

rule div_enum:
    input:
        "{OUT_FOLDER}reaction_weights_{condition}.csv",
        "{CLUSTER_FILES}full_rxn_enum_solutions_{condition}.csv"
    output:
        "{CLUSTER_FILES}div_enum_stats_{condition}_{parallel}.csv",
        "{CLUSTER_FILES}div_enum_solutions_{condition}_{parallel}.csv"
    params: 
        dist_anneal = lambda w: (1 - 1 / ({BATCH_NUM} * 2 * ({DIV_SOLS} / 10))) ** int(w.parallel)
    log:
        "logs/rxn_enum_{condition}_{parallel}.log"
    shell:
        "python utilities_cluster/cluster_div_enum.py -c {wildcards.condition} -p {wildcards.parallel} -d {params.dist_anneal} -i {DIV_SOLS}"

rule concat_rxn_sols:
    input:
        expand("{CLUSTER_FILES}rxn_enum_solutions_{condition}_{parallel}.csv", condition=get_conditions(), parallel=get_parallel())
    output:
        expand("{CLUSTER_FILES}full_rxn_enum_solutions_{condition}.csv", condition=get_conditions())
    log:
        "logs/concat_rxn_enum.log"
    shell:
        "python utilities_cluster/cluster_concat_rxn_solutions.py"

rule rxn_enum:
    input:
        "{OUT_FOLDER}reaction_weights_{condition}.csv",
        "{OUT_FOLDER}imat_solution_{condition}.csv"
    output:
        "{CLUSTER_FILES}rxn_enum_solutions_{condition}_{parallel}.csv"
    params: 
        rxn_range = lambda w: str({RXN_SOLS}*int(w.parallel)) + '_' + str({RXN_SOLS}*(int(w.parallel)+1))
    log:
        "logs/rxn_enum_{condition}_{parallel}.log"
    shell:
        "python utilities_cluster/cluster_rxn_enum.py -c {wildcards.condition} -p {wildcards.parallel} -r {params.rxn_range}"

rule weights_imat:
    input: 
        "{MODEL}",
        "{EXPRESSIONFILE}"
    output:
        "{OUT_FOLDER}reaction_weights_{condition}.csv",
        "{OUT_FOLDER}imat_solution_{condition}.csv"
    log:
        "logs/weights_imat_{condition}.log"
    shell:
        "python utilities_cluster/cluster_weights_imat.py -c {wildcards.condition}"