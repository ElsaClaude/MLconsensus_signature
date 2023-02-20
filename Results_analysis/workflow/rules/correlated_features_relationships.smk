rule correlated_features_relationships:
    params:
        prefix = config['prefix'],
        inputpath = config['input_path'],
        runname = config['run_name'],
        avgmcc = config['avg_mcc_threshold'],
        stdmcc = config['std_mcc']  
    input:
        expand('{inputpath}{prefix}/INFOGAIN/INFOGAIN.csv',
        inputpath = config['input_path'],
        prefix = config['prefix']),

        expand('{inputpath}{prefix}/PEARSON/pearson_corr_CONCAT.csv',
        inputpath = config['input_path'],
        prefix = config['prefix']),

        expand('{inputpath}{prefix}/SPEARMAN/spearman_corr_CONCAT.csv',
        inputpath = config['input_path'],
        prefix = config['prefix']),

        expand('{inputpath}{prefix}/NEO4J/{runname}_{prefix}_{suffix}.csv',
        inputpath = config['input_path'],
        prefix = config['prefix'],
        runname = config['run_name'],
        # suffix=["METHODS", "METHODS_FILTERED_MCC07_STD01","MODELS","MODELS_FILTERED_MCC07_STD01","FEATURES","FEATURES_FILTERED_MCC07_STD01","RELATIONSHIPS","RELATIONSHIPS_FILTERED_MCC07_STD01"])
        suffix=["RELATIONSHIPS","RELATIONSHIPS_FILTERED_MCC"+str(config['avg_mcc_threshold']).replace('.','')+"_STD"+str(config['std_mcc']).replace('.','')])

    output:
        expand('{inputpath}{prefix}/NEO4J/{runname}_{prefix}_{suffix}.csv',
        inputpath = config['input_path'],
        prefix = config['prefix'],
        runname = config['run_name'],
        # suffix=["METHODS", "METHODS_FILTERED_MCC07_STD01","MODELS","MODELS_FILTERED_MCC07_STD01","FEATURES","FEATURES_FILTERED_MCC07_STD01","RELATIONSHIPS","RELATIONSHIPS_FILTERED_MCC07_STD01"])
        suffix=["RELCORRFEAT","RELCORRFEAT_FILTERED_MCC"+str(config['avg_mcc_threshold']).replace('.','')+"_STD"+str(config['std_mcc']).replace('.','')])
    script:
        "../scripts/correlated_features_relationships.py"