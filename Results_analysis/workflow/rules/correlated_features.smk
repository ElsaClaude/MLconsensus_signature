rule correlated_features:
    params:
        prefix = config['prefix'],
        inputpath = config['input_path'],
        runname = config['run_name'],
        avgmcc = config['avg_mcc_threshold'],
        stdmcc = config['std_mcc'],
    input:
        expand('{inputpath}{prefix}/NEO4J/{runname}_{prefix}_{suffix}.csv',
        inputpath = config['input_path'],
        prefix = config['prefix'],
        runname = config['run_name'],
        suffix=["FEATURES","FEATURES_FILTERED_MCC"+str(config['avg_mcc_threshold']).replace('.','')+"_STD"+str(config['std_mcc']).replace('.',''),"RELCORRFEAT","RELCORRFEAT_FILTERED_MCC"+str(config['avg_mcc_threshold']).replace('.','')+"_STD"+str(config['std_mcc']).replace('.','')])

    output:
        expand('{inputpath}{prefix}/NEO4J/{runname}_{prefix}_{suffix}.csv',
        inputpath = config['input_path'],
        prefix = config['prefix'],
        runname = config['run_name'],
        suffix=["CORRFEAT","CORRFEAT_FILTERED_MCC"+str(config['avg_mcc_threshold']).replace('.','')+"_STD"+str(config['std_mcc']).replace('.','')])
    script:
        "../scripts/correlated_features.py"