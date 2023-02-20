rule models_features_relationships_neo4j:
    params:
        prefix = config['prefix'],
        inputpath = config['input_path'],
        nbrun = config['prefix'],
        nbclassif = config['classifiers'],
        nbclassifoptions = config['classif_options'],
        genomerefdir = config['genome_ref_dir'],
        runname = config['run_name'],
        avgmcc = config['avg_mcc_threshold'],
        stdmcc = config['std_mcc']
    input:
        expand('{inputpath}{prefix}/{runname}_{prefix}_{classifiers}/{runname}_{prefix}_{classifiers}_c.classification.results.STANDARDIZED.csv',
        inputpath = config['input_path'],
        prefix = config['prefix'],
        runname = config['run_name'],
        classifiers = config['classifiers']),

        expand('{inputpath}{prefix}/{runname}_{prefix}_{classifier}/{runname}_{prefix}_{classifier}_b.featureSelection.infoGain.csv',
        inputpath = config['input_path'],
        prefix = config['prefix'],
        runname = config['run_name'],
        classifier = config['classifiers'][0]),

        expand('{inputpath}{prefix}/INFOGAIN/INFOGAIN.csv',
        inputpath = config['input_path'],
        prefix = config['prefix'])

    output:
        expand('{inputpath}{prefix}/NEO4J/{runname}_{prefix}_{suffix}.csv',
        inputpath = config['input_path'],
        prefix = config['prefix'],
        runname = config['run_name'],
        suffix=["METHODS", "METHODS_FILTERED_MCC"+str(config['avg_mcc_threshold']).replace('.','')+"_STD"+str(config['std_mcc']).replace('.',''),"MODELS","MODELS_FILTERED_MCC"+str(config['avg_mcc_threshold']).replace('.','')+"_STD"+str(config['std_mcc']).replace('.',''),"FEATURES","FEATURES_FILTERED_MCC"+str(config['avg_mcc_threshold']).replace('.','')+"_STD"+str(config['std_mcc']).replace('.',''),"RELATIONSHIPS","RELATIONSHIPS_FILTERED_MCC"+str(config['avg_mcc_threshold']).replace('.','')+"_STD"+str(config['std_mcc']).replace('.','')])

    script:
        "../scripts/models_neo4j.py"