rule stats:
    params:
        prefix = config['prefix'],
        outputpath = config['output_path'],
        runname = config['run_name'],
        avgmcc = config['avg_mcc_threshold'],
        stdmcc = config['std_mcc'],
        classifnames = config['classifiers_weka_name']  
    input:
        redundantfeat = expand('{inputpath}{prefix}/NEO4J/{runname}_{prefix}_{suffix}.csv',
        inputpath = config['input_path'],
        prefix = config['prefix'],
        runname = config['run_name'],
        suffix=["RELCORRFEAT","RELCORRFEAT_FILTERED_MCC"+str(config['avg_mcc_threshold']).replace('.','')+"_STD"+str(config['std_mcc']).replace('.','')]),

        infogainorder = expand('{inputpath}{prefix}/{classifier}/{runname}_{prefix}_{classifier}_b.featureSelection.infoGain.csv',
        inputpath = config['input_path'],
        prefix = config['prefix'],
        runname = config['run_name'],
        classifier = config['classifiers'][0]),

        resultsstand = expand('{inputpath}{prefix}/{classifiers}/{runname}_{prefix}_{classifiers}_c.classification.results.STANDARDIZED.csv',
        inputpath = config['input_path'],
        prefix = config['prefix'],
        runname = config['run_name'],
        classifiers = config['classifiers'])

    output:
        expand('{outputpath}finished.txt', outputpath = config['output_path'])
    script:
        "../scripts/stats.R"