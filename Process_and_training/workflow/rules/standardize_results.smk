rule standardize:
    params:
        prefix = config['prefix'],
        inputpath = config['input_path'],
        outputpath = config['output_path']
    input:
        expand('{inputpath}{prefix}/{runname}_{prefix}_{classifiers}/{runname}_{prefix}_{classifiers}_c.classification.results.csv',
        inputpath = config['input_path'],
        prefix = config['prefix'],
        runname = config['run_name'],
        classifiers = config['classifiers'])
    output:
        expand('{inputpath}{prefix}/{runname}_{prefix}_{classifiers}/{runname}_{prefix}_{classifiers}_c.classification.results.STANDARDIZED.csv',
        inputpath = config['input_path'],
        prefix = config['prefix'],
        runname = config['run_name'],
        classifiers = config['classifiers'])
    script:
        "../scripts/standardize_results.py"