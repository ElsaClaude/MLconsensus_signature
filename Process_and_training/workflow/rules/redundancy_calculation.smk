rule redundancy_calculation:
    params:
        input_path = config['input_path'],
        nbrun = config['prefix']
    input:
        expand('{inputpath}{prefix}/{runname}_{prefix}_{classifier}/{runname}_{prefix}_{classifier}_a.classification.data_to_train.csv',
        inputpath = config['input_path'],
        prefix = config['prefix'],
        runname = config['run_name'],
        classifier = config['classifiers'][0]),
        expand('{inputpath}{prefix}/{runname}_{prefix}_{classifiers}/{runname}_{prefix}_{classifiers}_c.classification.results.STANDARDIZED.csv',
        inputpath = config['input_path'],
        prefix = config['prefix'],
        runname = config['run_name'],
        classifiers = config['classifiers'])
    output:
        expand('{inputpath}{prefix}/INFOGAIN/INFOGAIN.csv',
        inputpath = config['input_path'],
        prefix = config['prefix']),
        expand('{inputpath}{prefix}/PEARSON/pearson_corr_CONCAT.csv',
        inputpath = config['input_path'],
        prefix = config['prefix']),
        expand('{inputpath}{prefix}/SPEARMAN/spearman_corr_CONCAT.csv',
        inputpath = config['input_path'],
        prefix = config['prefix'])
    script:
        "../scripts/redundancy_calculation.py"