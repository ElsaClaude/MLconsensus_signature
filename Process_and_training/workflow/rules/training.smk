rule training:
    params:
        output_directory = config['output_directory'],
        name = config['name'],
        nb_samplings = config['nb_samplings'],
        cpus = config['cpus'],
        memory = config['memory']
    input:
        expand('{output_directory}/Stratified_sampling_{name}/Pipeline_StratSamp_OK.txt',output_directory = config['output_directory'], name = config['name'])
    output:
        expand('{output_directory}/Training_{name}/Pipeline_LaunchTraining_OK.txt',output_directory = config['output_directory'], name = config['name'])
    script:
        "../scripts/training.py"