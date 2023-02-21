rule training:
    params:
        output_directory = config['output_directory'],
        name = config['name'],
        nb_samplings = config['nb_samplings'],
        cpus = config['cpus'],
        memory = config['memory']
    input:
        expand('{output_directory}/Stratified_sampling/Pipeline_StratSamp_OK.txt',output_directory = config['output_directory'])
    output:
        expand('{output_directory}/Training/Pipeline_LaunchTraining_OK.txt',output_directory = config['output_directory'])
    script:
        "../scripts/training.py"