rule training:
    params:
        output_directory = config['output_directory'],
        name = config['name']
    input:
        expand('{output_directory}/Stratified_sampling/Pipeline_StratSamp_OK.txt',output_directory = config['output_directory'])
        expand('{input_data}/{metadata_file}', input_data = config['input_data'], metadata_file = config['metadata_file']),
        expand('{input_data}/{data_file}', input_data = config['input_data'], data_file = config['data_file'])
    output:
        expand('{output_directory}/Training/Pipeline_LaunchTraining_OK.txt',output_directory = config['output_directory'])
    script:
        "../scripts/training.py"