rule stratified_sampling:
    params:
        output_directory = config['output_directory'],
        name = config['name'],
        nb_samplings = config['nb_samplings']
    input:
        expand('{input_data}/{metadata_file}', input_data = config['input_data'], metadata_file = config['metadata_file']),
        expand('{input_data}/{data_file}', input_data = config['input_data'], data_file = config['data_file'])
    output:
        expand('{output_directory}/Stratified_sampling/Pipeline_StratSamp_OK.txt',output_directory = config['output_directory'])
    script:
        "../scripts/stratified_sampling.py"