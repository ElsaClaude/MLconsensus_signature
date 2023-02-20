rule make_file_tree:
    # params:
    #     output_directory = config['output_directory']
    input:
        expand('{output_directory}', output_directory = config['output_directory'])
    output:
        expand('{output_directory}/finish.txt',output_directory = config['output_directory'])
    shell:
        ""