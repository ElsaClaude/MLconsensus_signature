# MaLCons
 MaLCons stands for "Machine Learning Consensus" and is a two steps pipeline designed to perform the building of a consensus between feature signature of various Machine Learning (ML) models. The pipeline has been developed based on Snakemake workflow management system.
 
 

## Pipeline global behavior

 1. **run Sub-pipeline** --> Process and training
    - Stratified Sampling: execution of a data stratified sampling step for machine learning training
    - Training: launch ML models training
 2. **run script** --> check_training_process.py: check training status
 3. **run script** --> create_config_Results_analysis: automatic set up of the second sub-pipeline config.yaml file
 4. **run Sub-pipeline** --> Results analysis
    - Standardize results: needed step to standardize results of various ML model training runs
    - Redundancy calculation: compute features identified as bearing redundant information during training (for each sampling run)
    - Models features relationships: retrieve associations between various ML models and the features name composing their signatures
    - Correlated features relationship: 
    - Correlated features: 
    - Statistics: Compute exploration of the global ML training through the various models, methods and runs.\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**Give a proposed consensus feature signature as final result**
 
 
## Pipeline set up
### 1. Sub-pipeline --> Process and training
Adjust the `config.yaml` file of the **Process and training** and **Results analysis** step in the config folder.

Help for **Process and training** `config.yaml` file:

```
#input data folder path
input_data: /home/elsa/projects/rrg-adroit-ab/elsa/Thesis/Data/data_test_paper_pipeline

#input metadata file name
#metadata file must contain an ID column where IDs are the same as in the full data file
#metadata file must contain a Label and a Cluster column
metadata_file: CRC_METADATA.csv 

#input full data file name
#full data file must contain the same ID column as in the metadata file
#features and associated values can not contain missing values
data_file: CRC_gene_data.csv 

#output folder path
output_directory: /home/elsa/projects/rrg-adroit-ab/elsa/Thesis/test_paper_pipeline/output

#name of your project
name: TestPipelinePaper2 

#number of wanted run
nb_samplings: 3 

#number of available cpus
cpus: 2 

#amount of available memory
memory: 20G 
```
**NB:**
- To chose which ML models you want to train, you must update the files in `MaLCons/Process_and_training/softs/classifiers_list/`.
- The training script build SLURM sbatch files to be able to run training jobs. It has been designed to run on Digital Research Alliance of Canada (the Alliance). You can adapt the sbatch file format by modifying the `MaLCons/Process_and_training/workflow/scripts/training.py`

