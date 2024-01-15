# MaLCons
 MaLCons stands for "Machine Learning Consensus" and is a two steps pipeline designed to perform the building of a stable between feature signatures of multiple wrappers Machine Learning (ML) models. The pipeline has been developed based on Snakemake workflow management system.
 
 

## Pipeline global behavior
The full pipeline should be executed in the following order :
 1. **run Sub-pipeline 1 : Process and training** --> ``` snakemake MaLCons/Process_and_Training/workflow/Snakefile ```
    - Stratified Sampling: execution of a data stratified sampling step for machine learning training
    - Training: launch ML models training
 2. **run script 1** --> ```python check_training_process.py```
    - check training status
 4. **run script 2** --> ```python create_config_Results_analysis.py```
    - automatic set up of the second sub-pipeline config.yaml file
 6. **run Sub-pipeline 2 : Results analysis** --> ``` snakemake MaLCons/Results_analysis/workflow/Snakefile ```
    - Standardize results: needed step to standardize results of various ML model training runs
    - Redundancy calculation: compute features identified as bearing redundant information during training (for each sampling run)
    - Models features relationships: retrieve associations between various ML models and the features name composing their signatures
    - Correlated features relationship: 
    - Correlated features: 
    - Statistics: Compute exploration of the global ML training through the various models, methods and runs.\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**Give a proposed consensus feature signature as final result**
 
 
## 1st sub-pipeline set up --> Process and training
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
- The training script build SLURM sbatch files to be able to run training jobs. It has been designed to run on Digital Research Alliance of Canada (the Alliance). You can adapt the sbatch file creation by modifying the `MaLCons/Process_and_training/workflow/scripts/training.py`.

## Output file tree structure
Here is an overview of the output file tree structure : 
```bash 
[ProjectName]/  
├── Results_[ProjectName]/  
│   └── STATS/  
│       ├── intermethod_intersampling/  
│       ├── intermethod_intrasampling/  
│       ├── intramethod_intrasampling/  
│       ├── unimethod_intersampling/  
│       ├── boxplot_MCC_and_STDMCC_params_sensi_onlyMCCopti_FB_DASH_LINE.png  
│       └── MCC07_STD01FILTERED_geomcountplot_proportionate_nbfeat.png  
├── Stratified_sampling_[ProjectName]/  
└── Training_[ProjectName]/  
    ├── S1x/  
    ├── .../  
    └── Snx/ (where n is the number of specified sampling runs)  
        ├── BAYES_A1DE/   
        ├── .../ (all types of trained classifiers)  
        ├── LAZY_kNN/  
        ├── INFOGAIN/  
        ├── NEO4J/  
        ├── PEARSON/  
        └── SPEARMAN/  
```
The ```Results_[ProjectName]``` directory contains all wrappers results and aggregation output as well as intermediate results : 
- STATS directory gathers all final and intermediate results where -->
   - intermethod_intersampling : highest level of aggregation as described in our paper. In the table, variables having an intersect of (1/2)n  sampling runs are composing the potential stable signature.
   - intermethod_intrasampling : Intermediate results, exploration of signature composition comparison between types of classifiers considering one sampling run at a time
   - intramethod_intrasampling : Intermediate results, exploration of signature composition comparison between models of a same type of classifier and considering one sampling run at a time
   - unimethod_intersampling : Intermediate results, exploration of signature composition comparison between one type of classifier at the time and through all sampling runs
   - boxplot is a general visualization of performance of all wrappers trained
   - geomcoutplot is an overview of all wrappers non redundant signature length
- Stratified_sampling_[ProjectName] directory gathers all pairs of train/test subsets for all sampling runs
- Training_[ProjectName] contains all sampling runs raw data
   - each Snx directory presents :
      - one directory/type of classifier with results of all associated trained models
      - an INFOGAIN directory : information gain scores to retrieve redundant features in the associated sampling run
      - a NEO4J directory : all data formated to be included in a NEO4J vizualisation if needed
      - PEARSON and SPEARMAN directories : all pearson and spearman correlation scores to retrieve additionnal redundant features in the associated sampling run 

## Requirements
```
- Python 3.8.13 (also tested with v 3.8.10):
   - pandas
   - sklearn
   - csv
   - os
   - shutil
   - Path
   - stat
- snakemake 7.12.0 (also tested with v 5.31.1)
```
