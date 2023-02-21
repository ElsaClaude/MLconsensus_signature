import pandas as pd
from sklearn.model_selection import train_test_split
import os
import shutil
from pathlib import Path
import stat

def absoluteFilePaths(directory):
    fullpath = []
    for dirpath,_,filenames in os.walk(directory):
        for f in filenames:
            fullpath.append(os.path.abspath(os.path.join(dirpath, f)))
    return fullpath

def create_tree_file(samplings,scriptpath,classifiers,outputdir, name,cpus,maxfeat,memory):
    for samp in range(1,samplings+1,1):
        samp = 'S'+str(samp)+'x'
        for classif in classifiers:
            classifname = '_'.join(classif.split('/')[-1].split('.')[0].split('_')[1:3])
            runpath = outputdir+'/Training/'+samp+'/'+classifname+'/'

            Path(runpath).mkdir(parents=True, exist_ok=True)

            # copy classifiers.conf in run path
            shutil.copy(classif,os.path.join(runpath, "classifiers.conf"))

            # create config.conf in run path
            config = open(runpath+'config.conf', "w")
            config.write('project='+name+'_'+samp+'_'+classifname+'\n')
            config.write('doClassification=true\nclassificationClassName=Label\nsampling=true\nseparator=,\nloocv=false\ncoptimizers=MCC\nsearchmodes=FB\ndebug=true\nwd=./\n')
            config.write('trainFile='+outputdir+'/Stratified_sampling/'+name+'_'+samp+'_'+'train.csv\n')
            config.write('validationFile='+outputdir+'/Stratified_sampling/'+name+'_'+samp+'_'+'test.csv\n')
            config.write('cpus='+str(cpus)+'\n')
            config.write('maxNumberOfSelectedFeatures='+maxfeat+'\n')
            config.close()

            # create training job runner AND run it
            job = open(runpath+'job_runner_'+name+'_'+samp+'_'+classifname+'.sh', "w")
            job.write('#!/bin/bash\n')
            job.write('#SBATCH --job-name='+samp+'_'+classifname+'_'+name+'\n')
            job.write('#SBATCH --cpus-per-task='+str(cpus)+'\n')
            job.write('#SBATCH --mem='+memory+'\n')
            job.write('#SBATCH --output='+outputdir+'/Training/'+samp+'/'+classifname+'training.out'+'\n')
            job.write('#SBATCH --error='+outputdir+'/Training/'+samp+'/'+classifname+'training.err'+'\n')

            ### optional
            job.write('#SBATCH --account=rrg-adroit-ab\n')

            ### Running time depending on classifier to train (POSSIBLE TO MODIFY)
            if classifname == 'BAYES_NB':
                job.write('#SBATCH --time=0-05:00:00\n')
            if classifname in ['BAYES_A1DE','TREES_C45']:
                job.write('#SBATCH --time=2-00:00:00\n')
            if classifname in ['TREES_CART','LAZY_kNN']:
                job.write('#SBATCH --time=5-00:00:00\n')
            if classifname in ['BAYES_BN','TREES_RF']:
                job.write('#SBATCH --time=7-00:00:00\n')
            if classifname == 'FUNCTIONS_SVM':
                job.write('#SBATCH --time=14-00:00:00\n')

            ### training command
            job.write('cd '+outputdir+'/Training/'+samp+'/'+classifname+';')
            job.write('module load java/1.8.0_192;')
            job.write('time java -Xmx'+memory+' -XX:+HeapDumpOnOutOfMemoryError -jar '+scriptpath+'/../../softs/BioDiscML/biodiscml_02_12_2022.jar -train -config config.conf')

            job.close()

            ### allows execution of the job runner
            st = os.stat(runpath+'job_runner_'+name+'_'+samp+'_'+classifname+'.sh')
            os.chmod(runpath+'job_runner_'+name+'_'+samp+'_'+classifname+'.sh', st.st_mode | stat.S_IEXEC)
    return

def main():
    outputdir = snakemake.params[0]
    name = snakemake.params[1]
    samplings = snakemake.params[2]
    cpus = snakemake.params[3]
    memory = snakemake.params[4]

    maxfeat = snakemake.input[0]
    with open(maxfeat,'r') as file:
        maxfeat = file.read().split('maxfeat=')[1]

    scriptpath = os.path.dirname(__file__)
    classifierspath = scriptpath+'/../../softs/classifiers_list/'
    classifiers = absoluteFilePaths(classifierspath)

    output = snakemake.output[0]

    create_tree_file(samplings,scriptpath,classifiers,outputdir, name,cpus,maxfeat, memory)

    f = open(output, "w")
    f.write("Training launched !")
    f.close()

if __name__ == "__main__":
    main()