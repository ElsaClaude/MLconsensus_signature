import os

def absoluteFilePaths(directory):
    fullpath = []
    for dirpath,_,filenames in os.walk(directory):
        for f in filenames:
            fullpath.append(os.path.abspath(os.path.join(dirpath, f)))
    return fullpath

def create_config(configtraining, classifiers,scriptpath):
    with open(configtraining,'r') as training:
        training = training.read()
        outputdir = training.split('output_directory: ')[1].split('\n')[0]
        samplings = int(training.split('nb_samplings: ')[1].split('\n')[0])
        name = training.split('name: ')[1].split('\n')[0]

    with open(scriptpath+'/Results_analysis/config/config.yaml','w') as results:
        results.write('prefix:\n')
        for samp in range(1,samplings+1,1):
            results.write('  - S'+str(samp)+'x\n')

        results.write('run_name: '+name+'\n')

        ## list of used ML algorithms
        results.write('classifiers:\n')
        for classif in classifiers :
            classif = '_'.join(classif.split('/')[-1].split('.')[0].split('_')[1:3])
            results.write('  - '+classif+'\n')

        ## write classifiers knownm name and equivalent in weka doc
        results.write('classifiers_weka_name:\n')
        results.write('  BAYES_A1DE:\n    - weka: bayes.AveragedNDependenceEstimators.A1DE\n    - change: A1DE\n')
        results.write('  BAYES_BN:\n    - weka: bayes.BayesNet\n    - change: Bayes Network\n')
        results.write('  BAYES_NB:\n    - weka: bayes.NaiveBayes\n    - change: Naive Bayes\n')
        results.write('  TREES_C45:\n    - weka: trees.J48\n    - change: C4.5\n')
        results.write('  TREES_RF:\n    - weka: trees.RandomForest\n    - change: Random Forest\n')
        results.write('  TREES_CART:\n    - weka: trees.SimpleCart\n    - change: Simple Cart\n')
        results.write('  LAZY_kNN:\n    - weka: lazy.IBk\n    - change: kNN\n')
        results.write('  FUNCTIONS_SVM:\n    - weka: functions.SMO\n    - change: SVM\n')

        ## AVG MCC threshold (by default : 0.7)
        results.write('avg_mcc_threshold: 0.7\n')

        ## STD MCC threshold (by default : 0.1)
        results.write('std_mcc: 0.1\n')

        ## path to training results files
        results.write('input_path: '+outputdir+'/Training_'+name+'/\n')
        results.write('output_path: '+outputdir+'/Results_'+name+'/')
    return

def main():
    scriptpath = os.path.dirname(os.path.abspath(__file__))
    configtraining = scriptpath+'/Process_and_training/config/config.yaml'

    classifiers = scriptpath+'/Process_and_training/softs/classifiers_list/'
    classifiers = absoluteFilePaths(classifiers)

    create_config(configtraining,classifiers,scriptpath)

if __name__ == "__main__":
    main()