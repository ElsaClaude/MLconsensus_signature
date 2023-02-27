import os

def absoluteFilePaths(directory):
    fullpath = []
    for dirpath,_,filenames in os.walk(directory):
        for f in filenames:
            fullpath.append(os.path.abspath(os.path.join(dirpath, f)))
    return fullpath

def get_status(configtraining, classifiers,scriptpath):
    with open(configtraining,'r') as file:
        file = file.read()
        outputdir = file.split('output_directory: ')[1].split('\n')[0]
        samplings = int(file.split('nb_samplings: ')[1].split('\n')[0])
        name = file.split('name: ')[1].split('\n')[0]

    for samp in range(1,samplings+1,1):
        samp = 'S'+str(samp)+'x'
        for classif in classifiers:
            classifname = '_'.join(classif.split('/')[-1].split('.')[0].split('_')[1:3])
            stdout = outputdir+'/Training_'+name+'/'+samp+'/'+classifname+'/training.out'

            ## count trained models
            with open(stdout,'r') as stdout:
                stdout = stdout.read()
                trainedml = int(stdout.split('[model] (')[-1].split('/')[0])
            trainedml = trainedml + 1

            ## count models TO train (origin)
            classifconf = outputdir+'/Training_'+name+'/'+samp+'/'+classifname+'/classifiers.conf'
            nbtotrain = 0
            with open(classifconf,'r') as classifconf:
                for line in classifconf.readlines() :
                    if line.startswith('ccmd='):
                        nbtotrain += 1

            ## check if trained models < models to train
            ## if yes, print out the run
            if trainedml < nbtotrain :
                print('[model(s) not submitted to training] Run = '+name+', Sampling = '+samp+', Classifier = '+classifname+'\tTrained/ToTrain = '+str(trainedml)+'/'+str(nbtotrain))


    return

def main():
    scriptpath = os.path.dirname(os.path.abspath(__file__))
    configtraining = scriptpath+'/Process_and_training/config/config.yaml'

    classifiers = scriptpath+'/Process_and_training/softs/classifiers_list/'
    classifiers = absoluteFilePaths(classifiers)

    get_status(configtraining,classifiers,scriptpath)

if __name__ == "__main__":
    main()