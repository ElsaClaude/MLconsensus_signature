from profile import run
from threading import main_thread
import pandas as pd
from pandas.io.parsers import read_csv
pd.options.mode.chained_assignment = None  # default='warn'
import numpy as np
from scipy.stats import iqr
from scipy import mean
from scipy import std
import os

def BDML_runname_models(BDML_results,runname):
    BDML_results['ID'] = runname +'_'+ BDML_results['ID']
    return(BDML_results)

def update_options_field(row):
    if row.str.contains('\"').any() :
        row = row.str.replace("\\","")
        row = row.str.replace("\"","\"\"")
        row = "\"" + row + "\""
    return(row)

def import_BDML_results(BDML_results_path,runname):
    listresults = []
    for file in BDML_results_path :
        filedf = pd.read_csv(file, index_col=None, header=0,sep=',')
        listresults.append(filedf)

    BDML_results_full = pd.concat(listresults, axis=0, ignore_index=True)
    ## Standardize options field
    BDML_results_full[['Options']] = BDML_results_full[['Options']].apply(lambda x: update_options_field(x) if(np.all(pd.notnull(x[0]))) else x, axis = 1)

    ## format test BDML results TO REMOVE ##
    BDML_results_full = BDML_runname_models(BDML_results_full,runname)
    BDML_results_full['Type_Of_Classifier'] = BDML_results_full['classifier'].str.rsplit('.').str[0]
    return BDML_results_full

def create_features_DF(BDML_info_gain,results):
    alllabels = pd.unique(BDML_info_gain['class'].tolist())

    BDML_info_gain = list(BDML_info_gain)
    BDML_info_gain = pd.DataFrame(BDML_info_gain, columns=['Unique_Feature_Name'])
    BDML_info_gain['count'] = BDML_info_gain.index + 1
    BDML_info_gain = BDML_info_gain.iloc[1:]
    BDML_info_gain = BDML_info_gain.iloc[:-1]

    FeaturesDF_tmp = list(results['AttributeList'])
    FeaturesDF = []
    for i in FeaturesDF_tmp :
        i = i.split(',')
        for y in i :
            FeaturesDF.append(y)

    FeaturesDF = list(dict.fromkeys(FeaturesDF))
    FeaturesDF = list(map(int, FeaturesDF))


    FeaturesDF = BDML_info_gain[BDML_info_gain['count'].isin(FeaturesDF)]

    FeaturesDF = FeaturesDF.drop(columns=['count'])

    results['Class_Label'] = '_VS_'.join(alllabels)

    return FeaturesDF,BDML_info_gain,results

def create_relationships(BDML_info_gain,results,fullinfogain):
    Relationships = results[['ID','AttributeList']]
    Relationships['AttributeList'] = Relationships['AttributeList'].str.split(',')
    Relationships = Relationships.explode('AttributeList')
    Relationships['AttributeList'] = pd.to_numeric(Relationships['AttributeList'])

    Relationships = Relationships.merge(BDML_info_gain,left_on='AttributeList',right_on='count')
    Relationships = Relationships.drop(columns=['count','AttributeList'])
    Relationships = Relationships.rename(columns={'ID': 'Model_ID_FROM', 'Unique_Feature_Name': 'Unique_Feature_Name_TO'})

    Relationships = Relationships.merge(fullinfogain, left_on = 'Unique_Feature_Name_TO', right_on = 'feature')
    Relationships = Relationships.drop(columns = ['origin_index','feature'])

    return Relationships

def classifiers_nodes(results,classifnames):
    results = pd.merge(results, classifnames, left_on='classifier', right_on='Method_WEKA_name').drop('classifier', axis='columns')

    classifnames['Mean_MCC'] = classifnames['Method_common_name'].map(results.groupby('Method_common_name')['AVG_MCC'].agg(mean))
    classifnames['STD_MCC'] = classifnames['Method_common_name'].map(results.groupby('Method_common_name')['AVG_MCC'].agg(std))
    classifnames['IQR_MCC'] = classifnames['Method_common_name'].map(results.groupby('Method_common_name')['AVG_MCC'].agg(iqr))

    classifnames['Mean_Signature_Length'] = classifnames['Method_common_name'].map(results.groupby('Method_common_name')['nbrOfFeatures'].agg(mean))
    classifnames['STD_Signature_Length'] = classifnames['Method_common_name'].map(results.groupby('Method_common_name')['nbrOfFeatures'].agg(std))
    classifnames['IQR_Signature_Length'] = classifnames['Method_common_name'].map(results.groupby('Method_common_name')['nbrOfFeatures'].agg(iqr))

    return classifnames

def output_csv(inputpath,run,dftocsv,run_name,name):
    dftocsv.to_csv(inputpath+run+'/NEO4J/'+run_name+'_'+run+'_'+name+'.csv',sep=',',header=True,index=False)
    return

def main():
    prefix = snakemake.params[0]
    inputpath = snakemake.params[1]

    fullinfogainfiles= snakemake.input[-len(snakemake.params[2]):]
    nbexp = len(snakemake.params[2]) * len(snakemake.params[3])
    BDML_results_path = snakemake.input[:nbexp]
    BDML_info_gain_path = snakemake.input[nbexp:-len(snakemake.params[2])]

    classifnames = snakemake.params[4]
    runname = snakemake.params[5]

    thresholdavgmcc = snakemake.params[6]
    thresholdavgmcc_name = str(thresholdavgmcc).replace('.','')

    thresholdstdmcc = snakemake.params[7]
    thresholdstdmcc_name = str(thresholdstdmcc).replace('.','')

    classifdict = {'Method_WEKA_name':[],'Method_common_name':[]}
    for classif in classifnames.items():
        weka = classif[1][0]['weka']
        classifdict['Method_WEKA_name'].append(weka)

        common = classif[1][1]['change']
        classifdict['Method_common_name'].append(common)
    
    classifnames = pd.DataFrame.from_dict(classifdict)

    for run in prefix :
        if not os.path.exists(inputpath+run+'/NEO4J') :
            print('NEO4J folder does not exist : '+run)
            os.makedirs(inputpath+run+'/NEO4J')
            print('NEO4J folder created')

        runbdmlresults = list(filter(lambda k: run in k.split('/'), BDML_results_path))

        runbdmlig = ''.join(list(filter(lambda k: run in k.split('/'), BDML_info_gain_path)))
        runfullig = ''.join(list(filter(lambda k: run in k.split('/'), fullinfogainfiles)))

        BDML_results = import_BDML_results(runbdmlresults,runname)
        BDML_thresholds = BDML_results.loc[(BDML_results['AVG_MCC'] >= thresholdavgmcc) & (BDML_results['STD_MCC'] <= thresholdstdmcc)]
        
        BDML_info_gain = pd.read_csv(runbdmlig, sep=',', header=0)
        fullinfogain = pd.read_csv(runfullig, sep=',', header=0)

        FeaturesDF_notfilt,BDML_info_gain_notfilt, BDML_results  = create_features_DF(BDML_info_gain, BDML_results)
        FeaturesDF_filt,BDML_info_gain_filt, BDML_thresholds = create_features_DF(BDML_info_gain, BDML_thresholds)

        Relationships_notfilt = create_relationships(BDML_info_gain_notfilt, BDML_results,fullinfogain)
        Relationships_filt = create_relationships(BDML_info_gain_filt, BDML_thresholds,fullinfogain)

        methods_notfilt = classifiers_nodes(BDML_results,classifnames)
        methods_filt = classifiers_nodes(BDML_thresholds,classifnames)

        output_csv(inputpath,run,methods_notfilt,runname,'METHODS')
        output_csv(inputpath,run,methods_filt,runname,'METHODS_FILTERED_MCC'+thresholdavgmcc_name+'_STD'+thresholdstdmcc_name)

        output_csv(inputpath,run,BDML_results,runname,'MODELS')
        output_csv(inputpath,run,BDML_thresholds,runname,'MODELS_FILTERED_MCC'+thresholdavgmcc_name+'_STD'+thresholdstdmcc_name)

        output_csv(inputpath,run,FeaturesDF_notfilt,runname,'FEATURES')
        output_csv(inputpath,run,FeaturesDF_filt,runname,'FEATURES_FILTERED_MCC'+thresholdavgmcc_name+'_STD'+thresholdstdmcc_name)

        output_csv(inputpath,run,Relationships_notfilt,runname,'RELATIONSHIPS')
        output_csv(inputpath,run,Relationships_filt,runname,'RELATIONSHIPS_FILTERED_MCC'+thresholdavgmcc_name+'_STD'+thresholdstdmcc_name)

if __name__ == "__main__":
    main()

