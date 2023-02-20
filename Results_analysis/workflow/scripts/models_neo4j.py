from profile import run
from threading import main_thread
import pandas as pd
from pandas.io.parsers import read_csv
pd.options.mode.chained_assignment = None  # default='warn'
import numpy as np
import mygene
import argparse
from scipy.stats import iqr
from scipy import mean
from scipy import std
import os


## TO-DO list : 
#           - automate line 115 (FeaturesDF = convert_ensemblID_GENES_to_symbol(FeaturesDF))

def convert_ensemblID_GENES_to_symbol(FeaturesDF,refdir):
    FeaturesDF.loc[FeaturesDF['TypeOfFeature'] == 'Gene','EnsemblID'] = FeaturesDF['Name']

    FeaturesDF['EnsemblID'] = FeaturesDF['EnsemblID'].str.replace('(\.).*', '', regex=True)
    listtoconvert = FeaturesDF.loc[FeaturesDF.TypeOfFeature == 'Gene','EnsemblID'].values.tolist()
    mg = mygene.MyGeneInfo()
    convertedlist = mg.querymany(listtoconvert, scopes="ensembl.gene", fields="symbol", species="human",returnall=True)
    del convertedlist['missing']
    del convertedlist['dup']

    convertedlist = pd.DataFrame(convertedlist['out'])
    convertedlist = convertedlist.drop(columns=['_id','_score','notfound'],errors='ignore')
    convertedlist.loc[convertedlist['symbol'].isnull(),'symbol'] = ''
    # convertedlist['symbol'] = convertedlist['symbol'].fillna(convertedlist['query'], inplace=True)

    FeaturesDF = FeaturesDF.merge(convertedlist,left_on='EnsemblID', right_on='query')
    FeaturesDF = FeaturesDF.drop(columns=['query'])
    FeaturesDF = FeaturesDF.rename(columns={'symbol': 'GeneSymbolName'})

    genebiotype = read_csv(refdir+'human_ensembl_id_genes_AND_gene_biotype.csv',delimiter = ',',header=0)

    FeaturesDF = FeaturesDF.merge(genebiotype, left_on='EnsemblID', right_on='ensembl_gene_id',how='left')
    FeaturesDF = FeaturesDF.drop(columns=['ensembl_gene_id'])

    return FeaturesDF

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
    
    ## Keep only interesting models (THREHSOLDS TO MODIFY)
    # BDML_results = BDML_results.loc[(BDML_results['AVG_MCC'] >= 0.90) & (BDML_results['STD_MCC'] <= 0.05)]

    ## format test BDML results TO REMOVE ##
    BDML_results_full = BDML_runname_models(BDML_results_full,runname)
    BDML_results_full['Type_Of_Classifier'] = BDML_results_full['classifier'].str.rsplit('.').str[0]
    return BDML_results_full

def create_features_DF(BDML_info_gain,results,refdir):
    print(BDML_info_gain.columns[-2:])
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

    ###### WARNING : modify to suit better final nomenclature [?]
    FeaturesDF['Origin'] = FeaturesDF['Unique_Feature_Name'].str.split('_').str[0]
    FeaturesDF['TypeOfFeature'] = FeaturesDF['Unique_Feature_Name'].str.split('_').str[1]
    FeaturesDF['Name'] = FeaturesDF['Unique_Feature_Name'].str.split('_').str[3]
    FeaturesDF = FeaturesDF.drop(columns=['count'])
    ### automate this line to be able to apply a function of search no matter if feature is ensemblid, clinical data or else

    if 'Gene' in FeaturesDF[['TypeOfFeature']].values :
        FeaturesDF = convert_ensemblID_GENES_to_symbol(FeaturesDF,refdir)

    FeaturesDF = FeaturesDF.drop(columns=['Name'])

    usedtypeoffeatures = pd.unique(FeaturesDF['Origin'].tolist())
    usedtypeoffeatures = ' / '.join(usedtypeoffeatures)

    results['Data_Origin'] = usedtypeoffeatures
    results['Class_Label'] = '_VS_'.join(alllabels)

    return FeaturesDF,BDML_info_gain,results

def create_relationships(BDML_info_gain,results,fullinfogain):
    Relationships = results[['ID','AttributeList']]
    Relationships['AttributeList'] = Relationships['AttributeList'].str.split(',')
    Relationships = Relationships.explode('AttributeList')
    Relationships['AttributeList'] = pd.to_numeric(Relationships['AttributeList'])

    Relationships = Relationships.merge(BDML_info_gain,left_on='AttributeList',right_on='count')
    # Relationships = Relationships.merge(old_new,left_on='Unique_Feature_Name',right_on='Old_Unique')
    Relationships = Relationships.drop(columns=['count','AttributeList'])
    Relationships = Relationships.rename(columns={'ID': 'Model_ID_FROM', 'Unique_Feature_Name': 'Unique_Feature_Name_TO'})

    Relationships = Relationships.merge(fullinfogain, left_on = 'Unique_Feature_Name_TO', right_on = 'feature')
    Relationships = Relationships.drop(columns = ['origin_index','feature'])

    return Relationships

def methods_nodes(results,nboptions):
    methods = pd.read_csv(nboptions)
    results = pd.merge(results, methods[['Method_WEKA_name','Method_common_name']], left_on='classifier', right_on='Method_WEKA_name').drop('classifier', axis='columns')

    methods['Mean_MCC'] = methods['Method_common_name'].map(results.groupby('Method_common_name')['AVG_MCC'].agg(mean))
    methods['STD_MCC'] = methods['Method_common_name'].map(results.groupby('Method_common_name')['AVG_MCC'].agg(std))
    methods['IQR_MCC'] = methods['Method_common_name'].map(results.groupby('Method_common_name')['AVG_MCC'].agg(iqr))

    methods['Mean_Signature_Length'] = methods['Method_common_name'].map(results.groupby('Method_common_name')['nbrOfFeatures'].agg(mean))
    methods['STD_Signature_Length'] = methods['Method_common_name'].map(results.groupby('Method_common_name')['nbrOfFeatures'].agg(std))
    methods['IQR_Signature_Length'] = methods['Method_common_name'].map(results.groupby('Method_common_name')['nbrOfFeatures'].agg(iqr))

    return methods

def output_csv(inputpath,run,dftocsv,run_name,name):
    dftocsv.to_csv(inputpath+run+'/NEO4J/'+run_name+'_'+run+'_'+name+'.csv',sep=',',header=True,index=False)
    return

def main():
    # BDML_results_path,BDML_info_gain_path, nboptions, inputpath,refdir, runname , fullinfogain= get_arg()
    prefix = snakemake.params[0]
    inputpath = snakemake.params[1]

    # print(prefix)
    # print(inputpath)

    fullinfogainfiles= snakemake.input[-len(snakemake.params[2]):]
    nbexp = len(snakemake.params[2]) * len(snakemake.params[3])
    BDML_results_path = snakemake.input[:nbexp]
    BDML_info_gain_path = snakemake.input[nbexp:-len(snakemake.params[2])]
    # print(BDML_info_gain_path)



    nboptions = snakemake.params[4]
    refdir = snakemake.params[5]
    runname = snakemake.params[6]

    thresholdavgmcc = snakemake.params[7]
    thresholdavgmcc_name = str(thresholdavgmcc).replace('.','')

    thresholdstdmcc = snakemake.params[8]
    thresholdstdmcc_name = str(thresholdstdmcc).replace('.','')


    for run in prefix :
        print(run+' IS RUN')
        if not os.path.exists(inputpath+run+'/NEO4J') :
            print('NEO4J folder does not exist : '+run)
            os.makedirs(inputpath+run+'/NEO4J')
            print('NEO4J folder created')

        runbdmlresults = list(filter(lambda k: run in k.split('/'), BDML_results_path))
        # print(runbdmlresults)
        runbdmlig = ''.join(list(filter(lambda k: run in k.split('/'), BDML_info_gain_path)))
        print(runbdmlig)
        runfullig = ''.join(list(filter(lambda k: run in k.split('/'), fullinfogainfiles)))

        BDML_results = import_BDML_results(runbdmlresults,runname)
        BDML_thresholds = BDML_results.loc[(BDML_results['AVG_MCC'] >= thresholdavgmcc) & (BDML_results['STD_MCC'] <= thresholdstdmcc)]
        
        BDML_info_gain = pd.read_csv(runbdmlig, sep=',', header=0)
        fullinfogain = pd.read_csv(runfullig, sep=',', header=0)

        FeaturesDF_notfilt,BDML_info_gain_notfilt, BDML_results  = create_features_DF(BDML_info_gain, BDML_results,refdir)
        FeaturesDF_filt,BDML_info_gain_filt, BDML_thresholds = create_features_DF(BDML_info_gain, BDML_thresholds,refdir)

        Relationships_notfilt = create_relationships(BDML_info_gain_notfilt, BDML_results,fullinfogain)
        Relationships_filt = create_relationships(BDML_info_gain_filt, BDML_thresholds,fullinfogain)

        methods_notfilt = methods_nodes(BDML_results,nboptions)
        methods_filt = methods_nodes(BDML_thresholds,nboptions)

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

