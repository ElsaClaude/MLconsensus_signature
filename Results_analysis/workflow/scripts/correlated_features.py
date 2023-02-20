import pandas as pd
from pandas.io.parsers import read_csv
pd.options.mode.chained_assignment = None  # default='warn'
import mygene

def convert_ensemblID_GENES_to_symbol(FeaturesDF,refdir):
    FeaturesDF.loc[FeaturesDF['TypeOfFeature'] == 'Gene','EnsemblID'] = FeaturesDF['Name']
    FeaturesDF['EnsemblID'] = FeaturesDF['EnsemblID'].str.replace('(\.).*', '', regex=True)

    listtoconvert = FeaturesDF.loc[FeaturesDF.TypeOfFeature == 'Gene','EnsemblID'].values.tolist()
    mg = mygene.MyGeneInfo()
    convertedlist = mg.querymany(listtoconvert, scopes="ensembl.gene", fields="symbol", species="human",returnall=True)

    del convertedlist['missing']
    del convertedlist['dup']

    convertedlist = pd.DataFrame.from_dict(convertedlist['out'])
    convertedlist = convertedlist.drop(columns=['_id','_score','notfound'],errors='ignore')
    convertedlist.loc[convertedlist['symbol'].isnull(),'symbol'] = ''
    # convertedlist['symbol'] = convertedlist['symbol'].fillna(convertedlist['query'], inplace=True)

    FeaturesDF = FeaturesDF.merge(convertedlist,left_on='EnsemblID', right_on='query')
    FeaturesDF = FeaturesDF.drop(columns=['query'])
    FeaturesDF = FeaturesDF.rename(columns={'symbol': 'GeneSymbolName'})

    genebiotype = read_csv(refdir+'human_ensembl_id_genes_AND_gene_biotype.csv',delimiter = ',',header=0)

    FeaturesDF = FeaturesDF.merge(genebiotype, left_on='EnsemblID', right_on='ensembl_gene_id', how='left')
    FeaturesDF = FeaturesDF.drop(columns=['ensembl_gene_id'])
    return FeaturesDF

def create_corr_feat(NEO4Jcorr,NEO4Jfeat, refdir,inputpath , run,outputname,suffix):
    corr = pd.read_csv(NEO4Jcorr,header=0)
    shortfeat = pd.read_csv(NEO4Jfeat,header=0)

    toremove = pd.unique(shortfeat['Unique_Feature_Name'].tolist()).tolist()
    toadd = pd.unique(corr['target_feature'].tolist()).tolist() + pd.unique(corr['correlated_feature'].tolist()).tolist()
    toadd = list(set(toadd) - set(toremove))
    toadd = pd.DataFrame(toadd,columns=['Unique_Feature_Name'])

    toadd['tmp'] = toadd['Unique_Feature_Name'].str.split('__').str[0]
    toadd['Name'] = toadd['Unique_Feature_Name'].str.split('__').str[1]
    toadd[['Origin', 'TypeOfFeature']] = toadd['tmp'].str.split('_', 1, expand=True)
    toadd = toadd.drop(columns=['tmp'])

    if 'Gene' in toadd[['TypeOfFeature']].values :
        toadd = convert_ensemblID_GENES_to_symbol(toadd,refdir)

### Make also function for other feature type
    toadd = toadd.drop(columns=['Name'])

    toadd.to_csv(inputpath+run+'/NEO4J/'+outputname+'_'+run+'_CORRFEAT'+suffix+'.csv',sep=',',index=False)

def main():
    prefix = snakemake.params[0]
    nbrun = len(prefix)

    inputpath = snakemake.params[1]

    runname = snakemake.params[2]

    thresholdavgmcc = str(snakemake.params[3])
    thresholdavgmcc = thresholdavgmcc.replace('.','')

    thresholdstdmcc = str(snakemake.params[4])
    thresholdstdmcc = thresholdstdmcc.replace('.','')

    refdir = snakemake.params[5]

    NEO4Jcorr_files = list(filter(lambda k: "RELCORRFEAT" in k, snakemake.input))
    NEO4Jfeat_files = list(filter(lambda k: "FEATURES" in k, snakemake.input))

    for run in prefix :
        NEO4Jcorr = list(filter(lambda k: run in k.split('/'), NEO4Jcorr_files))
        NEO4Jcorr_notfilt = ''.join(list(filter(lambda k: "FILTERED" not in k.split('/'), NEO4Jcorr)))
        NEO4Jcorr_filt = ''.join(list(filter(lambda k: "FILTERED" in k.split('/'), NEO4Jcorr)))


        NEO4Jfeat = list(filter(lambda k: run in k.split('/'), NEO4Jfeat_files))
        NEO4Jfeat_notfilt = ''.join(list(filter(lambda k: "FILTERED" not in k.split('/'), NEO4Jfeat)))
        NEO4Jfeat_filt = ''.join(list(filter(lambda k: "FILTERED" in k.split('/'), NEO4Jfeat)))

        create_corr_feat(NEO4Jcorr_notfilt,NEO4Jfeat_notfilt, refdir,inputpath , run,runname,"")
        create_corr_feat(NEO4Jcorr_filt,NEO4Jfeat_filt, refdir,inputpath , run,runname,"_FILTERED_MCC"+thresholdavgmcc+'_STD'+thresholdstdmcc)

    # f = open(inputpath+"finished.txt", "a")
    # f.write("Now the file has more content!")
    # f.close()


if __name__ == "__main__":
    main()
