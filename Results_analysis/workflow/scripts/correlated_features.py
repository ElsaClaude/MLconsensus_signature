import pandas as pd
from pandas.io.parsers import read_csv
pd.options.mode.chained_assignment = None  # default='warn'

def create_corr_feat(NEO4Jcorr,NEO4Jfeat,inputpath , run,outputname,suffix):
    corr = pd.read_csv(NEO4Jcorr,header=0)
    shortfeat = pd.read_csv(NEO4Jfeat,header=0)

    toremove = pd.unique(shortfeat['Unique_Feature_Name'].tolist()).tolist()
    toadd = pd.unique(corr['target_feature'].tolist()).tolist() + pd.unique(corr['correlated_feature'].tolist()).tolist()
    toadd = list(set(toadd) - set(toremove))
    toadd = pd.DataFrame(toadd,columns=['Unique_Feature_Name'])

    print(toadd)

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

    NEO4Jcorr_files = list(filter(lambda k: "RELCORRFEAT" in k, snakemake.input))
    NEO4Jfeat_files = list(filter(lambda k: "FEATURES" in k, snakemake.input))

    for run in prefix :
        NEO4Jcorr = list(filter(lambda k: run in k.split('/'), NEO4Jcorr_files))
        NEO4Jcorr_notfilt = ''.join(list(filter(lambda k: "FILTERED" not in k.split('_'), NEO4Jcorr)))
        NEO4Jcorr_filt = ''.join(list(filter(lambda k: "FILTERED" in k.split('_'), NEO4Jcorr)))

        NEO4Jfeat = list(filter(lambda k: run in k.split('/'), NEO4Jfeat_files))
        NEO4Jfeat_notfilt = ''.join(list(filter(lambda k: "FILTERED" not in k.split('_'), NEO4Jfeat)))
        NEO4Jfeat_filt = ''.join(list(filter(lambda k: "FILTERED" in k.split('_'), NEO4Jfeat)))

        create_corr_feat(NEO4Jcorr_notfilt,NEO4Jfeat_notfilt,inputpath , run,runname,"")
        create_corr_feat(NEO4Jcorr_filt,NEO4Jfeat_filt,inputpath , run,runname,"_FILTERED_MCC"+thresholdavgmcc+'_STD'+thresholdstdmcc)

if __name__ == "__main__":
    main()
