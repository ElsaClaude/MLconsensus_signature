from pyexpat import features
import pandas as pd
from pandas.io.parsers import read_csv
pd.options.mode.chained_assignment = None  # default='warn'
import numpy as np
from functools import reduce
import argparse

def create_relationships(NEO4JmfRelationShips,fullIG,pears ,spear,inputpath,runname,run,suffix):
    NEO4JmfRelationShips = pd.read_csv(NEO4JmfRelationShips,header=0)
    fullIG = pd.read_csv(fullIG,header=0)
    pearson = pd.read_csv(pears,header=0)
    spearman = pd.read_csv(spear,header=0)

    targetfeatures = NEO4JmfRelationShips['Unique_Feature_Name_TO'].values
    targetfeatures = pd.unique(targetfeatures)

    targetIGscores = fullIG.loc[fullIG['feature'].isin(targetfeatures), ['feature','info_gain_score']]

    correlatedIG = targetIGscores.merge(fullIG, on = 'info_gain_score')
    correlatedIG = correlatedIG.drop(columns=['origin_index'])
    correlatedIG = correlatedIG.rename(columns={'feature_x': 'target_feature', 'feature_y': 'correlated_feature'})
    correlatedIG = correlatedIG[correlatedIG['target_feature'] != correlatedIG['correlated_feature']]

    correlatedIG['check_string'] = correlatedIG.apply(lambda row: ''.join(sorted([row['target_feature'], row['correlated_feature']])), axis=1)
    correlatedIG = correlatedIG.drop_duplicates('check_string')
    correlatedIG = correlatedIG.drop(columns=['check_string'])

    correlatedPearson = pearson.loc[(pearson['index_feat_1'].isin(targetfeatures)) | (pearson['index_feat_2'].isin(targetfeatures))]
    correlatedPearson['index_feat_1'],correlatedPearson['index_feat_2']=np.where(correlatedPearson['index_feat_2'].isin(targetfeatures),(correlatedPearson['index_feat_2'],correlatedPearson['index_feat_1']),(correlatedPearson['index_feat_1'],correlatedPearson['index_feat_2']))
    correlatedPearson = correlatedPearson.rename(columns={'index_feat_1': 'target_feature', 'index_feat_2': 'correlated_feature','corr_coeff' : 'Pearson_corr_coeff'})

    correlatedSpearman = spearman.loc[(spearman['index_feat_1'].isin(targetfeatures)) | (spearman['index_feat_2'].isin(targetfeatures))]
    correlatedSpearman['index_feat_1'],correlatedSpearman['index_feat_2']=np.where(correlatedSpearman['index_feat_2'].isin(targetfeatures),(correlatedSpearman['index_feat_2'],correlatedSpearman['index_feat_1']),(correlatedSpearman['index_feat_1'],correlatedSpearman['index_feat_2']))
    correlatedSpearman = correlatedSpearman.rename(columns={'index_feat_1': 'target_feature', 'index_feat_2': 'correlated_feature','rank_corr_coeff' : 'Spearman_rank_corr_coeff'})

    fullcorrelated = reduce(lambda  left,right: pd.merge(left,right,on=['target_feature','correlated_feature'],
                                            how='outer'), [correlatedIG,correlatedPearson,correlatedSpearman])
    fullcorrelated = fullcorrelated.fillna('')
    fullcorrelated['run_name'] = runname

    fullcorrelated.to_csv(inputpath+run+'/NEO4J/'+runname+'_'+run+'_RELCORRFEAT'+suffix+'.csv',sep=',',index=False)

def main():
    prefix = snakemake.params[0]
    nbruns = len(prefix)

    inputpath = snakemake.params[1]

    runname = snakemake.params[2]

    thresholdavgmcc = str(snakemake.params[3])
    thresholdavgmcc = thresholdavgmcc.replace('.','')

    thresholdstdmcc = str(snakemake.params[4])
    thresholdstdmcc = thresholdstdmcc.replace('.','')

    runs_fullig = snakemake.input[:nbruns]

    runs_pearson = snakemake.input[nbruns:nbruns*2]

    runs_spearman = snakemake.input[nbruns*2:nbruns*3]

    runs_relationshipsmf = snakemake.input[nbruns*3:]
    
    for run in prefix :
        NEO4JmfRelationShips = list(filter(lambda k: run in k.split('/'), runs_relationshipsmf))
        print(NEO4JmfRelationShips)
        NEO4JmfRelationShips_notfilt = ''.join(list(filter(lambda k: "FILTERED" not in k, NEO4JmfRelationShips)))
        print(NEO4JmfRelationShips_notfilt)

        NEO4JmfRelationShips_filt = ''.join(list(filter(lambda k: "FILTERED" in k, NEO4JmfRelationShips)))
        print(NEO4JmfRelationShips_filt)


        fullIG = ''.join(list(filter(lambda k: run in k.split('/'), runs_fullig)))
        print(fullIG)

        pears = ''.join(list(filter(lambda k: run in k.split('/'), runs_pearson)))

        spear = ''.join(list(filter(lambda k: run in k.split('/'), runs_spearman)))

        create_relationships(NEO4JmfRelationShips_notfilt,fullIG,pears ,spear,inputpath,runname, run, "")
        create_relationships(NEO4JmfRelationShips_filt,fullIG,pears ,spear,inputpath,runname,run,"_FILTERED_MCC"+thresholdavgmcc+'_STD'+thresholdstdmcc)
        
if __name__ == "__main__":
    main()
