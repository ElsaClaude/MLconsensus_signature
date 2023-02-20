import pandas as pd
from sklearn.model_selection import train_test_split
import csv

# def create_tree_file(metadatapath,datapath,outputdir, name):
#     metadatadf = pd.read_csv(metadatapath,sep=get_delimiter(metadatapath),header=0)
#     datadf = pd.read_csv(datapath,sep=get_delimiter(datapath),header=0)

#     fulldata = metadatadf.merge(datadf, on='ID')

#     clusterslist = fulldata.pop('Cluster')
#     truelabel = fulldata.pop('Label')

#     for i in range(1,11,1):
#         print(i)
#         clusters_train, clusters_test, labels_train, labels_test = train_test_split( fulldata, truelabel, test_size=1/3, stratify=clusterslist)
#         clusters_train = clusters_train.join(labels_train)
#         clusters_test = clusters_test.join(labels_test)

#         clusters_train.to_csv(outputdir+'/Stratified_sampling/'+name+'_StratSamp_'+str(i)+'x_train.csv',sep=',',header=True,index=False)
#         clusters_test.to_csv(outputdir+'/Stratified_sampling/'+name+'_StratSamp_'+str(i)+'x_test.csv',sep=',',header=True, index=False)

#     return

def main():
    metadatapath = snakemake.input[0]
    datapath = snakemake.input[1]
    outputdir = snakemake.params[0]
    name = snakemake.params[1]

    output = snakemake.output[0]

    stratified_sampling(metadatapath, datapath, outputdir, name)

    f = open(output, "a")
    f.write("Stratified sampling done !")
    f.close()

if __name__ == "__main__":
    main()