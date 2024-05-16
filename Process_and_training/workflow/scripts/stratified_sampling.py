import pandas as pd
from sklearn.model_selection import train_test_split
import csv

# def get_delimiter(file_path, bytes = 4096):
#     sniffer = csv.Sniffer()
#     data = open(file_path, "r").read(bytes)
#     delimiter = sniffer.sniff(data).delimiter
#     print(delimiter)
#     return delimiter

def stratified_sampling(metadatapath,datapath,samplings,outputdir, name):
    metadatadf = pd.read_csv(metadatapath,sep=',',header=0)
    datadf = pd.read_csv(datapath,sep=',',header=0)

    fulldata = metadatadf.merge(datadf, on='ID')

    maxfeat=len(fulldata.columns)

    clusterslist = fulldata.pop('Cluster')
    truelabel = fulldata.pop('Label')

    for i in range(1,samplings+1,1):
        clusters_train, clusters_test, labels_train, labels_test = train_test_split( fulldata, truelabel, test_size=1/3, stratify=clusterslist)
        clusters_train = clusters_train.join(labels_train)
        clusters_test = clusters_test.join(labels_test)

        clusters_train.to_csv(outputdir+'/Stratified_sampling_'+name+'/'+name+'_S'+str(i)+'x_train.csv',sep=',',header=True,index=False)
        clusters_test.to_csv(outputdir+'/Stratified_sampling_'+name+'/'+name+'_S'+str(i)+'x_test.csv',sep=',',header=True, index=False)

    return maxfeat

def main():
    metadatapath = snakemake.input[0]
    datapath = snakemake.input[1]
    outputdir = snakemake.params[0]
    name = snakemake.params[1]
    samplings = snakemake.params[2]

    output = snakemake.output[0]

    maxfeat = stratified_sampling(metadatapath, datapath, samplings,outputdir, name)

    f = open(output, "w")
    f.write("Stratified sampling done !\n")
    f.write("maxfeat="+str(maxfeat))
    f.close()

if __name__ == "__main__":
    main()
