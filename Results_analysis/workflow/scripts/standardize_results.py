import pandas as pd
from pandas.io.parsers import read_csv
pd.options.mode.chained_assignment = None  # default='warn'
import numpy as np

def update_file(input,outputpath):
    for inputfile in input :
        filename = inputfile.split('.')[:-1]
        filename = '.'.join(filename)
        inputdf = read_csv(inputfile, delimiter='\t')

        inputdf = inputdf.drop(inputdf[inputdf['ID'] == 'Resumed here'].index)
        print(inputdf.ID)
        inputdf = inputdf[~inputdf.ID.str.contains("error", na=True)]

        ## Standardize ID field to remove duplicates
        indexconcat = inputdf.groupby(inputdf['ID'])\
            .cumcount()\
            .astype(str)
        indexconcat = ['_' + s for s in indexconcat]
        indexconcat = np.array(indexconcat, dtype=str)
        indexconcat = np.where(indexconcat=='_0','',indexconcat)
        
        inputdf['ID'] = (inputdf['ID'].astype(str) + indexconcat.astype(str))

        inputdf.to_csv(filename+'.STANDARDIZED.csv',sep=',',header=True,index=False)
    return input

def main():
    outputpath = snakemake.params[2]
    input = snakemake.input

    update_file(input,outputpath)

if __name__ == "__main__":
    main()
