import pandas as pd
from pandas.io.parsers import read_csv
pd.options.mode.chained_assignment = None  # default='warn'
import numpy as np
import os
import math
import multiprocessing as mp
import operator as op
from functools import reduce

import re
import weka.core.jvm as jvm
from weka.attribute_selection import ASSearch, ASEvaluation, AttributeSelection
import weka.core.converters as converters
from weka.core.dataset import create_instances_from_lists,create_instances_from_matrices
from weka.filters import Filter


def ncr(n: int, r: int) -> int:
    """
    Calculates the number of different, unordered combination 
    of r objects from a set of n objects.
    nCr(n,r) = nPr(n,r) / r!

    Args:
        n (int): Number of objects in set
        r (int): Number of objects to combine
    Returns:
        int: Number of combinations
    """
    r = min(r, n-r)
    numer = reduce(op.mul, range(n, n-r, -1), 1)
    denom = reduce(op.mul, range(1, r+1), 1)
    return numer // denom

def processor_fun_pearson(*,outputfile : str, i_start: int, i_end: int, datam, datass) -> None:
    """
    Calculate correlation of rows in the range of i_start and i_end
    and the rows below these indices.
Args:
        i_start (int): Index of the start row
        i_end (int): Index of the end row
    Returns:
        None
    """    
    for i in range(i_start, i_end):
        temp = np.dot(datam[i:], datam[i].T)
        rs = temp / (datass[i:] * datass[i])
        rs = np.delete(rs, [0])
        # Create nparray that contains information about the indices
        res = np.zeros((rs.shape[0], 3))
        res[:, 0] = i
        res[:, 1] = list(range(i+1, i + rs.shape[0] + 1))
        res[:, 2] = rs

        # delete row where absolute value of pearson corr coeff is less than threshold (0.985)
        mask = (abs(res[:, 2]) > 0.985)
        res = res[mask]

        # append results to the final csv file
        with open(outputfile, "ab") as f:
            np.savetxt(f, res,
                   delimiter=',',
                   newline='\n', 
                   fmt=['%d','%d','%0.13f'])


def processor_fun_spearman(*,outputfile : str, i_start: int, i_end: int, denum, ranks) -> None:
    """
    Calculate correlation of rows in the range of i_start and i_end
    and the rows below these indices.
Args:
        i_start (int): Index of the start row
        i_end (int): Index of the end row
    Returns:
        None
    """    
    for i in range(i_start, i_end):
        d = np.subtract(ranks[i:], ranks[i].T)
        dsquare = d**2
        sumdsquare = np.sum(dsquare, axis=1)
        rs = 1-((6*sumdsquare)/denum)
        rs = np.delete(rs, [0])
        # Create nparray that contains information about the indices
        res = np.zeros((rs.shape[0], 3))
        res[:, 0] = i
        res[:, 1] = list(range(i+1, i + rs.shape[0] + 1))
        res[:, 2] = rs

        # delete row where absolute value of pearson corr coeff is less than threshold (0.985)
        mask = (abs(res[:, 2]) > 0.985)
        res = res[mask]

        # append results to the final csv file
        with open(outputfile, "ab") as f:
            np.savetxt(f, res,
                   delimiter=',',
                   newline='\n', 
                   fmt=['%d','%d','%0.13f'])

def listener(q,file):
    '''listens for messages on the q, writes to file. '''

    with open(file, 'ab') as f:
        while 1:
            m = q.get()
            if m == 'kill':
                break
            f.write(str(m) + '\n')
            f.flush()

def map_index_to_feature_name(fulldata,file):
    header = pd.read_csv(fulldata, index_col=0, nrows=0,sep='\t').columns.tolist()
    indextoname = pd.DataFrame(header, columns=['feature_name'])
    indextoname['origin_index'] = range(0, len(indextoname))

    toreplace = pd.read_csv(file,header=0)
    m = indextoname.set_index('origin_index')['feature_name'].to_dict()
    v = toreplace.filter(items=['index_feat_1','index_feat_2'])
    toreplace[v.columns] = v.replace(m)
    toreplace.to_csv(file,sep=',',header=True,index=False)

def get_pearson_corr_matrices(fulldata,inputpath,run):
    concatfilepearson = inputpath+run+'/PEARSON/pearson_corr_CONCAT.csv'

    file = open(concatfilepearson, "w")
    file.write('index_feat_1,index_feat_2,corr_coeff\n')
    file.close()

    items = pd.read_csv(fulldata, header = 0,sep='\t')
    items = items.drop(columns=['Instance','class'])
    items = items.T
    numpy_items = items.to_numpy()

    ## Chunck dataset
    n_chunks = 4 # Number of chunks
    n = numpy_items.shape[0] # Number items
    nr_pairings = ncr(numpy_items.shape[0], 2) # Number of all pairings
    # Minimum nr. of pairings per chunk
    per_chunk = int(math.ceil(nr_pairings/(n_chunks - 1)))
    split_indices = [] # Array containing the indices to split at
    t = 0
    for i in range(n + 1):
        # There are n - i pairings at index i
        t += n - i
        # If the current chunk has enough pairings 
        if t >= per_chunk:
            split_indices.append(i)
            t = 0

    s_indices = [0] + split_indices + [n]
    pairings_chunks = []
    for i in range(len(s_indices)-1):
        start = s_indices[i]
        end = s_indices[i + 1]
        pairings_chunks.append(sum(range(n - end, n - start)))

    ms = numpy_items.mean(axis=1)[(slice(None,None,None),None)]
    datam = numpy_items - ms
    datass = np.sqrt(np.sum(datam*datam, axis=1))

    # Create pool of worker processes which will carry out the computation
    #must use Manager queue here, or will not work
    manager = mp.Manager()
    q = manager.Queue()
    n_cpus = mp.cpu_count()
    pool = mp.Pool(n_cpus)

    #put listener to work first
    watcher = pool.apply_async(listener, (q,concatfilepearson))

    # Start workers
    s_indices = [0] + split_indices
    workers = []
    for i in range(0, len(s_indices)):
        start = s_indices[i]
        end = s_indices[i+1] if i < len(s_indices)-1 else n-1
        workers.append(pool.apply_async(
            processor_fun_pearson, 
            kwds={'i_start': start, 'i_end': end,'outputfile' : concatfilepearson,'datam' : datam,'datass' : datass}))

    for r in workers:
        r.wait()

    # Close the pool and wait till all workers are done
    q.put('kill')
    pool.close()
    pool.join()
    pool.terminate()

    return concatfilepearson

def get_spearman_rank_corr_matrices(fulldata,inputpath,run):
    concatfilespearman = inputpath+run+'/SPEARMAN/spearman_corr_CONCAT.csv'

    file = open(concatfilespearman, "w")
    file.write('index_feat_1,index_feat_2,rank_corr_coeff\n')
    file.close()

    items = pd.read_csv(fulldata, header = 0,sep='\t')
    items = items.drop(columns=['Instance','class'])
    items = items.assign(**items.iloc[:, 0:].rank(axis = 0, ascending = False,method='average').astype(float))
    nsamples = len(items.index)
    denum = nsamples*((nsamples*nsamples)-1)
    items = items.T
    ranks = items.to_numpy()

    ## Chunck dataset
    n_chunks = 50 # Number of chunks
    n = ranks.shape[0] # Number items
    nr_pairings = ncr(ranks.shape[0], 2) # Number of all pairings
    # Minimum nr. of pairings per chunk
    per_chunk = int(math.ceil(nr_pairings/(n_chunks - 1)))
    split_indices = [] # Array containing the indices to split at
    t = 0
    for i in range(n + 1):
        # There are n - i pairings at index i
        t += n - i
        # If the current chunk has enough pairings 
        if t >= per_chunk:
            split_indices.append(i)
            t = 0

    s_indices = [0] + split_indices + [n]
    pairings_chunks = []
    for i in range(len(s_indices)-1):
        start = s_indices[i]
        end = s_indices[i + 1]
        pairings_chunks.append(sum(range(n - end, n - start)))

    # Create pool of worker processes which will carry out the computation
    #must use Manager queue here, or will not work
    manager = mp.Manager()
    q = manager.Queue()
    n_cpus = mp.cpu_count()
    pool = mp.Pool(n_cpus)

    #put listener to work first
    watcher = pool.apply_async(listener, (q,concatfilespearman))

    # Start workers
    s_indices = [0] + split_indices
    workers = []
    for i in range(0, len(s_indices)):
        start = s_indices[i]
        end = s_indices[i+1] if i < len(s_indices)-1 else n-1
        workers.append(pool.apply_async(
            processor_fun_spearman, 
            kwds={'i_start': start, 'i_end': end,'outputfile' : concatfilespearman,'denum':denum,'ranks':ranks}))

    for r in workers:
        r.wait()

    # Close the pool and wait till all workers are done
    q.put('kill')
    pool.close()
    pool.join()
    pool.terminate()

    return concatfilespearman

def get_info_gain(fulldata, inputpath, run):
    fulldata = pd.read_csv(fulldata,sep='\t',header=0,dtype=object)

    featnames = fulldata.columns.values.tolist()
    featnames.remove('Instance')
    featnames = pd.DataFrame(featnames, columns=['feature'])
    featnames['origin_index'] = range(1, len(featnames) +1)
    featnames['origin_index'] = featnames['origin_index'].apply(str)

    fulldata = fulldata.drop(columns=['Instance'])
    fulldata = fulldata.apply(pd.to_numeric, errors='ignore')

    data = fulldata.to_numpy()

    data = create_instances_from_lists(x=data)
    num_to_nom = Filter(classname="weka.filters.unsupervised.attribute.StringToNominal", options=["-R", "last"])
    num_to_nom.inputformat(data)      #data is the weka dataset whose last column is numeric.
    data=num_to_nom.filter(data)   #newData is the weka dataset whose last column is nominal.

    search = ASSearch(classname="weka.attributeSelection.Ranker", options=["-N", "-1"])
    evaluator = ASEvaluation("weka.attributeSelection.InfoGainAttributeEval")
    attsel = AttributeSelection()
    attsel.search(search)
    attsel.evaluator(evaluator)
    attsel.select_attributes(data)

    res = attsel.results_string
    resdict = {}

    m = re.search('Ranked attributes:\n(.*)\nSelected attributes:', res, re.S)
    if m:
        found = m.group(1)

    found = found.splitlines()
    for line in found :
        splittedline = line.split(' ')
        resdict[splittedline[-1]] = {}
        resdict[splittedline[-1]]['origin_index'] = splittedline[-2]
        resdict[splittedline[-1]]['info_gain_score'] = splittedline[1]

    resdf = pd.DataFrame(resdict)
    resdf = resdf.T
    resdf = resdf.merge(featnames, on = 'origin_index')

    resdf.to_csv(inputpath+run+'/INFOGAIN/INFOGAIN.csv',index=False,sep=',',header=True)

def redundancy(input,inputpath,runs):
    ## CHECK IF FOLDER PEARSON AND SPEARMAN EXIST :
    for run in runs :
        if not os.path.exists(inputpath+run+'/PEARSON') :
            print('PEARSON folder does not exist : '+run)
            os.makedirs(inputpath+run+'/PEARSON')
            print('PEARSON folder created')

        if not os.path.exists(inputpath+run+'/SPEARMAN') :
            print('SPEARMAN folder does not exist : '+run)
            os.makedirs(inputpath+run+'/SPEARMAN')
            print('SPEARMAN folder created')

        if not os.path.exists(inputpath+run+'/INFOGAIN') :
            print('INFOGAIN folder does not exist : '+run)
            os.makedirs(inputpath+run+'/INFOGAIN')
            print('INFOGAIN folder created')

        for inputfile in input :
            if '_'+run+'_' in inputfile :
                concatfilepearson = get_pearson_corr_matrices(inputfile, inputpath,run)
                map_index_to_feature_name(inputfile,concatfilepearson)  
                print('PEARSON OK : '+run)

                concatfilespearman = get_spearman_rank_corr_matrices(inputfile, inputpath,run)
                map_index_to_feature_name(inputfile,concatfilespearman)  
                print('SPEARMAN OK : '+run)

    jvm.start(system_cp=True, packages=True,system_info=True)
    for run in runs :
        for inputfile in input :
            if run in inputfile :
                get_info_gain(inputfile, inputpath, run)
    jvm.stop()

def main():
    inputpath = snakemake.params[0]
    runs = snakemake.params[1]
    nbrun = len(snakemake.params[1])
    input = snakemake.input[:nbrun]

    redundancy(input,inputpath,runs)

if __name__ == "__main__":
    main()
