# import packages and files
import json
import numpy as np
import time
import multiprocessing as mp
import collections
import sys
import os

os.system("taskset -p 0xff %d" % os.getpid())
splits = 60

filename = sys.argv[1]
threads = int(sys.argv[2])

clustlist = json.load(open(filename+'_clusters.json','r'))
clustsizes = json.load(open(filename+'_clustersizes.json','r'))

def findhead_nochanges(r):
    stop = 0
    rcurrent = r
    while stop == 0:
        if isinstance(clustlist[rcurrent], int):
            clust = clustlist[rcurrent]
            stop = 1
        else:
            rcurrent = clustlist[rcurrent]
    return (rcurrent,clust)

def clustering(d):
    clusters_sub = {}
    for val in np.unique(list(d.values())):
        clusters_sub[str(val)] = [k for k,v in list(d.items()) if v==val]
    return clusters_sub

# the number of reads in the same cluster
clusters = {}
nreads = len(list(clustlist.keys()))
starttime = time.time()
count = 0
p_nreads = nreads/100
percent = 1
for key in list(clustlist.keys()):
    rcurrent, clust = findhead_nochanges(key)
    clustlist[key] = clust
    count += 1
    if count > percent*p_nreads:
        percent += 1
json.dump(clustlist,open(filename+'_clusters_unchained.json','w'))

# split into small group
starttime = time.time()
clustlist_large = dict((k, v) for k, v in list(clustlist.items()) if clustsizes[str(v)] >= 20)
clustlist = clustlist_large
clustlist_large = None
starttime = time.time()
dictsize = int(len(list(clustlist.items()))/threads)
dcl = collections.defaultdict(dict)
for i in range(threads):
    start_idx = i*dictsize
    end_idx = min((i+1)*dictsize,len(list(clustlist.items())))
    dcl[i] = dict(list(clustlist.items())[start_idx:end_idx])

# extract the chain according the value of dict
if __name__ == "__main__":
    starttime = time.time()
    maplist = []
    po = mp.Pool(threads)
    result = po.map(clustering,[dcl[0],dcl[1],dcl[2],dcl[3],dcl[4],dcl[5],dcl[6],dcl[7],dcl[8],dcl[9],dcl[10],dcl[11],dcl[12],dcl[13],dcl[14],dcl[15],dcl[16],dcl[17],dcl[18],dcl[19],dcl[20],dcl[21],dcl[22],dcl[23],dcl[24],dcl[25],dcl[26],dcl[27],dcl[28],dcl[29],dcl[30],dcl[31],dcl[32],dcl[33],dcl[34],dcl[35],dcl[36],dcl[37],dcl[38],dcl[39],dcl[40],dcl[41],dcl[42],dcl[43],dcl[44],dcl[45],dcl[46],dcl[47],dcl[48],dcl[49],dcl[50],dcl[51],dcl[52],dcl[53],dcl[54],dcl[55],dcl[56],dcl[57],dcl[58],dcl[59]])
    po.close()
    po.join()
    
# merge the result of same chain and then write them
    d = result[0]
    for i in range(1,len(result)):
        for k,v in list(result[i].items()):
            if k in list(d.keys()):
                d[k].extend(v)
            else:
                d[k] = v
    json.dump(d,open(filename+'_clusters_grouped.json','w'))

