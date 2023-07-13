import numpy as np
import csv
import json
import time
import sys
import os
import multiprocessing
import math
from multiprocessing import Value
import collections
import os
import gzip
from pathlib import Path

#print 'Run this script in Python 2'
starttime = time.time()
ovlpfile = sys.argv[1] #'/export/scratch2/xiongbin/project/meta_genome/cami/desman/2.orger_illuma_20200430/test_ovlp.txt' 
readnamefile = sys.argv[2] #'/export/scratch2/xiongbin/project/meta_genome/cami/desman/2.orger_illuma_20200430/readnames.txt'
maxsize = sys.argv[3] #33000 #float('inf')
rID = sys.argv[4]
run_ID = rID+'_max'+maxsize
threads = int(sys.argv[5])

if maxsize == 'inf':
    maxsize = float('inf')
else:
    maxsize = int(maxsize)
chunksize = 2600000 #This is approximately the number of bytes used for 100000 lines


def chunkify(fname):
    fileEnd = os.path.getsize(fname)
    with gzip.open(fname,'rb') as f: #must be binary rb
        chunkEnd = f.tell()
        while True:
            chunkStart = chunkEnd
            f.seek(chunksize,1)
            l = f.readline()
            chunkEnd = f.tell()
            yield chunkStart, chunkEnd - chunkStart
            if not l: break

# returen the keys and values of dict (read names and number)
def findhead_nochanges(r):
    stop = 0
    rcurrent = r
    while stop == 0:
        if isinstance(clusters[rcurrent], int):
            clust = clusters[rcurrent]
            stop = 1
        else:
            rcurrent = clusters[rcurrent]
    return (rcurrent,clust)

def findhead_lim(r):
    stop = 0
    rcurrent = r
    pathlen = 1
    adapt = []
    while stop == 0:
        if isinstance(clusters[rcurrent], int):
            clust = clusters[rcurrent]
#            if pathlen >= 3:
#                for node in adapt[:-1]:
#                    clusters[node] = rcurrent
            stop = 1
        else:
            adapt.append(rcurrent)
            rcurrent = clusters[rcurrent]
            pathlen += 1
    return (rcurrent,clust,pathlen)

def clusteralgorithm(fileID):#(fileID,clusters,clusterlist):
    linecount = 0
    with open(fileID,'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            matches = row[1]
            cllim_r1_head, clust1, pathlen1 = findhead_lim(row[0][:-2])
            cllim_r2_head, clust2, pathlen2 = findhead_lim(row[1][:-2])
            if clust1 != clust2:
                clustsize = clusterlist[clust1]+clusterlist[clust2]
                if clustsize <= maxsize:
                    if pathlen2 < pathlen1:
                        clusters[cllim_r2_head] = cllim_r1_head
                        clusterlist[clust1] = clustsize
                        clusterlist.pop(clust2)
                    else:
                        clusters[cllim_r1_head] = cllim_r2_head
                        clusterlist[clust2] = clustsize
                        clusterlist.pop(clust1)
            linecount += 1
    return linecount, matches

def getchunkfile(a):#(chunknum_,chunkstart_,chunksize_):
    chunknum_ = a[0]
    chunkstart_ = a[1]
    chunksize_ = a[2]
    with open('Chunkfile_'+run_ID+'_'+str(chunknum_),'w') as fw:
        with gzip.open(ovlpfile,'rt') as fo: #Must be binary?
            fo.seek(chunkstart_)
            lines = fo.read(chunksize_).splitlines()
            for line in lines:
                linesplit = line.rstrip().split('\t')
                cllim_r1_head, clust1 = findhead_nochanges(linesplit[0][:-2])
                cllim_r2_head, clust2 = findhead_nochanges(linesplit[1][:-2])
                if clust1 != clust2 and clusterlist[clust1]+clusterlist[clust2] < maxsize:
                    fw.write(line+'\n')
    return

nmatches = 150
if __name__ == '__main__':
    clusters = {}
    clusterlist = {}
    i = 1
    with open(readnamefile,'r') as f:
        for line in f:
            clusters[line.rstrip()] = i
            clusterlist[i] = 1
            i += 1

# Make chunks and sessions
    starttime = time.time()
    sessiondict = {}
    session = 1
    chunknum = 1
    sessionlist = []
    for chunkStart,chunkSize in chunkify(ovlpfile):
        sessiondict[session]=[session,chunkStart,chunkSize]
        session += 1
        
# print('Chunks made' )
    numsessions = len(list(sessiondict.keys()))
    starttime = time.time()
    sess = 1
    
    while sess <= numsessions:
        po = multiprocessing.Pool(threads)
        for th in range(threads):
            nfiles = 0
            if sess <= numsessions:
                session = sessiondict[sess]
                sess += 1
                nfiles += 1
                po.apply_async(getchunkfile, args=([th+1,session[1],session[2]],))
                
        po.close()
        po.join()
 
# merge chunks
        with open('Chunkfile_'+run_ID,'w') as fcat:
            for ch in range(1,threads+1):
                my_file = Path('Chunkfile_'+run_ID+'_'+str(ch))
                if my_file.exists():
                    with open('Chunkfile_'+run_ID+'_'+str(ch),'r') as infile:
                        fcat.write(infile.read())
# cluster the reads                   
        if os.stat('Chunkfile_'+run_ID).st_size != 0:

            linecount, matches = clusteralgorithm('Chunkfile_'+run_ID)
            if sess % threads == 100:
                with gzip.open(ovlpfile,'rt') as fo:
                    fo.seek(sessiondict[sess][1])
                    for line in fo:
                        ln = line.split('\t')
                        threshold = ln[2].rstrip()
                        break
                with open(run_ID+'_round'+str(sess)+'_threshold'+str(threshold)+'_clusters.json','w') as fp:
                    json.dump(clusters,fp)
                with open(run_ID+'_round'+str(sess)+'_threshold'+str(threshold)+'_clustersizes.json','w') as fp:
                    json.dump(clusterlist,fp)

# remove redundant file
        my_file = Path('Chunkfile_'+run_ID)
        if my_file.exists():
            os.remove('Chunkfile_'+run_ID)
        for ch in range(1,nfiles+1):
            my_file = Path('Chunkfile_'+run_ID+'_'+str(ch))
            if my_file.exists():
                os.remove('Chunkfile_'+run_ID+'_'+str(ch))
        nextline = (sess-1)*100000

    with open(run_ID+'_final_clusters.json','w') as fp:
        json.dump(clusters,fp)
    with open(run_ID+'_final_clustersizes.json','w') as fp:
        json.dump(clusterlist,fp)


