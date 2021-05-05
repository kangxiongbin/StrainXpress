#!/usr/bin/env python

import json
import re
import random
import sys,os
import pandas as pd

clusters_groupe = sys.argv[1] # CAMI_mine_low_max15000_final_clusters_grouped2.json
fq = sys.argv[2] # E_coli_alb.fq
outdir = sys.argv[3] # output folder 
os.system("mkdir -p "+outdir)

cluster2read =json.load(open(clusters_groupe,'r'))
read_count_list = [len(i) for i in cluster2read.values()]
cluster2size={key:len(val) for key,val in cluster2read.items()}
read2cluster = {read:key for key,val in cluster2read.items() for read in val}
cluster2size_se=pd.Series(cluster2size)
tar_clusters=list(cluster2read.keys())
cluster_len=len(tar_clusters)

def get_fq4cluster(tar_clusters,file,outdir):
    for cluster in tar_clusters:
        outdir4c=outdir + "/" + cluster + "/"
        os.system("mkdir -p " + outdir4c)
        locals()["fq1w%s"%cluster] = open(outdir4c + cluster + '.1.fq', 'w')
        locals()["fq2w%s"%cluster] = open(outdir4c + cluster + '.2.fq', 'w')
    #
    pre_record = []
    pre_read = ''
    pre_header = ''
    with open(file, 'r') as fq:

        for i, line in enumerate(fq):
            if i%1000000 == 0: print('this is the: '+str(i)+' for 100w lines')
            if i % 4 == 0:
                if (i and (pre_read in read2cluster)) and (read2cluster[pre_read] in tar_clusters):
                    if re.search(r'/1$', pre_header):
                        locals()["fq1w%s"%read2cluster[pre_read]].writelines(pre_record)
                    else:
                        locals()["fq2w%s"%read2cluster[pre_read]].writelines(pre_record)
                pre_record = []
                pre_record.append(line)
                pre_read = re.split('[@/]', line)[-2]
                pre_header = line
            else:
                pre_record.append(line)
                
        if (i and (pre_read in read2cluster)) and (read2cluster[pre_read] in tar_clusters):
            if re.search(r'/1$', pre_header):
                locals()["fq1w%s"%read2cluster[pre_read]].writelines(pre_record)
            else:
                locals()["fq2w%s"%read2cluster[pre_read]].writelines(pre_record)
    fq.close()
    for cluster in tar_clusters:
        locals()["fq1w%s"%cluster].close()
        locals()["fq2w%s"%cluster].close()

print("begin...")
print("#"*50)

k=420
split_n = int(cluster_len/k)+1
for i in range(split_n):
    print("the "+str(i+1)+"/"+str(split_n)+" part start...")
    start=k*i
    end=k*(i+1) if (k*(i+1))<cluster_len else  cluster_len
    cc=tar_clusters[start:end]
    get_fq4cluster(cc,fq,outdir)
    print("the "+str(i+1)+"/"+str(split_n)+" part finished...\n")
    print("#"*50)
