#!/usr/bin/env python

import sys, os, re
from argparse import ArgumentParser
from collections import defaultdict
from itertools import islice

__author__ = "Xiongbin Kang, Luo Xiao"


usage = """%prog -fq <input fq file>

This program is used to cluster the reads that stem from identical species.

# need to install panda

"""

def main():
    parser = ArgumentParser(description = usage)
    parser.add_argument('-fq', dest = 'fq', type = str, help = "The input fq file.")
    parser.add_argument('-t', dest = 'threads', default = 10, type = int, help = 'The number of threads when run strainxpress. Default is 10.')
    parser.add_argument('-size', dest = 'size', default = 15000, type = int, help = "The maximum size of the cluster. Default is 15000. ")
    parser.add_argument('-insert_size', dest = 'insert_size', default = 450, type = int, help = "The length of insert size of reads. Default is 450. ")
    parser.add_argument('-average_read_len', dest = 'average_read_len', default = 250, type = int, help = "The average length of reads. Default is 250. ")
    parser.add_argument('-split_nu', dest = 'split_nu', default = 8, type = int, help = "Split the fq file into several files. Default is 8.")
    parser.add_argument('-fast', dest = 'fast', action='store_true', help = "When fq file is big or have some high coverage bacteria, suggest use the fast model. Defualt is don't use fast.")
    args = parser.parse_args()
    
    (fq, threads, size, insert_size, average_read_len, split_nu) = (args.fq, args.threads, args.size, args.insert_size, args.average_read_len, args.split_nu)
    
    if not (fq):
        print("Specify min over length and min identical.")
        parser.print_help()
        sys.exit()

    id = "strainxpress"
    folder = os.getcwd()
    bin = os.path.split(os.path.realpath(__file__))[0]  
    k = os.path.split(os.path.realpath(__file__))[0]
    nu = int(os.popen("wc -l " + fq ).readline().split()[0])
    nu_sub = int(nu/(8*split_nu)+1)*8
    folder_name = "fq_"+ str(size)

    split_line = "split {} -l {} -d -a 2 sub".format(fq, nu_sub)
    execute(split_line)

    cmd = ["minimap2 -t %s -c --sr -X -k 21 -w 11 -s 60 -m 30 -n 2 -r 0 -A 4 -B 2 --end-bonus=100 %s sub0%s 2> /dev/null | python %s/filter_trans_ovlp_inline_v3.py -len 30 -iden 0.9 -oh 1 > sub0%s.map " %(threads, fq, i, bin,i) for i in range(0, split_nu )]

    cmd2 = "\n".join(cmd)
    with open('cmd_overlap.sh', 'w') as fa:
        fa.write(cmd2)
        fa.write("\n")

    cmd_minimap = "cat cmd_overlap.sh | xargs -i -P %s bash -c \"{}\";" %threads
    execute(cmd_minimap) # run the minimap2 and get the overlap file

    if args.fast:
        execute("cat sub*.map > all_reads_sort.map")
    else:
        execute("for X in sub*.map; do sort  -k3 -nr < $X > sorted-$X; done;")
        execute("sort -k3 -nr -m sorted-sub*.map > all_reads_sort.map;")

    execute("rm *sub*;") 
    cmd_fq_name = "python %s/get_readnames.py %s readnames.txt" %(bin, fq) # get the name of reads
    execute(cmd_fq_name)

    cmd_cluster1 =  "python %s/bin_pointer_limited_filechunks_shortpath.py all_reads_sort.map readnames.txt %s %s %s" %(bin, size, id, threads)
    cmd_rm = "rm -rf Chunkfile*; rm %s_max%s_final_clustersizes.json %s_max%s_final_clusters_unchained.json %s_max%s_final_clusters.json" %(id, size, id, size, id, size)
    cmd_cluster2 = "python %s/getclusters.py %s_max%s_final %s" %(bin, id, size, threads)
    cmd_cluster3 = "python %s/get_fq_cluster.py %s_max%s_final_clusters_grouped.json %s %s/%s" %(bin, id, size, fq, folder, folder_name)
    execute(cmd_cluster1)
    execute(cmd_cluster2)
    execute(cmd_cluster3)
    execute(cmd_rm)
    
    fq_id = os.popen("ls %s/%s/" %(folder, folder_name))

    for i in fq_id:
        i = i.strip()
        fq1 = "%s/%s/%s/%s.1.fq" %(folder, folder_name, i, i)
        fq2 = "%s/%s/%s/%s.2.fq" %(folder, folder_name, i, i)
        read_folder = "%s/%s/%s" %(folder, folder_name, i)

        cmd_polyte = "cd %s; python %s/polyte.tune_params.py -p1 %s -p2 %s  -m 50 -m_EC 77 -t 1 \
        --hap_cov 10 --insert_size  %s --stddev 27  --mismatch_rate 0 --min_clique_size 2 \
        --average_read_len %s --edge1 0.93 --edge2 1.0 > %s/log.txt 2>&1" \
        %(read_folder, bin, fq1, fq2, insert_size, average_read_len, read_folder)

        with open('cmd_polyte.sh', 'a') as fa:
            fa.write(cmd_polyte)
            fa.write("\n")
    
    cmd_polyte = "cat cmd_polyte.sh | xargs -i -P %s bash -c \"{}\";" %threads
    execute(cmd_polyte)
    
    fq_id = os.popen("ls %s/%s/" %(folder, folder_name))

    for i in fq_id:
        i = i.strip()
        read_folder = "%s/%s/%s/" %(folder, folder_name, i)
        file_contigs = read_folder + "contigs.fasta"
        
        if not os.path.isfile(file_contigs):
            print("You need to rerun polyte in %s or decrease the size of cluster and rerun the whole steps" %read_folder)

    # collect all contigs from individual group.
    cmd_merge_contigs = "cat %s/%s/*/contigs.fasta > all.contigs_%s.fasta" %(folder, folder_name, size)
    execute(cmd_merge_contigs)
    
    # extend the length of contigs into master contigs
    
    #1.filter contigs with length < 150bp
    if os.path.exists('contigs_b.fastq'):
        execute("rm contigs_b.fastq")
        
    fname = "all.contigs_%s.fasta" %size
    if os.path.exists(fname):
        with open(fname, 'r') as f:
            n = 0
            for i in f:
                i = i.strip()
                if re.search(r"^>",i):
                    con_id = i
                else:
                    lenq = len(i)
                    if lenq > 150:
                        n +=1
                        seq = ("@"+str(n), i, "+", "="*lenq)
                        con_fq = "\n".join(seq)

                        with open('contigs_b.fastq', 'a') as con_fq_w:
                            con_fq_w.write(con_fq)
                            con_fq_w.write("\n")                        
                            
    contig_counts = n                      
    execute("mkdir -p stageb")
    cmd_get_overlap = "cd stageb; minimap2 -t %s --sr -X -c -k 21 -w 11 -s 60 -m 30 -n 2 -r 0 \
     -A 4 -B 2 --end-bonus=100 ../contigs_b.fastq ../contigs_b.fastq | python %s/filter_trans_ovlp_inline_v3.py \
      -len 100 -iden 0.99 -oh 2 -sfo > sfoverlaps.out;" %(threads, bin)
    execute(cmd_get_overlap)
    
    # convert minimap format to savage format
    cmd_convert = "cd stageb; python %s/sfo2overlaps.py --in sfoverlaps.out \
    --out sfoverlap.out.savage --num_singles %s --num_pairs 0; mkdir -p fastq;\
    cp ../contigs_b.fastq ./fastq/singles.fastq;" %(bin, contig_counts)
    execute(cmd_convert)
    
    cmd_extend = "cd stageb; python %s/pipeline_per_stage.v3.py --no_error_correction --remove_branches true \
     --stage b --min_overlap_len 100 --min_overlap_perc 0 --edge_threshold 1 --overlaps ./sfoverlap.out.savage \
     --fastq ./fastq --max_tip_len 1000 --num_threads %s; python %s/fastq2fasta.py ./singles.fastq \
     ./final_contigs.fasta; rm -r contigs_b.fasta  fastq graph* p* s* removed_tip_sequences.fastq;" %(bin, threads, bin)
    
    execute(cmd_extend)
    
    

def execute(cmd):
    te = os.system(cmd + " 2>output.txt")
    if te:
        with open("output.txt","r") as file:
            print("Don't execute the command: %s " %cmd, end='')
            print(file.read())
    else:
        print("successfully execute: %s" %cmd)
        



if __name__ == '__main__':
        sys.exit(main())
