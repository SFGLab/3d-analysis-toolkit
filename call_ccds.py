import pandas as pd
import numpy as np
import itertools
import argparse
import os
import sys
from operator import add, itemgetter
import time
import multiprocessing as mp

min_huge_ccd = 25_000_000

def check_huge_ccds(results):
    if (len(results) == 0):
        return True
    for chrom, result in results:
        for ccd in result:
            if(ccd[1]-ccd[0] >= min_huge_ccd):
                return True
    return False

def get_ccds(filename, output_filename, min_loops, min_diff, min_length):
    loops = pd.read_csv(filename, sep="\t", header=None).iloc[:,0:7]
    loops.columns=["chr1", "pos1", "end1", "chr2", "pos2", "end2", "pet"]
    
    loops = loops.loc[loops["chr1"] == loops["chr2"]].loc[loops["end2"]-loops["pos1"] < min_huge_ccd] # take only loops where chr1==chr2
    loops = loops[["chr1", "pos1", "end2"]]
    loops.columns=["chr", "start", "end"]
    loops = loops.sort_values(["chr", "start", "end"])
    results = list()
    if(len(loops) > 0):
        while(check_huge_ccds(results)):
            results = list()
            loops_parsed = dict()
            for chromosome in loops["chr"].unique():
                loops_current_chr = list()
                for key, loop in loops.loc[loops["chr"] == chromosome].iterrows():
                    loops_current_chr.append((loop["start"], loop["end"]))
                loops_parsed[chromosome] = loops_current_chr

            task_list = [(loops, chromosome, min_loops, min_diff) for chromosome, loops in loops_parsed.items()]

            threads = min(len(task_list), 24)
            pool = mp.Pool(threads)
            results = pool.map(parse_chrom, task_list)
            pool.close()
            pool.join()

            print("Calculated for: " + str(min_loops), flush=True) 
            min_loops += 1
    else:
        results = list()

    with open(output_filename, 'w') as file:
        for chrom, result in results:
            for ccd in result:
                if(ccd[1]-ccd[0] >= min_length):
                    file.write(chrom+"\t"+str(ccd[0])+"\t"+str(ccd[1])+"\n")

def parse_chrom(args):
    loops_current_chr, chromosome, min_loops, min_diff = args

    ccds_current_chr = list()

    temp = [0]*(max(loops_current_chr, key=itemgetter(1))[1]+1)

    for loop in loops_current_chr:
        temp[loop[0]:loop[1]+1] = map(add, temp[loop[0]:loop[1]+1], ([1]*(loop[1]-loop[0]+1)))
    current_ccd_left = 0
    for i in range(0, len(temp)):
        diff = abs(temp[i]-temp[i-1])
        if(current_ccd_left == 0):
            if(temp[i] >= min_loops or diff >= min_diff):
                current_ccd_left = i
        else:
            if(temp[i] < min_loops or diff > temp[i]):
                ccds_current_chr.append((current_ccd_left, i-1))
                current_ccd_left = 0

    if(current_ccd_left > 0):
        ccds_current_chr.append((current_ccd_left, len(temp)-1))
        current_ccd_left = 0

    return (chromosome, ccds_current_chr)

def main():
    parser = argparse.ArgumentParser(prog='CCD caller')
    parser.add_argument('input_file', help='Input bedpe file.')
    parser.add_argument('-o', '--output_file', default='ccds.bed', help='Output file.')
    parser.add_argument('--min_loops', type=int, default=3)
    parser.add_argument('--min_diff', type=int, default=2)
    parser.add_argument('--min_length', type=int, default=10_000)
    args = parser.parse_args()

    if not os.path.exists(args.input_file):
        print("Input does not exist")
        sys.exit(1)

    start_time_total = time.time()        
    get_ccds(
        args.input_file, args.output_file,
        args.min_loops, args.min_diff, args.min_length
    )
    print("--- Executed in %s seconds ---" % (time.time() - start_time_total))

if __name__ == "__main__":
    main()
