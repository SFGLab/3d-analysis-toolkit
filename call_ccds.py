import pandas as pd
import numpy as np
import itertools
from operator import add
import time

def get_ccds(filename, output_filename, min_loops, min_diff, min_length):
    loops = pd.read_csv(filename, sep="\t", names=["chr1", "pos1", "end1", "chr2", "pos2", "end2", "pet"])
    loops = loops.loc[loops["chr1"] == loops["chr2"]] # take only loops where chr1==chr2
    loops = loops[["chr1", "pos1", "end2"]]
    loops.columns=["chr", "start", "end"]
    loops = loops.sort_values(["chr", "start", "end"])
    ccds = list()

    for chromosome in loops["chr"].unique():
        temp = [0]*(loops.loc[loops["chr"] == chromosome]["end"].max()+1)
        loops_current_chr = list()
        for key, loop in loops.loc[loops["chr"] == chromosome].iterrows():
            loops_current_chr.append((loop["start"], loop["end"]))
        for loop in loops_current_chr:
            temp[loop[0]:loop[1]+1] = map(add, temp[loop[0]:loop[1]+1], ([1]*(loop[1]-loop[0]+1)))
        current_ccd_left = 0
        for i in range(0, len(temp)):
            diff = abs(temp[i]-temp[i-1])
            if(current_ccd_left == 0):
                if(temp[i] >= min_loops or diff >= min_diff):
                    current_ccd_left = i
            else:
                if(temp[i] < min_diff or diff > temp[i]):
                    ccds.append((chromosome, current_ccd_left, i-1))
                    current_ccd_left = 0
        break
    with open(output_filename, 'w') as file:
        for ccd in ccds:
            if(ccd[2]-ccd[1] >= min_length):
                file.write(ccd[0]+"\t"+str(ccd[1])+"\t"+str(ccd[2])+"\n")
            



filename = "hg19/GSM1872886_GM12878_CTCF_PET_clusters_fltered_cut.bedpe"

start_time_total = time.time()
get_ccds(filename, "hg19/ccds.bed", 3, 2, 10000)
print("--- Executed in %s seconds ---" % (time.time() - start_time_total))