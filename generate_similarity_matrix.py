# script for generating histograms for anchor & loop (interaction) lengths
import subprocess
import numpy as np
import math
import pandas as pd
import shutil
import shlex
import os
import time
from matplotlib import pyplot as plt
from os import listdir
from os.path import isfile, join
from collections import defaultdict
import multiprocessing as mp


def run_comparison(files):
    file1, file2 = files
    cmd = "cat " + file1 + " | wc -l"
    reference_count = int(subprocess.getoutput(cmd))
    cmd = "pairToPair -a "+file1+" -b "+file2+" | wc -l"
    common_count = int(subprocess.getoutput(cmd))
    return str(round(common_count/reference_count*100, 1))+"%"

def get_counts(file):
    cmd = "cat " + file + " | wc -l"
    reference_count = int(subprocess.getoutput(cmd))
    return str(reference_count)

def create_loops(file, folder, peaks=False):
    peaks_line = ""
    if(peaks):
        peaks_line = "--peaks_filename " + os.path.splitext(file)[0] + ".bed"
    fileName = folder+file.split("/")[-1]
    loop_command = '/home/mateuszchilinski/.pyenv/shims/python /mnt/raid/ctcf_prediction_anal/cluster-paired-end-tags/cluster_pets/cluster_PETs.py --pets_filename '+file+' '+peaks_line+' --clusters_filename '+fileName
    if(peaks):
        loop_command += "_2"
    output = subprocess.getoutput(loop_command)
    if(peaks):
        output = subprocess.getoutput("cut -f1-7 " + fileName+"_2 > " + fileName)
        os.remove(fileName+"_2")



def generate_matrix(folder_to_compare):
    threads = 16
    files_to_compare = [folder_to_compare+f for f in listdir(folder_to_compare) if isfile(join(folder_to_compare, f)) and f.split(".")[-1] == "bedpe"]
    matrix = defaultdict(dict)
    task_list = list()
    for file1 in files_to_compare:
        for file2 in files_to_compare:
            task_list.append((file1, file2))
    pool = mp.Pool(threads)
    results = pool.map(run_comparison, task_list)
    pool.close()
    pool.join() 
    i = 0
    for result in results:
        file1, file2 = task_list[i]
        matrix[file1.split("/")[-1].split(".")[0]][file2.split("/")[-1].split(".")[0]] = result
        i += 1

    df = pd.DataFrame(matrix).T
    df = pd.concat(
        [pd.concat(
            [df],
            keys=['File'], axis=1)],
        keys=['Reference']
    )

    df = df.sort_index().sort_index(axis = 1)
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
        print(df)

    print("")
    print("Counts of interactions in files:")

    for file in files_to_compare:
        print(file.split("/")[-1] + ": " + get_counts(file))

start_time = time.time()

folder_to_compare = '/mnt/raid/ctcf_prediction_anal/trios_new_ctcf/ctcf_named/output/'
#folder_to_compare = '/mnt/raid/ctcf_prediction_anal/trios_new_ctcf/ctcf_named/output/'
print("===== INTERACTIONS =====")
generate_matrix(folder_to_compare)

if os.path.exists(folder_to_compare+"temp") and os.path.isdir(folder_to_compare+"temp"):
    shutil.rmtree(folder_to_compare+"temp")
os.mkdir(folder_to_compare+"temp")

if os.path.exists(folder_to_compare+"temp2") and os.path.isdir(folder_to_compare+"temp2"):
    shutil.rmtree(folder_to_compare+"temp2")
os.mkdir(folder_to_compare+"temp2")

files_to_compare = [folder_to_compare+f for f in listdir(folder_to_compare) if isfile(join(folder_to_compare, f)) and f.split(".")[-1] == "bedpe"]
for file in files_to_compare:
    create_loops(file, folder_to_compare+"temp/")
print("===== LOOPS (NO PEAKS) =====")
generate_matrix(folder_to_compare+"temp/")
for file in files_to_compare:
    if(os.path.isfile(os.path.splitext(file)[0]+".bed")): 
        create_loops(file, folder_to_compare+"temp2/", True)
print("===== LOOPS (PEAKS) =====")
generate_matrix(folder_to_compare+"temp2/")

print("--- Executed in %s seconds ---" % (time.time() - start_time))

#shutil.rmtree(folder_to_compare+"temp")