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
from common import Interaction, createFolder, removeFolder, loadInteractions, run_comparison, run_comparison_bed, get_counts, create_loops, enlarge_anchors

def generate_matrix(folder_to_compare, enlargeAnchors=0, func_to_use=run_comparison, ext="bedpe"):
    if enlargeAnchors > 0:
        createFolder(folder_to_compare+"enlarged/")
        enlarge_anchors(folder_to_compare, enlargeAnchors)
        folder_to_compare += "enlarged/"
    threads = 16
    files_to_compare = [folder_to_compare+f for f in listdir(folder_to_compare) if isfile(join(folder_to_compare, f)) and f.split(".")[-1] == ext]
    matrix = defaultdict(dict)
    task_list = list()
    for file1 in files_to_compare:
        for file2 in files_to_compare:
            task_list.append((file1, file2))
    pool = mp.Pool(threads)
    results = pool.map(func_to_use, task_list)
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

def createRandomSample(file, folder, size=150000):
    fileName = folder+file.split("/")[-1]
    interactions = loadInteractions(file)
    interactions.sort(reverse=True, key=lambda x: x.pet)
    lastPet = interactions[size].pet
    sureInteractions = [interaction for interaction in interactions if interaction.pet > lastPet]
    if(len(sureInteractions) < size):
        toRandomInteractions = [interaction for interaction in interactions if interaction.pet == lastPet]
        chosenInteractions = np.random.choice(toRandomInteractions,size-len(sureInteractions),replace=False)
        sureInteractions = np.append(chosenInteractions, sureInteractions)
    with open(fileName, 'w') as f:
        for interaction in sureInteractions:
            f.write(interaction.generateLine())
    return

start_time = time.time()

folder_to_compare = '/mnt/raid/ctcf_prediction_anal/GM_comparisons_tries/'
randomSampling = False
enlargeAnchors = 500 # 0 = disabled
#folder_to_compare = '/mnt/raid/ctcf_prediction_anal/trios_new_ctcf/ctcf_named/output/'
rs_temp = ""

print("===== PEAKS =====")

generate_matrix(folder_to_compare,0,run_comparison_bed, "bed")

if randomSampling:
    rs_temp = "rs_temp/"
    print("===== SAMPLING IS ON =====")
    createFolder(folder_to_compare+rs_temp)
    files_to_copy = [folder_to_compare+f for f in listdir(folder_to_compare) if isfile(join(folder_to_compare, f)) and f.split(".")[-1] == "bed"]
    for file in files_to_copy:
        shutil.copy2(file, folder_to_compare+rs_temp+file.split("/")[-1])

    files_to_compare = [folder_to_compare+f for f in listdir(folder_to_compare) if isfile(join(folder_to_compare, f)) and f.split(".")[-1] == "bedpe"]
    for file in files_to_compare:
        createRandomSample(file, folder_to_compare+rs_temp, 200000)
print("===== INTERACTIONS =====")
generate_matrix(folder_to_compare+rs_temp, enlargeAnchors)

createFolder(folder_to_compare+rs_temp+"temp")
createFolder(folder_to_compare+rs_temp+"temp2")

files_to_compare = [folder_to_compare+rs_temp+f for f in listdir(folder_to_compare+rs_temp) if isfile(join(folder_to_compare+rs_temp, f)) and f.split(".")[-1] == "bedpe"]

for file in files_to_compare:
    create_loops(file, folder_to_compare+rs_temp+"temp/", False)
    if(randomSampling):
        createRandomSample(file, folder_to_compare+rs_temp+"temp/", 20000)
print("===== LOOPS (NO PEAKS) =====")
generate_matrix(folder_to_compare+rs_temp+"temp/", enlargeAnchors)
for file in files_to_compare:
    if(os.path.isfile(os.path.splitext(file)[0]+".bed")): 
        create_loops(file, folder_to_compare+rs_temp+"temp2/", True)
        if(randomSampling):
            createRandomSample(file, folder_to_compare+rs_temp+"temp2/", 10000)
print("===== LOOPS (PEAKS) =====")
generate_matrix(folder_to_compare+rs_temp+"temp2/", enlargeAnchors)

print("--- Executed in %s seconds ---" % (time.time() - start_time))

removeFolder(folder_to_compare+rs_temp+"temp")
removeFolder(folder_to_compare+rs_temp+"temp2")

if randomSampling:
    removeFolder(folder_to_compare+rs_temp)