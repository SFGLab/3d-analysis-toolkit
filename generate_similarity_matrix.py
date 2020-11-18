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
from multiprocessing.pool import ThreadPool
from common import Interaction, createFolder, removeFolder, loadInteractions, run_comparison, run_comparison_bed, get_counts, create_loops, enlarge_anchors, saveFile, removeOverlapping
from find_motifs import filterInteractionsByMotifs

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
    pool = ThreadPool(threads)
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
    count_of_what = "interactions"
    if(ext=="bed"):
        count_of_what = "peaks"
    print("Counts of "+count_of_what+" in files:")

    for file in files_to_compare:
        print(file.split("/")[-1] + ": " + get_counts(file))

def createRandomSample(fileName, interactions, size=150000):
    interactions.sort(reverse=True, key=lambda x: x.pet)
    lastPet = interactions[size].pet
    sureInteractions = [interaction for interaction in interactions if interaction.pet > lastPet]
    if(len(sureInteractions) < size):
        toRandomInteractions = [interaction for interaction in interactions if interaction.pet == lastPet]
        chosenInteractions = np.random.choice(toRandomInteractions,size-len(sureInteractions),replace=False)
        sureInteractions = np.append(chosenInteractions, sureInteractions)
    return (fileName, sureInteractions)

def createRandomSampleFile(file, folder, size=150000):
    fileName = folder+file.split("/")[-1]
    interactions = loadInteractions(file)
    sureInteractions = createRandomSample(fileName, interactions, size)
    saveFile(fileName, sureInteractions)
    return

def filterInteractionsByMotifsParallel(args):
    fileName, interactions = args
    return (fileName, filterInteractionsByMotifs(interactions))
start_time_total = time.time()

#folder_to_compare = '/mnt/raid/ctcf_prediction_anal/GM_comparisons_tries/'
randomSampling = False
filterMotifs = True
enlargeAnchors = 500 # 0 = disabled
maxLength = 500000
folder_to_compare = '/mnt/raid/ctcf_prediction_anal/trios_new_ctcf/ctcf_named/output/'
rs_temp = ""

print("===== PEAKS =====")

generate_matrix(folder_to_compare,0,run_comparison_bed, "bed")
if randomSampling or filterMotifs or maxLength > 0:
    rs_temp = "modified_temp/"
    createFolder(folder_to_compare+rs_temp)
    files_to_copy = [folder_to_compare+f for f in listdir(folder_to_compare) if isfile(join(folder_to_compare, f)) and f.split(".")[-1] == "bed"]
    for file in files_to_copy:
        shutil.copy2(file, folder_to_compare+rs_temp+file.split("/")[-1])

    files_to_compare = [folder_to_compare+f for f in listdir(folder_to_compare) if isfile(join(folder_to_compare, f)) and f.split(".")[-1] == "bedpe"]

    samples = dict()

    for file in files_to_compare:
        samples[folder_to_compare+rs_temp+file.split("/")[-1]] = loadInteractions(file, maxLength=maxLength)

if randomSampling:
    print("===== SAMPLING IS ON =====")
    samples_2 = dict()
    for sample, interactions in samples.items():
        samples_2[sample] = createRandomSample(interactions, 200000)
    samples = samples_2

if filterMotifs:
    start_time = time.time()
    print("===== FILTERING MOTIFS IS ON =====")
    
    samples_2 = dict()
    #for sample, interactions in samples.items():
    #    samples_2[sample] = filterInteractionsByMotifs(interactions)

    
    task_list = [(sample, interactions) for sample, interactions in samples.items()]
    threads = 16
    pool = ThreadPool(threads)
    results = pool.map(filterInteractionsByMotifsParallel, task_list)
    pool.close()
    pool.join() 

    for result in results:
        sample, interactions = result
        samples_2[sample] = interactions
    samples = samples_2
    print("--- Executed in %s seconds ---" % (time.time() - start_time))


if randomSampling or filterMotifs or maxLength > 0:
    for sample, interactions in samples.items():
        saveFile(sample, interactions)

print("===== INTERACTIONS =====")
generate_matrix(folder_to_compare+rs_temp, enlargeAnchors)

createFolder(folder_to_compare+rs_temp+"temp")
createFolder(folder_to_compare+rs_temp+"temp2")

files_to_compare = [folder_to_compare+rs_temp+f for f in listdir(folder_to_compare+rs_temp) if isfile(join(folder_to_compare+rs_temp, f)) and f.split(".")[-1] == "bedpe"]

for file in files_to_compare:
    create_loops(file, folder_to_compare+rs_temp+"temp/", False)
    if(randomSampling):
        createRandomSampleFile(file, folder_to_compare+rs_temp+"temp/", 20000)
print("===== LOOPS (NO PEAKS) =====")
generate_matrix(folder_to_compare+rs_temp+"temp/", enlargeAnchors)

for file in files_to_compare:
    if(os.path.isfile(os.path.splitext(file)[0]+".bed")): 
        create_loops(file, folder_to_compare+rs_temp+"temp2/", True)
        if(randomSampling):
            createRandomSampleFile(file, folder_to_compare+rs_temp+"temp2/", 10000)
print("===== LOOPS (PEAKS) =====")
generate_matrix(folder_to_compare+rs_temp+"temp2/", enlargeAnchors)



print("--- Executed in %s seconds ---" % (time.time() - start_time_total))

removeFolder(folder_to_compare+rs_temp+"temp")
removeFolder(folder_to_compare+rs_temp+"temp2")

if randomSampling:
    removeFolder(folder_to_compare+rs_temp)