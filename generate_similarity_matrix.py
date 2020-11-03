# script for generating histograms for anchor & loop (interaction) lengths
import subprocess
import numpy as np
import math
import pandas as pd
import shutil
import shlex
import os
from matplotlib import pyplot as plt
from os import listdir
from os.path import isfile, join
from collections import defaultdict
def run_comparison(file1, file2):
    cmd = "cat " + file1 + " | wc -l"
    reference_count = int(subprocess.getoutput(cmd))
    cmd = "pairToPair -a "+file1+" -b "+file2+" | wc -l"
    common_count = int(subprocess.getoutput(cmd))
    return str(round(common_count/reference_count*100, 1))+"%"

def get_counts(file):
    cmd = "cat " + file + " | wc -l"
    reference_count = int(subprocess.getoutput(cmd))
    return str(reference_count)

def create_loops(file, folder):
    loop_command = '/home/mateuszchilinski/.pyenv/shims/python /mnt/raid/ctcf_prediction_anal/cluster-paired-end-tags/cluster_pets/cluster_PETs.py --pets_filename '+file+' --clusters_filename '+folder+file.split("/")[-1]
    output = subprocess.getoutput(loop_command)


def generate_matrix(folder_to_compare):

    files_to_compare = [folder_to_compare+f for f in listdir(folder_to_compare) if isfile(join(folder_to_compare, f))]
    matrix = defaultdict(dict)
    for file1 in files_to_compare:
        for file2 in files_to_compare:
            matrix[file1.split("/")[-1]][file2.split("/")[-1]] = run_comparison(file1, file2)

    df = pd.DataFrame(matrix).T
    df = pd.concat(
        [pd.concat(
            [df],
            keys=['File'], axis=1)],
        keys=['Reference']
    )

    df = df.sort_index().sort_index(axis = 1)
    print(df)

    print("")
    print("Counts of interactions in files:")

    for file in files_to_compare:
        print(file.split("/")[-1] + ": " + get_counts(file))

folder_to_compare = '/mnt/raid/ctcf_prediction_anal/trios_new_ctcf/ctcf_named_2/'
generate_matrix(folder_to_compare)

if os.path.exists(folder_to_compare+"temp") and os.path.isdir(folder_to_compare+"temp"):
    shutil.rmtree(folder_to_compare+"temp")
os.mkdir(folder_to_compare+"temp")

files_to_compare = [folder_to_compare+f for f in listdir(folder_to_compare) if isfile(join(folder_to_compare, f))]
for file in files_to_compare:
    create_loops(file, folder_to_compare+"temp/")
generate_matrix(folder_to_compare+"temp/")

#shutil.rmtree(folder_to_compare+"temp")