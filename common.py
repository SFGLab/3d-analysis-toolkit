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
from sortedcontainers import SortedList

class Interaction:
    def __init__(self, chr1, pos1, end1, chr2, pos2, end2, pet):
        self.chr1 = chr1
        self.chr2 = chr2
        self.pos1 = pos1
        self.pos2 = pos2
        self.end1 = end1
        self.end2 = end2
        self.pet = pet
    def __eq__(self, other):
        """Overrides the default implementation"""
        if isinstance(other, Interaction):
            return self.chr1 == other.chr1 and self.pos1 == other.pos1 and self.end1 == other.end1 and self.chr2 == other.chr2 and self.pos2 == other.pos2 and self.end2 == other.end2
        return False
    def generateLine(self):
        return self.chr1+"\t"+str(self.pos1)+"\t"+str(self.end1)+"\t"+self.chr2+"\t"+str(self.pos2)+"\t"+str(self.end2)+"\t"+str(self.pet)+"\n"

def createFolder(folder):
    if os.path.exists(folder) and os.path.isdir(folder):
        shutil.rmtree(folder)
    os.mkdir(folder)
def removeFolder(folder):
    if os.path.exists(folder) and os.path.isdir(folder):
        shutil.rmtree(folder)

def loadInteractions(fileName):
    interactions = SortedList([], key=lambda x: (x.chr1, x.pos1, x.end1, x.chr2, x.pos2, x.end2))
    with open(fileName, 'r') as f: #open the file
        lines = f.readlines()
        for line in lines:
            values = line.split("\t")
            interaction = Interaction(values[0], int(values[1]), int(values[2]), values[3], int(values[4]), int(values[5]), int(values[6]))
            interactions.add(interaction)
    return interactions
def checkOverlap(interaction1, interaction2):
    extension = -1 # -1 added for actuall 1bp overlap
    if(interaction1.chr1 != interaction2.chr1 or interaction1.chr2 != interaction2.chr2):
        return False
    if(max(interaction1.pos1, interaction2.pos1)-extension <= min(interaction1.end1, interaction2.end1)):
        if(max(interaction1.pos2, interaction2.pos2)-extension <= min(interaction1.end2, interaction2.end2)):
            return True
    return False

def run_comparison(files):
    file1, file2 = files
    cmd = "cat " + file1 + " | wc -l"
    reference_count = int(subprocess.getoutput(cmd))
    cmd = "pairToPair -a "+file1+" -b "+file2+" | wc -l"
    common_count = int(subprocess.getoutput(cmd))
    return str(round(common_count/reference_count*100, 1))+"%"

def run_comparison_bed(files):
    file1, file2 = files
    cmd = "cat " + file1 + " | wc -l"
    reference_count = int(subprocess.getoutput(cmd))
    cmd = "cat " + file2 + " | wc -l"
    reference_count2 = int(subprocess.getoutput(cmd))
    reference_count = min(reference_count, reference_count2)
    cmd = "bedtools intersect -wa -a "+file1+" -b "+file2+" | wc -l"
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

def enlarge_anchors(folder_to_compare, enlargeAnchors):
    files_to_enlarge = [folder_to_compare+f for f in listdir(folder_to_compare) if isfile(join(folder_to_compare, f)) and f.split(".")[-1] == "bedpe"]
    for file in files_to_enlarge:
        interactions = loadInteractions(file)
        for interaction in interactions:
            interaction.pos1 -= enlargeAnchors
            interaction.pos2 -= enlargeAnchors
            interaction.end1 += enlargeAnchors
            interaction.end2 += enlargeAnchors
        fileName = folder_to_compare+"enlarged/"+file.split("/")[-1]
        with open(fileName, 'w') as f:
            for interaction in interactions:
                f.write(interaction.generateLine())