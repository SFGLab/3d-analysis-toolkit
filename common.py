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
from sortedcontainers import SortedList

class Peak:
    def __init__(self, chr, pos, end, name, score, strand, signalValue, pValue, qValue):
        self.chr = chr
        self.pos = pos
        self.end = end
        self.name = name
        self.score = score
        self.strand = strand
        self.signalValue = signalValue
        self.pValue = pValue # only this field is usually set
        self.qValue = qValue
    def __eq__(self, other):
        """Overrides the default implementation"""
        if isinstance(other, Interaction):
            return self.chr == other.chr and self.pos == other.pos and self.end == other.end
        return False
    def generateLine(self, add_new_line=True):
        return self.chr+"\t"+str(self.pos)+"\t"+str(self.end)+"\t"+self.name+"\t"+str(self.score)+"\t"+self.strand+"\t"+str(self.signalValue)+"\t"+str(self.pValue)+"\t"+str(self.qValue)+"\n"
    def __hash__(self):
        return hash((self.chr, self.pos, self.end))

class Interaction:
    def __init__(self, chr1, pos1, end1, chr2, pos2, end2, pet, samples=""):
        self.chr1 = chr1
        self.chr2 = chr2
        self.pos1 = pos1
        self.pos2 = pos2
        self.end1 = end1
        self.end2 = end2
        self.pet = pet
        self.samples = list()
        if(samples != ""):
            self.samples.extend(samples)
    def __eq__(self, other):
        """Overrides the default implementation"""
        if isinstance(other, Interaction):
            return self.chr1 == other.chr1 and self.pos1 == other.pos1 and self.end1 == other.end1 and self.chr2 == other.chr2 and self.pos2 == other.pos2 and self.end2 == other.end2
        return False
    def generateLine(self, add_sample_names=False, add_new_line=True):
        sample_names = ""
        if(add_sample_names):
            sample_names = ','.join(self.samples)
        line = self.chr1+"\t"+str(self.pos1)+"\t"+str(self.end1)+"\t"+self.chr2+"\t"+str(self.pos2)+"\t"+str(self.end2)+"\t"+str(self.pet)
        if(add_sample_names):
            line +="\t"+sample_names
        if(add_new_line):
            line += "\n"
        return line
    def __hash__(self):
        return hash((self.chr1, self.chr2, self.pos1, self.pos2, self.end1, self.end2))

class SV:
    def __init__(self, chrm, pos, end, svlen, svtype, line="", samples=""):
        if("chr" in chrm):
            self.chr = chrm
        else:
            self.chr = "chr"+chrm
        self.pos = pos
        self.svlen = svlen
        self.end = end
        self.svtype = svtype
        self.line = line
        if(samples):
            self.samples = samples.split(",")
    def __hash__(self):
        return hash((self.chr, self.pos, self.end, self.svtype))
    def __eq__(self, other):
        """Overrides the default implementation"""
        if isinstance(other, Interaction):
            return self.chr == other.chr and self.pos == other.pos and self.end == other.end and sv.svtype == other.svtype
        return False
    def getGTs(self):
        return self.samples

def createFolder(folder):
    if os.path.exists(folder) and os.path.isdir(folder):
        shutil.rmtree(folder)
    os.mkdir(folder)
def removeFolder(folder):
    if os.path.exists(folder) and os.path.isdir(folder):
        shutil.rmtree(folder)

def orderInt(x): # allows pickling
    return (x.chr1, x.pos1, x.end1, x.chr2, x.pos2, x.end2)

def loadSVs(fileName, samples=""):
    svs = SortedList([], key=lambda x: (x.chr, x.pos))
    with open(fileName, 'r') as f: #open the file
        lines = f.readlines()
        for line in lines:
            if("#" in line): continue
            values = line.split("\t")
            if(len(values) > 9 and "0/0" in values[9]): continue
            if("ALGORITHM" in values[7]):
                algorithms = values[7].split(";ALGORITHM=")[1].split(";")[0].split(",")
                if(len(algorithms) < 3): continue
            svlen = int(values[7].split("SVLEN=")[1].split(";")[0])
            if(";END=" in values[7]):
                end = int(values[7].split(";END=")[1].split(";")[0])
            else:
                end = int(values[7].split("END=")[1].split(";")[0])
            sv = SV(values[0], int(values[1]), end, svlen, values[4], line, ",".join(values[9:]).split("\n")[0])
            svs.add(sv)
    return svs    

def loadInteractions(fileName, classToMake=Interaction, maxLength=0):
    interactions = SortedList([], key=orderInt)
    with open(fileName, 'r') as f: #open the file
        lines = f.readlines()
        for line in lines:
            values = line.split("\t")
            if(maxLength > 0 and int(values[5])-int(values[1]) > maxLength): continue # to remove if needed
            if(len(values) > 7):
                sample_list = values[7].split("\n")[0].split(",")
            else:
                sample_list = list()
                sample_list.append(fileName.split("/")[-1].split(".")[0])
            interaction = classToMake(values[0], int(values[1]), int(values[2]), values[3], int(values[4]), int(values[5]), int(values[6]), sample_list)
            interactions.add(interaction)
    return interactions

def loadPeaks(fileName):
    macs = ["WTC11"]
    is_mac = False
    if any(mac in fileName for mac in macs):
        is_mac = True
    peaks = SortedList([], key=lambda x: (x.chr, x.pos, x.end))
    with open(fileName, 'r') as f: #open the file
        lines = f.readlines()
        for line in lines:
            values = line.split("\t")
            if not is_mac:
                peak = Peak(values[0], int(values[1]), int(values[2]), values[3], int(values[4]), values[5], float(values[6]), float(values[7]), float(values[8]))
            else:
                peak = Peak(values[0], int(values[1]), int(values[2]), values[3], 0, ".", "0", "-1", "-1")
            peaks.add(peak)
    return peaks

def saveFile(fileName, items):
    with open(fileName, 'w') as f:
        for item in items:
            f.write(item.generateLine())

def checkOverlapPeak(peak1, peak2, minOverlap=1):
    if(peak1.chr != peak2.chr):
        return False
    if(min(peak1.end, peak2.end)-max(peak1.pos, peak2.pos) >= minOverlap):
        return True
    return False

def checkOverlapSV(anchor, sv, minOverlap=1):
    anchorChr, anchorPos, anchorEnd = anchor
    if(anchorChr != sv.chr):
        return False
    if(min(anchorEnd, sv.end)-max(anchorPos, sv.pos) >= minOverlap):
        return True
    return False

def checkOverlap(interaction1, interaction2, minOverlap=1):
    if(interaction1.chr1 != interaction2.chr1 or interaction1.chr2 != interaction2.chr2):
        return False
    if(min(interaction1.end1, interaction2.end1)-max(interaction1.pos1, interaction2.pos1) >= minOverlap):
        if(min(interaction1.end2, interaction2.end2)-max(interaction1.pos2, interaction2.pos2) >= minOverlap):
            return True
    return False

def run_comparison(files):
    file1, file2 = files
    interactions1 = loadInteractions(file1)
    interactions2 = loadInteractions(file2)
    if(len(interactions1) == 0):
        if(len(interactions2) == 0):
            return 100.0
        return 0.0
    interactions_overlap = getOverlapping(interactions1, interactions2)
    return round(len(set(interactions_overlap))/len(set(interactions1))*100, 1)

def getOverlappingSV(interactions, svs): # intersection
    results = list()
    if(len(interactions) == 0):
        return list()
    for interaction in interactions:
        beforeSV = svs.bisect_left(SV(interaction.chr1, interaction.pos1-20000, interaction.pos1-20000, 0, ""))
        afterSV = svs.bisect_right(SV(interaction.chr2, interaction.end2+20000, interaction.end2+20000, 0, ""))

        for sv in svs[beforeSV:min(afterSV, len(svs))]:
            if(sv.chr < interaction.chr1):
                continue
            if(sv.chr > interaction.chr1):
                break
            if checkOverlapSV((interaction.chr1, interaction.pos1, interaction.end1), sv):
                results.append((interaction, sv))
            if checkOverlapSV((interaction.chr2, interaction.pos2, interaction.end2), sv):
                results.append((interaction, sv))

            #interactions.add(Interaction(interaction1.chr1, min(interaction1.pos1, interaction2.pos1), max(interaction1.end1, interaction2.end1),
            #interaction1.chr2, min(interaction1.pos2, interaction2.pos2), max(interaction1.end2, interaction2.end2), interaction1.pet+interaction2.pet))
    return results

def getOverlapping(interactions1, interactions2, toRemove=False): # intersection
    interactions = SortedList([], key=lambda x: (x.chr1, x.pos1, x.end1, x.chr2, x.pos2, x.end2))
    if(toRemove):
        interactions = list()
    if(len(interactions1) == 0):
        return interactions2
    for interaction1 in interactions1:
        beforeInteraction = interactions2.bisect_left(Interaction(interaction1.chr1, interaction1.pos1-20000, interaction1.pos1, interaction1.chr1, interaction1.pos1, interaction1.pos1, interaction1.pet))
        afterInteraction = interactions2.bisect_right(Interaction(interaction1.chr1, interaction1.end1+20000, interaction1.end1, interaction1.chr1, interaction1.end1, interaction1.end1, interaction1.pet))

        for interaction2 in interactions2[beforeInteraction:min(afterInteraction, len(interactions2))]:
            if(interaction2.chr1 < interaction1.chr1):
                continue
            if(interaction2.chr1 > interaction1.chr1):
                break
            if not checkOverlap(interaction1, interaction2):
                continue
            if not(toRemove):
                interactions.add(Interaction(interaction1.chr1, interaction1.pos1, interaction1.end1,
            interaction1.chr2, interaction1.pos2, interaction1.end2, interaction1.pet+interaction2.pet))
            else:
                interactions.append(interaction1)
            #interactions.add(Interaction(interaction1.chr1, min(interaction1.pos1, interaction2.pos1), max(interaction1.end1, interaction2.end1),
            #interaction1.chr2, min(interaction1.pos2, interaction2.pos2), max(interaction1.end2, interaction2.end2), interaction1.pet+interaction2.pet))
    return interactions

def removeOverlapping(interactions1, interactions2): # subtract I1-I2
    interactions = SortedList([], key=lambda x: (x.chr1, x.pos1, x.end1, x.chr2, x.pos2, x.end2))
    for interaction in interactions1:
        interactions.add(interaction)
    to_remove = getOverlapping(interactions, interactions2, True)
    for interaction in to_remove:
        interactions.discard(interaction)
    return interactions

def run_comparison_bed(files, min_overlap=-1):
    file1, file2 = files
    cmd = "cat " + file1 + " | wc -l"
    reference_count = int(subprocess.getoutput(cmd))
    cmd = "cat " + file2 + " | wc -l"
    reference_count2 = int(subprocess.getoutput(cmd))
    reference_count = min(reference_count, reference_count2)
    if(reference_count == 0):
        if(reference_count2 == 0):
            return 100.0
        return 0.0
    if(min_overlap < 0):
        cmd = "bedtools intersect -wa -a "+file1+" -b "+file2+" | cut -f1-3 | uniq | wc -l"
    else:
        cmd = "bedtools intersect -wa -a "+file1+" -b "+file2+" -f "+str(min_overlap)+" | cut -f1-3 | uniq | wc -l"
    common_count = int(subprocess.getoutput(cmd))
    return round(common_count/reference_count*100, 1)

def get_counts(file):
    cmd = "cat " + file + " | uniq | wc -l"
    reference_count = int(subprocess.getoutput(cmd))
    return reference_count

def get_stats(file, ext):
    if(ext == "bed"):
        return (get_counts(file), None, None)
    interactions = loadInteractions(file)
    count_all = len(interactions)
    avg_loop = 0
    avg_anchor = 0
    for interaction in interactions:
        avg_loop += (interaction.end2-interaction.pos1)/count_all
        avg_anchor += ((interaction.end1-interaction.pos1)+(interaction.end2-interaction.pos2))/(2*count_all)
    return (count_all, avg_loop, avg_anchor)
def create_loops(file, folder, peaks=False):
    settings = ""
    settings = '--pet_cutoff 1 --cluster_cutoff 1 --extension 50'
    peaks_line = ""
    if(peaks):
        if(os.path.exists(os.path.splitext(file)[0]+".bed")):
            fileName = os.path.splitext(file)[0]
        else:
            fileName = os.path.splitext(file)[0].split("_R")[0]
        peaks_line = "--peaks_filename " + fileName + ".bed"
    fileName = folder+file.split("/")[-1]
    loop_command = '/home/mateuszchilinski/.pyenv/shims/python /mnt/raid/ctcf_prediction_anal/cluster-paired-end-tags/cluster_pets/cluster_PETs.py '+settings+' --pets_filename '+file+' '+peaks_line+' --clusters_filename '+fileName
    if(peaks):
        loop_command += "_2"
    
    output = subprocess.getoutput(loop_command)
    if(peaks):
        output = subprocess.getoutput("cut -f1-7 " + fileName+"_2 | uniq > " + fileName)
        os.remove(fileName+"_2")

def enlarge_anchors(folder_to_compare, enlargeAnchors):
    files_to_enlarge = [folder_to_compare+f for f in listdir(folder_to_compare) if isfile(join(folder_to_compare, f)) and (f.split(".")[-1] == "bedpe" or f.split(".")[-1] == "BE3")]
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