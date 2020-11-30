import argparse
import progressbar
import shutil
import shlex
import os
from common import Interaction, Peak, createFolder, removeFolder, loadInteractions, loadPeaks, checkOverlap, checkOverlapPeak
from os import listdir
from os.path import isfile, join
import math 
def mergeFilesBedpe(interactions_files, output_folder, sample_name):
    interactionsReplicates = list()

    for file in interactions_files:
        interactions = loadInteractions(file)
        interactionsReplicates.append(interactions)

    allInteractions = interactionsReplicates[0]
    if(len(interactions_files) < 2):
        with open(output_folder+'/'+sample_name+'.bedpe', 'w') as f:
            for interaction in allInteractions:
                f.write(interaction.generateLine())
        return
    fileId = 2
    for interactionsReplicate in interactionsReplicates[1:]:
        print("Merging interactions from file #" + str(fileId))
        bar = progressbar.ProgressBar(maxval=len(interactionsReplicate), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
        bar.start()
        i = 0
        for interaction in interactionsReplicate:
            bar.update(i + 1)
            i = i + 1
            beforeInteraction = allInteractions.bisect_left(Interaction(interaction.chr1, interaction.pos1-5000, interaction.pos1, interaction.chr1, interaction.pos1, interaction.pos1, interaction.pet))
            if(beforeInteraction > len(allInteractions)-1): # is after all, just to be sure
                allInteractions.add(interaction)
                continue
            afterInteraction = allInteractions.bisect_right(Interaction(interaction.chr1, interaction.end1+5000, interaction.end1, interaction.chr1, interaction.end1, interaction.end1, interaction.pet))
            
            found_interaction = False
            for interactionId in range(beforeInteraction, min(len(allInteractions)-1, afterInteraction+1)):
                if(checkOverlap(allInteractions[interactionId], interaction)): # check if overlaps on both sides
                    allInteractions[interactionId].pet += interaction.pet
                    # merge anchors
                    allInteractions[interactionId].pos1 = min(allInteractions[interactionId].pos1, interaction.pos1)
                    allInteractions[interactionId].end1 = max(allInteractions[interactionId].end1, interaction.end1)

                    allInteractions[interactionId].pos2 = min(allInteractions[interactionId].pos2, interaction.pos2)
                    allInteractions[interactionId].end2 = max(allInteractions[interactionId].end2, interaction.end2)

                    found_interaction = True
                    break
            if not found_interaction: # it does not overlap with any, lets add it as independent one
                allInteractions.add(interaction)
        bar.finish()
    with open(output_folder+'/'+sample_name+'.bedpe', 'w') as f:
        for interaction in allInteractions:
            f.write(interaction.generateLine())

    print("Finished parsing!")

def mergeFilesBed(peaks_files, output_folder, sample_name):
    peaksReplicates = list()

    for file in peaks_files:
        peaks = loadPeaks(file)
        peaksReplicates.append(peaks)

    allPeaks = peaksReplicates[0]
    if(len(peaks_files) < 2):
        with open(output_folder+'/'+sample_name+'.bedpe', 'w') as f:
            for peak in allPeaks:
                f.write(peak.generateLine())
        return
    fileId = 2
    for peaksReplicate in peaksReplicates[1:]:
        print("Merging peaks from file #" + str(fileId))
        bar = progressbar.ProgressBar(maxval=len(peaksReplicate), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
        bar.start()
        i = 0
        for peak in peaksReplicate:
            bar.update(i + 1)
            i = i + 1
            beforePeak = allPeaks.bisect_left(Peak(peak.chr, peak.pos-5000, peak.end, peak.name, peak.score, peak.strand, peak.signalValue, peak.pValue, peak.qValue))
            if(beforePeak > len(allPeaks)-1): # is after all, just to be sure
                allPeaks.add(peak)
                continue
            afterPeak = allPeaks.bisect_right(Peak(peak.chr, peak.pos+5000, peak.end, peak.name, peak.score, peak.strand, peak.signalValue, peak.pValue, peak.qValue))
            
            found_peak = False
            for peakId in range(beforePeak, min(len(allPeaks)-1, afterPeak+1)):
                if(checkOverlapPeak(allPeaks[peakId], peak)): # check if overlaps on both sides
                    allPeaks[peakId].pValue = max(allPeaks[peakId].pValue, peak.pValue)
                    # merge anchors
                    allPeaks[peakId].pos = min(allPeaks[peakId].pos, peak.pos)
                    allPeaks[peakId].end = max(allPeaks[peakId].end, peak.end)

                    found_peak = True
                    break
            if not found_peak: # it does not overlap with any, lets add it as independent one
                allPeaks.add(peak)
        bar.finish()
        fileId += 1
    with open(output_folder+'/'+sample_name+'.bed', 'w') as f:
        for peak in allPeaks:
            f.write(peak.generateLine())

    print("Finished parsing!")

#parser = argparse.ArgumentParser()
#parser.add_argument('interactions_files', metavar='file', help="Files containing ", nargs='+')

#args = parser.parse_args()
input_folder = "/mnt/raid/ctcf_prediction_anal/trios_new_ctcf/ctcf_named/"
output_folder = "/mnt/raid/ctcf_prediction_anal/trios_new_ctcf/ctcf_named/output/"

if os.path.exists(output_folder) and os.path.isdir(output_folder):
    shutil.rmtree(output_folder)
os.mkdir(output_folder)

samples = {}

files_interactions = [input_folder+f for f in listdir(input_folder) if isfile(join(input_folder, f)) and (f.split(".")[-1] == "bedpe" or f.split(".")[-1] == "BE3")]

for file in files_interactions:
    sample_name = file.split("/")[-1]
    if("_" in sample_name):
        sample_name = sample_name.split("_")[0]
    else:
        sample_name = sample_name.split(".")[0]
    if(sample_name not in samples):
        samples[sample_name] = list()
    samples[sample_name].append(file)

files_peaks = [input_folder+f for f in listdir(input_folder) if isfile(join(input_folder, f)) and (f.split(".")[-1] == "bed")]

samples_bed = {}

for file in files_peaks:
    sample_name = file.split("/")[-1]
    if("_" in sample_name):
        sample_name = sample_name.split("_")[0]
    else:
        sample_name = sample_name.split(".")[0]
    if(sample_name not in samples_bed):
        samples_bed[sample_name] = list()
    samples_bed[sample_name].append(file)

#interactions_files = args.interactions_files

for sample_name, interactions_files in samples.items():
    mergeFilesBedpe(interactions_files, output_folder, sample_name)

for sample_name, peaks_files in samples_bed.items():
    mergeFilesBed(peaks_files, output_folder, sample_name)


print("done")
