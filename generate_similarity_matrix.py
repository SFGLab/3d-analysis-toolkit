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
from common import Interaction, createFolder, removeFolder, loadInteractions, run_comparison, run_comparison_bed, get_counts, create_loops, enlarge_anchors, saveFile, removeOverlapping
from call_ccds import get_ccds
from find_motifs import filterInteractionsByMotifs
from pandas_profiling import ProfileReport
from jinja2 import Environment, FileSystemLoader
from weasyprint import HTML
import seaborn as sns
import time

def applyColouring(df, low=0, high=100, percentage=True):
    df = df.apply(pd.to_numeric).style.background_gradient("RdYlGn", axis=None, subset=df.columns).set_table_styles([
            {'selector': 'th', 'props': [('border-style','solid'), ('border-color','black'),('border-width','1px')]},
            {'selector': 'td', 'props': [('border-style','solid'), ('border-color','black'),('border-width','1px')]},
            {'selector': 'table', 'props': [('border-style','solid'), ('border-color','black'),('border-width','1px')]},
            {'selector': 'tr', 'props': [('border-style','solid'), ('border-color','black'),('border-width','1px')]}
            ])
    if(percentage):
        df = df.format('{:.2f}%'.format)
    return df.render()

def generateReportSection(toReport):

    (count, regular, maximum, avg) = toReport
    
    code = "<h2>Counts</h2>\n"
    code += applyColouring(count, count.min().min(), count.max().max(), False)
    with(pd.option_context("display.float_format", '{:.2f}%'.format)):
        code += "<h2>Similarity matrix</h2>\n"
        code += applyColouring(regular)

        code += "<h2>Similarity matrix - maximum</h2>\n"
        code += applyColouring(maximum)

        code += "<h2>Similarity matrix - average</h2>\n"
        code += applyColouring(avg)
    return code

def generateHTMLReport(options, peaks, interactions, loops_no_peaks, loops_peaks, ccds):
    env = Environment(loader=FileSystemLoader('.'))
    template = env.get_template("template.html")
    (filterMotifs, maxLength, enlargeAnchors, randomSampling) = options
    content = "<h1>Settings</h1>\n"
    content += "<table border=1><tr><th>Setting</th><th>Value</th></tr>"
    content += "<tr><td>Filter Motifs</td><td>"+str(filterMotifs)+"</td></tr>"
    content += "<tr><td>Max Length</td><td>"+str(maxLength)+"</td></tr>"
    content += "<tr><td>Enlarge Anchors</td><td>"+str(enlargeAnchors)+"</td></tr>"
    content += "<tr><td>Anchors Apart By</td><td>"+str(2*enlargeAnchors)+"</td></tr>"
    content += "<tr><td>Random Sampling</td><td>"+str(randomSampling)+"</td></tr></table>"

    content += "<h1>Peaks</h1>\n"
    content += generateReportSection(peaks)

    content += "<h1>Interactions</h1>\n"
    content += generateReportSection(interactions)

    content += "<h1>Loops (no peaks)</h1>\n"
    content += generateReportSection(loops_no_peaks)

    content += "<h1>Loops (peaks)</h1>\n"
    content += generateReportSection(loops_peaks)

    content += "<h1>CCDs (based on loops with peaks, enlarging always false)</h1>\n"
    content += generateReportSection(ccds)

    template_vars = {"content": content}

    html_out = template.render(template_vars)
    ts = str(int(time.time()))
    with open('report'+ts+'.html', 'w') as f:
        f.write(html_out)
    #HTML(string=html_out).write_pdf("report.pdf")

def generate_matrix(folder_to_compare, enlargeAnchors=0, func_to_use=run_comparison, ext="bedpe", getSimilarityMatrices=True, generateReport=False):
    if(getSimilarityMatrices):
        if enlargeAnchors > 0:
            createFolder(folder_to_compare+"enlarged/")
            enlarge_anchors(folder_to_compare, enlargeAnchors)
            folder_to_compare += "enlarged/"
        threads = 16
        if(ext=="bedpe"):
            files_to_compare = [folder_to_compare+f for f in listdir(folder_to_compare) if isfile(join(folder_to_compare, f)) and (f.split(".")[-1] == "bedpe" or f.split(".")[-1] == "BE3")]
        else:
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

        df = df.sort_index().sort_index(axis = 1)
        with pd.option_context('display.max_rows', None, 'display.max_columns', None, 'display.width', None, "display.float_format", '{:.2f}%'.format):  # more options can be specified also
            df_full = pd.concat(
                [pd.concat(
                    [df],
                    keys=['File'], axis=1)],
                keys=['Reference']
            )
            df_max = df.where(df > df.transpose(), df.transpose()).fillna(df)
            df_avg = ((df+df.transpose())/2.0).fillna(df)
            
            if(generateReport == False):
                print(df_full)
                print(df_max)
                print(df_avg)
    else:
        if(ext=="bedpe"):
            files_to_compare = [folder_to_compare+f for f in listdir(folder_to_compare) if isfile(join(folder_to_compare, f)) and (f.split(".")[-1] == "bedpe" or f.split(".")[-1] == "BE3")]
        else:
            files_to_compare = [folder_to_compare+f for f in listdir(folder_to_compare) if isfile(join(folder_to_compare, f)) and f.split(".")[-1] == ext]
    print("")
    count_of_what = "interactions"
    if(ext=="bed"):
        count_of_what = "peaks"

    counts_df = pd.DataFrame(columns=['Count'])
    for file in files_to_compare:
        counts_df = counts_df.append(pd.Series({'Count': get_counts(file)}, name=file.split("/")[-1].split(".")[0]))
    counts_df = pd.concat(
                [pd.concat(
                    [counts_df],
                    keys=['File'], axis=0)]
            ).transpose()
    if(generateReport == False):
        print("Counts of "+count_of_what+" in files:")
        print(counts_df)
    else:
        return (counts_df, df_full, df_max, df_avg)

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

def filterMotifsFunc(samples):
    start_time = time.time()
    print("===== FILTERING MOTIFS IS ON =====")
    
    samples_2 = dict()
    
    task_list = [(sample, interactions) for sample, interactions in samples.items()]
    threads = 16
    pool = mp.Pool(threads)
    results = pool.map(filterInteractionsByMotifsParallel, task_list)
    pool.close()
    pool.join() 

    for result in results:
        sample, interactions = result
        samples_2[sample] = interactions
    print("--- Executed in %s seconds ---" % (time.time() - start_time))
    return samples_2

def run_comparison_bed_ccd(files):
    return run_comparison_bed(files, 0.5)

start_time_total = time.time()

#folder_to_compare = '/mnt/raid/ctcf_prediction_anal/GM_comparisons_tries/'
randomSampling = False
filterMotifs = True
getSimilarityMatrices = True
generateReport = True
enlargeAnchors = 1000 # 0 = disabled
maxLength = 500000
folder_to_compare = '/mnt/raid/ctcf_prediction_anal/trios_new_ctcf/ctcf_named/output/output2/'
#folder_to_compare = '/mnt/raid/ctcf_prediction_anal/trios_new_ctcf/ctcf_named/'
rs_temp = ""

print("===== PEAKS =====")

peaks_matrix = generate_matrix(folder_to_compare,0,run_comparison_bed, "bed", getSimilarityMatrices, generateReport=True)

if(generateReport):
    print("Generated, added to report.")

if(getSimilarityMatrices):
    if randomSampling or filterMotifs or maxLength > 0:
        rs_temp = "modified_temp/"
        createFolder(folder_to_compare+rs_temp)
        files_to_copy = [folder_to_compare+f for f in listdir(folder_to_compare) if isfile(join(folder_to_compare, f)) and f.split(".")[-1] == "bed"]
        for file in files_to_copy:
            shutil.copy2(file, folder_to_compare+rs_temp+file.split("/")[-1])

        files_to_compare = [folder_to_compare+f for f in listdir(folder_to_compare) if isfile(join(folder_to_compare, f)) and (f.split(".")[-1] == "bedpe" or f.split(".")[-1] == "BE3")]

        samples = dict()

        for file in files_to_compare:
            sample = file.split("/")[-1]
            if("BE3" in sample):
                sample = sample.replace("BE3", "bedpe")
            samples[folder_to_compare+rs_temp+sample] = loadInteractions(file, maxLength=maxLength)

    if randomSampling:
        print("===== SAMPLING IS ON =====")
        samples_2 = dict()
        for sample, interactions in samples.items():
            samples_2[sample] = createRandomSample(interactions, 200000)
        samples = samples_2

    if filterMotifs:
        samples = filterMotifsFunc(samples)


    if randomSampling or filterMotifs or maxLength > 0:
        for sample, interactions in samples.items():
            saveFile(sample, interactions)

print("===== INTERACTIONS =====")
interactions_matrix = generate_matrix(folder_to_compare+rs_temp, enlargeAnchors, getSimilarityMatrices=getSimilarityMatrices, generateReport=True)

if(generateReport):
    print("Generated, added to report.")

createFolder(folder_to_compare+rs_temp+"temp")
createFolder(folder_to_compare+rs_temp+"temp2")
createFolder(folder_to_compare+rs_temp+"temp3")

files_to_compare = [folder_to_compare+rs_temp+f for f in listdir(folder_to_compare+rs_temp) if isfile(join(folder_to_compare+rs_temp, f)) and (f.split(".")[-1] == "bedpe" or f.split(".")[-1] == "BE3")]

for file in files_to_compare:
    create_loops(file, folder_to_compare+rs_temp+"temp/", False)
    if(randomSampling):
        createRandomSampleFile(file, folder_to_compare+rs_temp+"temp/", 20000)
print("===== LOOPS (NO PEAKS) =====")
loops_no_peaks_matrix = generate_matrix(folder_to_compare+rs_temp+"temp/", enlargeAnchors, getSimilarityMatrices=getSimilarityMatrices, generateReport=True)


if(generateReport):
    print("Generated, added to report.")

for file in files_to_compare:
    if(os.path.isfile(os.path.splitext(file)[0]+".bed") or os.path.isfile(os.path.splitext(file)[0].split("_R")[0]+".bed")): 
        create_loops(file, folder_to_compare+rs_temp+"temp2/", True)
        if(randomSampling):
            createRandomSampleFile(file, folder_to_compare+rs_temp+"temp2/", 10000)
print("===== LOOPS (PEAKS) =====")
loops_peaks_matrix = generate_matrix(folder_to_compare+rs_temp+"temp2/", enlargeAnchors, getSimilarityMatrices=getSimilarityMatrices, generateReport=True)

loops_location = folder_to_compare+rs_temp+"temp2/"

files_to_compare = [loops_location+f for f in listdir(loops_location) if isfile(join(loops_location, f)) and (f.split(".")[-1] == "bedpe" or f.split(".")[-1] == "BE3")]

for file in files_to_compare:
    fileName = os.path.splitext(file.split("/")[-1])[0]
    get_ccds(file, folder_to_compare+rs_temp+"temp3/"+fileName+".bed", 3, 2, 10000)

print("===== CCDs (based on loops with peaks, enlarging always false) =====")
ccds_peaks_matrix = generate_matrix(folder_to_compare+rs_temp+"temp3/", 0, run_comparison_bed_ccd, "bed", getSimilarityMatrices=getSimilarityMatrices, generateReport=True)

if(generateReport):
    print("Generated, added to report.")
    generateHTMLReport((filterMotifs, maxLength, enlargeAnchors, randomSampling), peaks_matrix, interactions_matrix, loops_no_peaks_matrix, loops_peaks_matrix, ccds_peaks_matrix)
    print("Report generated, saved.")

print("--- Executed in %s seconds ---" % (time.time() - start_time_total))



#removeFolder(folder_to_compare+rs_temp+"temp")
#removeFolder(folder_to_compare+rs_temp+"temp2")

if randomSampling:
    removeFolder(folder_to_compare+rs_temp)