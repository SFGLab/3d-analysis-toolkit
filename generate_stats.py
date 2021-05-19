# script for generating histograms for anchor & loop (interaction) lengths

import numpy as np
import math
from matplotlib import pyplot as plt
import os

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

fileNames = ['/mnt/raid/legolas/ctcf_prediction_anal/cremins_data/HFF_main_dot_calls.bedpe',
'/mnt/raid/legolas/ctcf_prediction_anal/cremins_data/HFFC6.bedpe'] # first file determines if it's bed or bedpe!
legendNames = [fileName.split("/")[-1].split(".")[0] for fileName in fileNames]

limit = 20000
limit2 = 2000000
limit3 = 100
xLogLoops = True
yLog = True
saveFig = True

datas1 = []
datas2 = []
datas3 = []
main_ext = fileNames[0].split("/")[-1].split(".")[1]

if(main_ext == "bedpe"):
    for fileName in fileNames:
        data = list() # anchors
        data2 = list() # loops
        data3 = list() # pet count
        with open(fileName) as f:
            for line in f:
                values = line.split("\t")
                data.append(int(values[2])-int(values[1]))
                data.append(int(values[5])-int(values[4]))
                data2.append(int(values[5])-int(values[1]))
                data3.append(int(values[6]))
        datas1.append(data)
        datas2.append(data2)
        datas3.append(data3)
else:
    for fileName in fileNames:
        data = list() # ccds
        with open(fileName) as f:
            for line in f:
                values = line.split("\t")
                data.append(int(values[2])-int(values[1]))
        datas2.append(data)
if(main_ext == "bedpe"):
    bins = np.arange(0, limit, 25) # fixed bin size

    plt.xlim([0, limit])

    if(yLog):
        plt.yscale('log', nonposy='clip')

    plt.hist(datas1, bins=bins, alpha=0.5)
    plt.title('Histogram of anchor lengths')

    if(len(datas1) == 1):
        plt.xlabel('Length (avg: ' + str(np.round(np.average(data), 2)) + ', std:' + str(np.round(np.std(data), 2)) + ', shown: ' + str(np.round(sum(x <= limit for x in data)/len(data)*100, 2)) + '%)')
    plt.legend(legendNames)
    plt.ylabel('Count')

    plt.show()

bins = np.arange(0, limit2, 5000) # fixed bin size

hist, bins, _ = plt.hist(datas2, bins=bins, alpha=0.5)

logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))

if(xLogLoops):
    plt.hist(datas2, bins=logbins)
    plt.xscale('log')
plt.xlim(left=10000)

if(yLog):
    plt.yscale('log', nonposy='clip')

if(main_ext == "bedpe"):
    plt.title('Histogram of loops lengths')
else:
    plt.title('Histogram of CCDs/TADs lengths')
if(len(datas2) == 1):
    plt.xlabel('Length (avg: ' + str(np.round(np.average(data2), 2)) + ', std:' + str(np.round(np.std(data2), 2)) + ', shown: ' + str(np.round(sum(x <= limit2 for x in data2)/len(data2)*100, 2)) + '%)')
plt.legend(legendNames)
plt.ylabel('Count')

plt.show()

if(main_ext == "bedpe"):
    bins = np.arange(0, limit3, 1) # fixed bin size

    plt.xlim([0, limit3])

    if(yLog):
        plt.yscale('log', nonposy='clip')

    plt.hist(datas3, bins=bins, alpha=0.5)
    plt.title('Histogram of pet counts')
    if(len(datas3) == 1):
        plt.xlabel('PET (avg: ' + str(np.round(np.average(data3), 2)) + ', std:' + str(np.round(np.std(data3), 2)) + ', shown: ' + str(np.round(sum(x <= limit3 for x in data3)/len(data3)*100, 2)) + '%)')
    plt.legend(legendNames)
    plt.ylabel('Count')

    plt.show()