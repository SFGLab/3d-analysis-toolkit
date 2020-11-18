# script for generating histograms for anchor & loop (interaction) lengths

import numpy as np
import math
from matplotlib import pyplot as plt

fileName = '/mnt/raid/ctcf_prediction_anal/GM_comparisons_tries/motif_filtered/temp2/enlarged/GM12878_R2_uniq.bedpe'
limit = 4000
limit2 = 2000000
limit3 = 300

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
bins = np.arange(0, limit, 25) # fixed bin size

plt.xlim([0, limit])

plt.hist(data, bins=bins, alpha=0.5)
plt.title('Histogram of anchor lengths in ' + fileName)
plt.xlabel('Length (avg: ' + str(np.round(np.average(data), 2)) + ', std:' + str(np.round(np.std(data), 2)) + ', shown: ' + str(np.round(sum(x <= limit for x in data)/len(data)*100, 2)) + '%)')
plt.ylabel('Count')

plt.show()

bins = np.arange(0, limit2, 5000) # fixed bin size

plt.xlim([0, limit2])

plt.hist(data2, bins=bins, alpha=0.5)
plt.title('Histogram of loops lengths in ' + fileName)
plt.xlabel('Length (avg: ' + str(np.round(np.average(data2), 2)) + ', std:' + str(np.round(np.std(data2), 2)) + ', shown: ' + str(np.round(sum(x <= limit2 for x in data2)/len(data2)*100, 2)) + '%)')
plt.ylabel('Count')

plt.show()

bins = np.arange(0, limit3, 1) # fixed bin size

plt.xlim([0, limit3])

plt.hist(data3, bins=bins, alpha=0.5)
plt.title('Histogram of pet counts in ' + fileName)
plt.xlabel('PET (avg: ' + str(np.round(np.average(data3), 2)) + ', std:' + str(np.round(np.std(data3), 2)) + ', shown: ' + str(np.round(sum(x <= limit3 for x in data3)/len(data3)*100, 2)) + '%)')
plt.ylabel('Count')

plt.show()