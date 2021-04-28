
import os
import pyBigWig
import math
import pandas as pd
import time
start_time = time.time()

window = 1000
n = window/25
folder = "/mnt/raid/ctcf_prediction_anal/3d_analysis_toolkit/encode_data/"


result = pd.DataFrame()

files_i = 1
all_files = len(os.listdir(folder))-1
for filename in os.listdir(folder):
    start_time_file = time.time()
    if not(filename.endswith(".bigWig")): continue 
    data = []
    full_filename = os.path.join(folder, filename)
    filename_short = filename.split(".")[0]
    bw = pyBigWig.open(full_filename)
    for chrm, total in bw.chroms().items(): # can be easily paralleled
        for i in range(0, math.floor(total/window)):
            data.append(
                {"chr": chrm, 
                "pos": window*i, 
                "end": window*(i+1), 
                filename_short: sum(bw.values(chrm, window*i, window*(i+1)))/n
                })
    if(len(result) == 0):
        result = pd.DataFrame(data)
    else:
        df = pd.DataFrame(data)
        result = pd.merge(df, result, how='outer', on=["chr", "pos", "end"])
    print("Processed file " + str(files_i) + " out of " + str(all_files) + " | Processing time: %s seconds" % (time.time() - start_time))
    files_i += 1

result.to_csv("xd.csv", index=False)
print("--- %s seconds ---" % (time.time() - start_time))