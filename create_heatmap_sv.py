import subprocess
import pandas as pd

import seaborn as sns
from matplotlib import pyplot as plt
import subprocess

def run_command(command):
    print("___COMMAND: " + command)
    process = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    print(process.stderr)
    return process.stdout

command = "bedtools intersect -wa -f 0.8 -r -a %s -b %s | uniq | wc -l"
sv_files = ["vcfs/consensuSV__HG00512.vcf", "vcfs/consensuSV__HG00513.vcf", "vcfs/consensuSV__HG00514.vcf",
            "vcfs/consensuSV__HG00731.vcf", "vcfs/consensuSV__HG00732.vcf", "vcfs/consensuSV__HG00733.vcf",
            "vcfs/consensuSV__NA19238.vcf", "vcfs/consensuSV__NA19239.vcf", "vcfs/consensuSV__NA19240.vcf"]

sv_files_results = {}


for sv_file in sv_files:
    sv_files_results[sv_file] = {}
    sv_files_results[sv_file]["total"] = int(str(run_command("grep -vc \"#\" %s" % (sv_file))).split("b'")[1].split("\\n")[0])
    sv_files_results[sv_file]["sample"] = str(run_command("bcftools query -l %s" % (sv_file))).split("b'")[1].split("\\n")[0]
    


for sv_file in sv_files:
    for sv_file2 in sv_files:
        sv_files_results[sv_file][sv_file2] = int(str(run_command(command % (sv_file, sv_file2))).split("b'")[1].split("\\n")[0])

sv_full_results = {}

for sv_file in sv_files:
    sv_full_results[sv_files_results[sv_file]["sample"]] = {}
    for sv_file2 in sv_files:
        sv_full_results[sv_files_results[sv_file]["sample"]][sv_files_results[sv_file2]["sample"]] = sv_files_results[sv_file][sv_file2]/sv_files_results[sv_file]["total"]
sv_full_results = pd.DataFrame(sv_full_results)
sns.set(rc={'figure.figsize':(8,8)})
sns.heatmap(sv_full_results, cmap="Reds", annot=True)
plt.savefig("heatmap.png", dpi=800)