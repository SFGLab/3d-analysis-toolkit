from merge_replicates import mergeFilesBedpe
from os import listdir
from os.path import isfile, join
from common import loadInteractions, loadSVs, getOverlappingSV
import shutil
import os

input_folder = "/mnt/raid/ctcf_prediction_anal/3d_analysis_toolkit/stats_svs/"
output_folder = "/mnt/raid/ctcf_prediction_anal/3d_analysis_toolkit/stats_svs/output/"
files_interactions = [input_folder+f for f in listdir(input_folder) if isfile(join(input_folder, f)) and (f.split(".")[-1] == "bedpe" or f.split(".")[-1] == "BE3")]

if os.path.exists(output_folder) and os.path.isdir(output_folder):
    shutil.rmtree(output_folder)
os.mkdir(output_folder)

mergeFilesBedpe(files_interactions, output_folder, "UNION", True)

sv_filename = "prediction/SVs_GTed_only.vcf"
interactions = loadInteractions(output_folder+"UNION.bedpe")

svs = loadSVs(sv_filename)

ovsvs = getOverlappingSV(interactions, svs)


with open(output_folder+'/database.bedpe', 'w') as f:
    for interaction in interactions:
        line = interaction.generateLine(interaction, False)
        for loop, sv in ovsvs:
            if(loop == interaction):
                line += "\t"+sv.getGTs()
        if("\n" not in line):
            line += "\n"
        f.write(line)


print("Done")