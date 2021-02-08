from merge_replicates import mergeFilesBedpe
from os import listdir
from os.path import isfile, join
from common import loadInteractions, loadSVs, getOverlappingSV
import shutil
import os
import sqlite3


samples = ["HG00512", "HG00513", "HG00514", "HG00731", "HG00732", "HG00733", "NA19238", "NA19239", "NA19240"]

samples_to_table = ""

for sample in samples:
    samples_to_table += '''
	"'''+sample+'''_anchor"	INTEGER,'''
for sample in samples:
    samples_to_table += '''
	"'''+sample+'''_sv"	TEXT'''
    if(sample != samples[-1]):
        samples_to_table += ''',
        '''
conn = sqlite3.connect('TestDB.db')  # You can create a new database by changing the name within the quotes
c = conn.cursor() # The database will be saved in the location where your 'py' file is saved

# Create table - CLIENTS
c.execute('''CREATE TABLE "anchor_database" (
	"Id"	INTEGER PRIMARY KEY AUTOINCREMENT,
	"Chr1"	TEXT,
	"Chr2"	TEXT,
	"Pos1"	INTEGER,
	"End1"	INTEGER,
	"Pos2"	INTEGER,
	"End2"	INTEGER,
	"Freq"	INTEGER,
    "SV_type" TEXT,
    '''+samples_to_table+'''
);''')

conn.commit()

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


for interaction in interactions:
    sql_line_base = "INSERT INTO \"anchor_database\" VALUES (NULL, '"+interaction.chr1+"', '"+interaction.chr2+"', "+str(interaction.pos1)+", "+str(interaction.end1)+", "+str(interaction.pos2)+", "+str(interaction.end2)+", "+str(interaction.pet)
    
    sql_line_later = ""
    for sample in samples:
        if(sample in interaction.samples):
            sql_line_later += ", 1"
        else:
            sql_line_later += ", 0"

    found = False
    for loop, sv in ovsvs:
        if(loop == interaction):
            sql_line = sql_line_base
            sql_line += ", '"+sv.svtype+"'"
            sql_line += sql_line_later
            for gt in sv.getGTs():
                sql_line += ", '"+gt+"'"
            sql_line += ");"
            c.execute(sql_line)
            found = True
    if not(found):
        sql_line_base += ", 'NONE'"+sql_line_later+"".join([", '0/0'" for sample in samples]) + ");"
        c.execute(sql_line_base)

conn.commit()

print("Done")