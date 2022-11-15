import subprocess
import pandas as pd
from io import StringIO

file_mother = "/mnt/raid/repos/hichip/YRI/GM19238.bedpe"
file_father = "/mnt/raid/repos/hichip/YRI/GM19239.bedpe"
file_child = "/mnt/raid/repos/hichip/YRI/GM19240.bedpe"

# file_mother = "/mnt/raid/repos/hichip/PUR/HG00732.bedpe"
# file_father = "/mnt/raid/repos/hichip/PUR/HG00731.bedpe"
# file_child = "/mnt/raid/repos/hichip/PUR/HG00733.bedpe"

# file_mother = "/mnt/raid/repos/hichip/CHS/HG00513.bedpe"
# file_father = "/mnt/raid/repos/hichip/CHS/HG00512.bedpe"
# file_child = "/mnt/raid/repos/hichip/CHS/HG00514.bedpe"

child = pd.read_csv(file_child, sep="\t", names=["chr1", "pos1", "end1", "chr2", "pos2", "end2"])
mother = pd.read_csv(file_mother, sep="\t", names=["chr1", "pos1", "end1", "chr2", "pos2", "end2"])
father = pd.read_csv(file_father, sep="\t", names=["chr1", "pos1", "end1", "chr2", "pos2", "end2"])

child_mother = pd.merge(child, mother, how="inner")
child_father = pd.merge(child, father, how="inner")

child_mother_father = pd.merge(child_mother, child_father, how="inner")

child_mother_uniq = pd.concat([child_mother,child_mother_father]).drop_duplicates(keep=False)
child_father_uniq = pd.concat([child_father,child_mother_father]).drop_duplicates(keep=False)

child_uniq = pd.concat([child,child_mother,child_father]).drop_duplicates(keep=False)

child_uniq["genotype"] = ".|."
child_mother_uniq["genotype"] = "1|0"
child_father_uniq["genotype"] = "0|1"
child_mother_father["genotype"] = "1|1"

child_phased = pd.concat([child_uniq,child_mother_uniq,child_father_uniq,child_mother_father]).sort_values(["chr1", "pos1", "end1", "chr2", "pos2", "end2"]).reset_index(drop=True)

print("Phased child sample: " + file_child)
print("Based on: " + file_mother + " (left)")
print("Based on: " + file_father + " (right)")
print("")
print("Results:")
all_count = child_phased.count()["chr1"]
unable_to_phase = child_phased[child_phased["genotype"] == ".|."].count()["chr1"]
print(".|. (unable to phase): " + str(unable_to_phase) + " (" + str(round(unable_to_phase/all_count*100.0, 2)) + "%)")
mother_only = child_phased[child_phased["genotype"] == "1|0"].count()["chr1"]
print("1|0 (only mother): " + str() + " (" + str(round(mother_only/all_count*100.0, 2)) + "%)")
father_only = child_phased[child_phased["genotype"] == "0|1"].count()["chr1"]
print("0|1 (only father): " + str() + " (" + str(round(father_only/all_count*100.0, 2)) + "%)")
both = child_phased[child_phased["genotype"] == "1|1"].count()["chr1"]
print("1|1 (both): " + str() + " (" + str(round(both/all_count*100.0, 2)) + "%)")
