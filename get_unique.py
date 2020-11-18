from common import saveFile, removeOverlapping, loadInteractions

file1 = "/mnt/raid/ctcf_prediction_anal/GM_comparisons_tries/motif_filtered/temp2/enlarged/GM12878_R2.bedpe"
file2 = "/mnt/raid/ctcf_prediction_anal/GM_comparisons_tries/motif_filtered/temp2/enlarged/GM12878_R1.bedpe"

file_out = "/mnt/raid/ctcf_prediction_anal/GM_comparisons_tries/motif_filtered/temp2/enlarged/"+file1.split("/")[-1].split(".")[0]+"_uniq.bedpe"

inter1 = loadInteractions(file1)
uniq_inter = removeOverlapping(inter1, loadInteractions(file2))
print("Found " + str(len(uniq_inter)) + " unique interactions, accounting for " + str(round(len(uniq_inter)/len(inter1)*100, 1))+"%"+" of the file. Saved.")
saveFile(file_out, uniq_inter)