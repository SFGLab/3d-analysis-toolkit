from common import Peak, loadInteractions, saveFile
from Bio import SeqIO, motifs
import numpy as np
from sortedcontainers import SortedList

interactionsFile = "/mnt/raid/ctcf_prediction_anal/GM_comparisons_tries/GM12878_R1.bedpe"
interactions = loadInteractions(interactionsFile)
genome = dict()

pfmMatrixFile = "/mnt/raid/ctcf_prediction_anal/3d_analysis_toolkit/MA0139.1.pfm"

for record in SeqIO.parse("/mnt/raid/trios_data/GRCh38_full_analysis_set_plus_decoy_hla.fa", "fasta"):
    if("_" in record.id or "*" in record.id or "EBV" in record.id):
        continue
    genome[record.id] = record.seq

extension = 0
with open(pfmMatrixFile) as handle:
    ctcf_motif = motifs.read(handle, "pfm")

pssm = ctcf_motif.pssm

print("All interactions:" + str(len(interactions)))
interactions_with_motif = SortedList([], key=lambda x: (x.chr1, x.pos1, x.end1, x.chr2, x.pos2, x.end2))
for interaction in interactions:
    sequence1 = genome[interaction.chr1][interaction.pos1-extension:interaction.end1+extension]
    sequence2 = genome[interaction.chr1][interaction.pos1-extension:interaction.end1+extension]

    search_results1 = [(x,y) for x,y in pssm.search(sequence1, threshold=7.0)]
    search_results2 = [(x,y) for x,y in pssm.search(sequence2, threshold=7.0)]

    full_results = list()
    added = False
    for r1 in search_results1:
        for r2 in search_results2:
            if(r1[0]*r2[0] >= 0):
                interactions_with_motif.add(interaction)
                added = True
                break
        if(added):
            break

print("Interactions with motif: " + str(len(interactions_with_motif)))

saveFile(interactionsFile.split(".")[0]+"_3.bedpe", interactions_with_motif)