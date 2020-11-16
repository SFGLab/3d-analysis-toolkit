from common import Peak, loadPeaks, saveFile
from Bio import SeqIO, motifs
import numpy as np
from sortedcontainers import SortedList


peakFile = "/mnt/raid/ctcf_prediction_anal/GM_comparisons_tries/GM12878_R2.bed"
peaks = loadPeaks(peakFile)
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

print("All peaks:" + str(len(peaks)))
peaks_with_motif = SortedList([], key=lambda x: -x.val)
for peak in peaks:
    sequence = genome[peak.chr][peak.pos-extension:peak.end+extension]
    search_results = pssm.search(sequence, threshold=4.0)
    #for position, score in search_results:
    #    print("Position %d: score = %5.3f" % (position, score))
    if(sum(1 for _ in search_results) > 0): # it's generator so its necessary...
        peaks_with_motif.add(peak)
print("Peaks with motif: " + str(len(peaks_with_motif)))

saveFile(peakFile.split(".")[0]+"_2.bed", peaks_with_motif)