import pandas as pd
import numpy as np
import time
from Bio import SeqIO, motifs
from common import loadInteractions, loadSVs, getOverlappingSV, checkOverlapSV
from find_motifs import decideOne

def parse_one(genome, pssm, anchor1, anchor2, which_anchor, sv):
    extension = 200

    if(which_anchor == 1):
        chrm, pos, end = anchor1
    else:
        chrm, pos, end = anchor2
        
    if(sv.pos >= pos and sv.end <= end): # sv all inside
        sequence = genome[chrm][pos-extension:sv.pos+extension] + genome[chrm][sv.end-extension:end+extension]
    elif(sv.pos <= pos and sv.end >= end): # sv all outside - no anchor at all
        sequence = ""
    elif(sv.pos < pos): # sv starts before anchor, ends in anchor
        sequence = genome[chrm][sv.end-extension:end+extension]
    elif(sv.pos >= pos): # sv starts in anchor, ends after anchor 
        sequence = genome[chrm][pos-extension:sv.pos+extension]
    # left anchor    with np.errstate(invalid='ignore'):
    if(which_anchor == 1):
        seq1 = sequence
        seq2 = genome[anchor2[0]][anchor2[1]-extension:anchor2[2]+extension]
    else:
        seq1 = genome[anchor1[0]][anchor1[1]-extension:anchor1[2]+extension]
        seq2 = sequence
    if(decideOne(pssm, seq1, seq2)):
        return True
    else:
        return False

def predict_loops(sv_filename, loops_filename):
    pfmMatrixFile = "/mnt/raid/ctcf_prediction_anal/3d_analysis_toolkit/MA0139.1.pfm"
    with open(pfmMatrixFile) as handle:
        ctcf_motif = motifs.read(handle, "pfm")
    pssm = ctcf_motif.pssm

    genome = dict()
    for record in SeqIO.parse("/mnt/raid/trios_data/GRCh38_full_analysis_set_plus_decoy_hla.fa", "fasta"):
        if("_" in record.id or "*" in record.id or "EBV" in record.id):
            continue
        genome[record.id] = record.seq
    
    interactions = loadInteractions(loops_filename)
    svs = loadSVs(sv_filename)
    overlapping = getOverlappingSV(interactions, svs)
    # candidates
    loops_altered = 0
    total_dup_intersect = 0
    sv_count = set()
    results = list()
    for (interaction, sv) in overlapping:
        if(sv.svtype == "<DEL>"):
            if checkOverlapSV((interaction.chr1, interaction.pos1, interaction.end1), sv):
                if(parse_one(genome, pssm, (interaction.chr1, interaction.pos1, interaction.end1), (interaction.chr2, interaction.pos2, interaction.end2), 1, sv)):
                    pass
                else:
                    loops_altered += 1
                    sv_count.add(sv)
            if checkOverlapSV((interaction.chr2, interaction.pos2, interaction.end2), sv):
                if(parse_one(genome, pssm, (interaction.chr1, interaction.pos1, interaction.end1), (interaction.chr2, interaction.pos2, interaction.end2), 2, sv)):
                    pass
                else:
                    loops_altered += 1
                    sv_count.add(sv)
        elif(sv.svtype == "<DUP>"):
            results.append((interaction, sv))
            total_dup_intersect += 1
    print(loops_filename + " " + sv_filename + " Total to remove: " + str(loops_altered) + " SVs total: " + str(len(sv_count)) + "DUP total:" + str(total_dup_intersect))
    pass

def main():
    #sv_filename = "prediction/gm_svs.vcf"
    sv_filename = "prediction/HG00512_full.vcf"
    loops_filename = "prediction/GM12878.bedpe"
    #loops_filename = "prediction/HG00512.bedpe"
    start_time_total = time.time()
    predict_loops(sv_filename, loops_filename)
    print("--- Executed in %s seconds ---" % (time.time() - start_time_total))

if __name__ == "__main__":
    main()
