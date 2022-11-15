from common import loadInteractions, saveFile, orderInt
from Bio import SeqIO, motifs
import numpy as np
from sortedcontainers import SortedList
import time 
import vcf
import math

class SNP():
    def __init__(self, record=None, chrm=None, pos=None):
        if(record):
            self.chr = record.CHROM.split("chr")[1]
            if(len(record.ALT) > 1):
                self.alleles = [record.ALT[0].sequence, record.ALT[1].sequence]
            elif(record.nucl_diversity > 0.99):
                self.alleles = [record.REF, record.ALT[0].sequence]
            elif(record.nucl_diversity < 0.01):
                self.alleles = [record.ALT[0].sequence, record.ALT[0].sequence]

            self.pos = record.POS
        else:
            self.chr = chrm.split("chr")[1]
            self.pos = pos
        try:
            self.chr = int(self.chr)
        except:
            self.chr = ord(self.chr)
    def __lt__(self, other):
        return self.chr < other.chr or (self.chr == other.chr and self.pos < other.pos)

def updateSequence(seq, snps, toRemove, id):
    for snp in snps:
        relativePos = snp.pos-toRemove
        seq = seq[0:relativePos] + snp.alleles[id] + seq[relativePos+1:]
    return seq

def decideOne(pssm, seq1, seq2, probabilistic=False):
    with np.errstate(invalid='ignore'):
        if not(probabilistic):
            search_results1 = [(x,y) for x,y in pssm.search(seq1, threshold=7.0)]
            search_results2 = [(x,y) for x,y in pssm.search(seq2, threshold=7.0)]
        else:
            search_results1 = [(x,y) for x,y in pssm.search(seq1, threshold=4.0)]
            search_results2 = [(x,y) for x,y in pssm.search(seq2, threshold=4.0)]
            res1_l = 0
            res1_r = 0
            res2_l = 0
            res2_r = 0
            for r1 in search_results1:
                if(r1[0] >= 0): # right
                    res1_r += math.pow(2, r1[1])
                else: # left
                    res1_l += math.pow(2, r1[1])
            for r2 in search_results2:
                if(r2[0] >= 0): # right
                    res2_r += math.pow(2, r2[1])
                else: # left
                    res2_l += math.pow(2, r2[1])

            if(res1_l+res1_r == 0):
                res1 = -1 # no motif
            else:
                res1 = res1_l / (res1_l+res1_r) # res1 = % of strength in motif to the left             
            
            if(res2_l+res2_r == 0):
                res2 = -1 # no motif
            else:
                res2 = res2_l / (res2_l+res2_r) # res2 = % of strength in motif to the left
                
            return (res1, res2)
    for r1 in search_results1:
        for r2 in search_results2:
            if(r1[0] >= 0 and r2[0] < 0): # >........<
                return True
            if(r1[1] >= 9.0 and r2[1] >= 9.0):
                if(r1[0]*r2[0] >= 0): # >.........> or <.........<
                    return True
            if(r1[1] >= 12.0 and r2[1] >= 12.0):
                if(r1[0] < 0 and r2[0] >= 0): # <.........>
                    return True
    return False
def decideInteraction(pssm, seq1, seq2, snps1=None, snps2=None, toRemove1=None, toRemove2=None, probabilistic=False):
    if not(snps1):
        decision = decideOne(pssm, seq1, seq2, probabilistic)
        if(probabilistic):
            return decision
        if(decision):
            return True
        return False
    else:
        # maternal / first -> 1|1
        seq1_left = updateSequence(seq1._data, snps1, toRemove1, 0)
        seq2_left = updateSequence(seq2._data, snps2, toRemove2, 0)
        decision = decideOne(pssm, seq1_left, seq2_left)
        
        if(probabilistic):
            return decision
        if(decision):
            return True
        # paternal / second 1|1 <-
        seq1_right = updateSequence(seq1._data, snps1, toRemove1, 1)
        seq2_right = updateSequence(seq2._data, snps2, toRemove2, 1)
        decision = decideOne(pssm, seq1_right, seq2_right)
        
        if(probabilistic):
            return decision
        if(decision):
            return True
        return False

def filterInteractionsByMotifs(interactions, probabilistic=False, vcfFile=None):
    genome = dict()

    pfmMatrixFile = "/mnt/raid/repos/3d-analysis-toolkit/MA0139.1.pfm"
    start_time = time.time()
    if(vcfFile):
        snps = SortedList()
        vcf_reader = vcf.Reader(filename=vcfFile)
        for record in vcf_reader:
            if "chr" in record.CHROM and not("_" in record.CHROM) and not("*" in record.CHROM) and not("EBV" in record.CHROM) and not(record.is_indel):
                snps.add(SNP(record))
        print("--- Loading SNPs executed in %s seconds ---" % (time.time() - start_time))

    for record in SeqIO.parse("/mnt/raid/GRCh38_full_analysis_set_plus_decoy_hla.fa", "fasta"):
        if("_" in record.id or "*" in record.id or "EBV" in record.id):
            continue
        genome[record.id] = record.seq

    extension = 0
    with open(pfmMatrixFile) as handle:
        ctcf_motif = motifs.read(handle, "pfm")

    pssm = ctcf_motif.pssm

    interactions_with_motif = SortedList([], key=orderInt)
    for interaction in interactions:
        sequence1 = genome[interaction.chr1][interaction.pos1-extension:interaction.end1+extension]
        sequence2 = genome[interaction.chr2][interaction.pos2-extension:interaction.end2+extension]
        if(vcfFile):
            all_snps_seq1 = snps[snps.bisect_left(SNP(None, interaction.chr1, interaction.pos1)):snps.bisect_right(SNP(None, interaction.chr1, interaction.end1))]
            all_snps_seq2 = snps[snps.bisect_left(SNP(None, interaction.chr2, interaction.pos2)):snps.bisect_right(SNP(None, interaction.chr2, interaction.end2))]
            decision = decideInteraction(pssm, sequence1, sequence2, all_snps_seq1, all_snps_seq2, interaction.pos1+1, interaction.pos2+1, probabilistic)
            if(probabilistic):
                interaction.prob1 = decision[0]
                interaction.prob2 = decision[1]
                interactions_with_motif.add(interaction)
                continue
            if(decision):
                interactions_with_motif.add(interaction)
        else:
            decision = decideInteraction(pssm, sequence1, sequence2, probabilistic=probabilistic)
            if(probabilistic):
                interaction.prob1 = decision[0]
                interaction.prob2 = decision[1]
                interactions_with_motif.add(interaction)
                continue
            if(decision):
                interactions_with_motif.add(interaction)

    return interactions_with_motif

def main():
    maxLength = 500000

    start_time = time.time()

    probabilistic = True
    
    interactionsFile = "/home/mchilinski/Downloads/CTCF_hg38_PeakSupport.bedpe"
    vcfFile = None # "/mnt/raid/trios_data/sv-callings/sv_callings_data/snp_trios/SNP_HG00512.vcf"

    interactions = loadInteractions(interactionsFile, maxLength=0)
    print("All interactions:" + str(len(interactions)))
    interactions_with_motif = filterInteractionsByMotifs(interactions, probabilistic, vcfFile)

    print("Interactions with motif: " + str(len(interactions_with_motif)))
    print("--- Executed in %s seconds ---" % (time.time() - start_time))
    saveFile(interactionsFile.split(".")[0]+"_2.bedpe", interactions_with_motif, add_prob=True)

if __name__ == "__main__":
    main()