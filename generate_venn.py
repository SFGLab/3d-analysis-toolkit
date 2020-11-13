import venn
import matplotlib
import matplotlib.pyplot as plt
import math
import time
from sortedcontainers import SortedList
import itertools
import asyncio
import multiprocessing as mp

class Interaction:
    def __init__(self, chr1, pos1, end1, chr2, pos2, end2, pet):
        self.chr1 = chr1
        self.chr2 = chr2
        self.pos1 = pos1
        self.pos2 = pos2
        self.end1 = end1
        self.end2 = end2
        self.pet = pet
    def generateLine(self):
        return self.chr1+"\t"+str(self.pos1)+"\t"+str(self.end1)+"\t"+self.chr2+"\t"+str(self.pos2)+"\t"+str(self.end2)+"\t"+str(self.pet)+"\n"
    def __eq__(self, other):
        """Overrides the default implementation"""
        if isinstance(other, Interaction):
            return checkOverlap(self, other)
        return False
    def __hash__(self):
        return hash((self.chr1, self.chr2, self.pos1, self.pos2, self.end1, self.end2))
def checkOverlap(interaction1, interaction2):
    extension = -1 # -1 added for actuall 1bp overlap
    if(interaction1.chr1 != interaction2.chr1 or interaction1.chr2 != interaction2.chr2):
        return False
    if(max(interaction1.pos1, interaction2.pos1)-extension <= min(interaction1.end1, interaction2.end1)):
        if(max(interaction1.pos2, interaction2.pos2)-extension <= min(interaction1.end2, interaction2.end2)):
            return True
    return False

def orderInt(x):
    return (x.chr1, x.pos1, x.end1, x.chr2, x.pos2, x.end2)

def loadInteractions(fileName):
    interactions = SortedList([], key=orderInt)
    with open(fileName, 'r') as f: #open the file
        lines = f.readlines()
        for line in lines:
            values = line.split("\t")
            interaction = Interaction(values[0], int(values[1]), int(values[2]), values[3], int(values[4]), int(values[5]), int(values[6]))
            interactions.add(interaction)
    return interactions

def getSet(args):
    interactions_all, comb = args
    interactions = SortedList([], key=lambda x: (x.chr1, x.pos1, x.end1, x.chr2, x.pos2, x.end2))
    i = 0
    # needed twice because first intersect, then removal
    for c_i in comb:
        if(c_i == 1):
            interactions = getOverlapping(interactions, interactions_all[i])
        i += 1
    i = 0
    for c_i in comb:
        if(c_i == 0):
            interactions = removeOverlapping(interactions, interactions_all[i])
        i += 1

    seen = set()
    seen_add = seen.add
    return ("".join(map(str, comb)), [x for x in interactions if not (x in seen or seen_add(x))])

def removeOverlapping(interactions1, interactions2): # subtract I1-I2
    interactions = SortedList([], key=lambda x: (x.chr1, x.pos1, x.end1, x.chr2, x.pos2, x.end2))
    for interaction in interactions1:
        interactions.add(interaction)
    to_remove = getOverlapping(interactions, interactions2, True)
    for interaction in to_remove:
        interactions.discard(interaction)
    return interactions

def getOverlapping(interactions1, interactions2, toRemove=False): # intersection
    interactions = SortedList([], key=lambda x: (x.chr1, x.pos1, x.end1, x.chr2, x.pos2, x.end2))
    if(toRemove):
        interactions = list()
    if(len(interactions1) == 0):
        return interactions2
    for interaction1 in interactions1:
        beforeInteraction = interactions2.bisect_left(Interaction(interaction1.chr1, interaction1.pos1-20000, interaction1.pos1, interaction1.chr1, interaction1.pos1, interaction1.pos1, interaction1.pet))
        afterInteraction = interactions2.bisect_right(Interaction(interaction1.chr1, interaction1.end1+20000, interaction1.end1, interaction1.chr1, interaction1.end1, interaction1.end1, interaction1.pet))

        for interaction2 in interactions2[beforeInteraction:min(afterInteraction, len(interactions2))]:
            if(interaction2.chr1 < interaction1.chr1):
                continue
            if(interaction2.chr1 > interaction1.chr1):
                break
            if not checkOverlap(interaction1, interaction2):
                continue
            if not(toRemove):
                interactions.add(Interaction(interaction1.chr1, interaction1.pos1, interaction1.end1,
            interaction1.chr2, interaction1.pos2, interaction1.end2, interaction1.pet+interaction2.pet))
            else:
                interactions.append(interaction1)
            #interactions.add(Interaction(interaction1.chr1, min(interaction1.pos1, interaction2.pos1), max(interaction1.end1, interaction2.end1),
            #interaction1.chr2, min(interaction1.pos2, interaction2.pos2), max(interaction1.end2, interaction2.end2), interaction1.pet+interaction2.pet))
    return interactions
start_time = time.time()


start_time = time.time()
interactions_all = dict()

interactions_all[0] = loadInteractions("/mnt/raid/ctcf_prediction_anal/trios_new_ctcf/ctcf_named/output/NA19238.bedpe") # 1
interactions_all[1] = loadInteractions("/mnt/raid/ctcf_prediction_anal/trios_new_ctcf/ctcf_named/output/NA19239.bedpe") # 2
interactions_all[2] = loadInteractions("/mnt/raid/ctcf_prediction_anal/trios_new_ctcf/ctcf_named/output/NA19240.bedpe") # 3
interactions_all[3] = loadInteractions("/mnt/raid/ctcf_prediction_anal/trios_new_ctcf/ctcf_named/output/GM12878.bedpe") # 4
print("--- Loaded in %s seconds ---" % (time.time() - start_time))
start_time = time.time()

all_combinations = list(itertools.product([0, 1], repeat=4))
all_combinations = sorted(all_combinations, key=lambda x: sum(list(x)))

threads = 16
task_list = list()
for combination in all_combinations:
    task_list.append((interactions_all, combination))
pool = mp.Pool(threads)
results = pool.map(getSet, task_list)
pool.close()
pool.join() 

labels = dict()
for result in results:
    labels[result[0]] = str(len(result[1]))

print("--- Combinations calculated in %s seconds ---" % (time.time() - start_time))
start_time = time.time()

fig, ax = venn.venn4(labels, names=['NA19238', 'NA19239', 'NA19240', 'GM12878'])
fig.savefig('venn4.png', bbox_inches='tight')
plt.close()

print("--- The rest done in %s seconds ---" % (time.time() - start_time))


#venn(interactions_all)
#plt.show()