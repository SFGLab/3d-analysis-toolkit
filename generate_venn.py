import venn
import matplotlib
import matplotlib.pyplot as plt
import math
import time
from sortedcontainers import SortedList
import itertools
import asyncio
import multiprocessing as mp
from common import Interaction, createFolder, removeFolder, loadInteractions, checkOverlap

class VennInteraction(Interaction):
    def __eq__(self, other):
        """Overrides the default implementation"""
        if isinstance(other, Interaction):
            return checkOverlap(self, other)
        return False
    def __hash__(self):
        return hash((self.chr1, self.chr2, self.pos1, self.pos2, self.end1, self.end2))

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

interactions_all[0] = loadInteractions("/mnt/raid/ctcf_prediction_anal/GM_comparisons_tries/GM12878_R1.bedpe", VennInteraction)
interactions_all[1] = loadInteractions("/mnt/raid/ctcf_prediction_anal/GM_comparisons_tries/GM12878_R2.bedpe", VennInteraction)

#interactions_all[0] = loadInteractions("/mnt/raid/ctcf_prediction_anal/trios_new_ctcf/ctcf_named/output/NA19238.bedpe") # 1
#interactions_all[1] = loadInteractions("/mnt/raid/ctcf_prediction_anal/trios_new_ctcf/ctcf_named/output/NA19239.bedpe") # 2
#interactions_all[2] = loadInteractions("/mnt/raid/ctcf_prediction_anal/trios_new_ctcf/ctcf_named/output/NA19240.bedpe") # 3
#interactions_all[3] = loadInteractions("/mnt/raid/ctcf_prediction_anal/trios_new_ctcf/ctcf_named/output/GM12878.bedpe") # 4
print("--- Loaded in %s seconds ---" % (time.time() - start_time))
start_time = time.time()

all_combinations = list(itertools.product([0, 1], repeat=len(interactions_all)))
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

fig, ax = venn.venn2(labels, names=['GM12878_R1', 'GM12878_R2'])
fig.savefig('venn2.png', bbox_inches='tight')
#fig, ax = venn.venn4(labels, names=['GM19238', 'GM19239', 'GM19240', 'GM12878'])
#fig.savefig('venn4.png', bbox_inches='tight')
plt.close()

print("--- The rest done in %s seconds ---" % (time.time() - start_time))


#venn(interactions_all)
#plt.show()