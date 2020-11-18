import venn
import matplotlib
import matplotlib.pyplot as plt
import math
import time
from sortedcontainers import SortedList
import itertools
import asyncio
from multiprocessing.pool import ThreadPool
from common import Interaction, createFolder, removeFolder, loadInteractions, checkOverlap, getOverlapping, removeOverlapping

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





start_time = time.time()


start_time = time.time()
interactions_all = dict()

samples = ["/mnt/raid/ctcf_prediction_anal/GM_comparisons_tries/GM12878_R1.bedpe", "/mnt/raid/ctcf_prediction_anal/GM_comparisons_tries/GM12878_R2.bedpe"]
#samples = list("/mnt/raid/ctcf_prediction_anal/trios_new_ctcf/ctcf_named/output/NA19238.bedpe", "/mnt/raid/ctcf_prediction_anal/trios_new_ctcf/ctcf_named/output/NA19239.bedpe", 
#"/mnt/raid/ctcf_prediction_anal/trios_new_ctcf/ctcf_named/output/NA19240.bedpe", "/mnt/raid/ctcf_prediction_anal/trios_new_ctcf/ctcf_named/output/GM12878.bedpe")

i = 0
for sample in samples:
    interactions_all[i] = loadInteractions(sample, VennInteraction)
    i += 1

print("--- Loaded in %s seconds ---" % (time.time() - start_time))
start_time = time.time()

all_combinations = list(itertools.product([0, 1], repeat=len(interactions_all)))
all_combinations = sorted(all_combinations, key=lambda x: sum(list(x)))

threads = 16
task_list = list()
for combination in all_combinations:
    task_list.append((interactions_all, combination))
pool = ThreadPool(threads)
results = pool.map(getSet, task_list)
pool.close()
pool.join() 

labels = dict()
for result in results:
    labels[result[0]] = str(len(result[1]))

print("--- Combinations calculated in %s seconds ---" % (time.time() - start_time))
start_time = time.time()


fig, ax = getattr(venn, 'venn'+str(len(interactions_all)))(labels, names=list(map(lambda x: x.split("/")[-1].split(".")[0], samples)))
fig.savefig('venn.png', bbox_inches='tight')
#fig, ax = venn.venn4(labels, names=['GM19238', 'GM19239', 'GM19240', 'GM12878'])
#fig.savefig('venn4.png', bbox_inches='tight')
plt.close()

print("--- The rest done in %s seconds ---" % (time.time() - start_time))


#venn(interactions_all)
#plt.show()