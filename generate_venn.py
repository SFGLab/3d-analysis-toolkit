import venn
import matplotlib
import matplotlib.pyplot as plt
import math
import time
from sortedcontainers import SortedList
import itertools

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
            return self.chr1 == other.chr1 and self.pos1 == other.pos1 and self.end1 == other.end1 and self.chr2 == other.chr2 and self.pos2 == other.pos2 and self.end2 == other.end2
        return False
def checkOverlap(interaction1, interaction2):
    extension = -1
    if(interaction1.chr1 != interaction2.chr1 or interaction1.chr2 != interaction2.chr2):
        return False
    if(max(interaction1.pos1, interaction2.pos1)-extension <= min(interaction1.end1, interaction2.end1)):
        if(max(interaction1.pos2, interaction2.pos2)-extension <= min(interaction1.end2, interaction2.end2)):
            return True
    return False

def loadInteractions(fileName):
    interactions = SortedList([], key=lambda x: (x.chr1, x.pos1, x.end1, x.chr2, x.pos2, x.end2))
    with open(fileName, 'r') as f: #open the file
        lines = f.readlines()
        for line in lines:
            values = line.split("\t")
            interaction = Interaction(values[0], int(values[1]), int(values[2]), values[3], int(values[4]), int(values[5]), int(values[6]))
            interactions.add(interaction)
    return interactions

def getOverlappingMultiple(argv):
    if(len(argv) == 0):
        return list()
    interactions = argv[0]
    if(len(argv) > 1):
        for arg in argv[1:]:
            interactions = getOverlapping(interactions, arg)
    return interactions

def getOverlapping(interactions1, interactions2):
    interactions = list()
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
            interactions.append(Interaction(interaction1.chr1, min(interaction1.pos1, interaction2.pos1), max(interaction1.end1, interaction2.end1),
            interaction1.chr2, min(interaction1.pos2, interaction2.pos2), max(interaction1.end2, interaction2.end2), interaction1.pet+interaction2.pet))
    return interactions
start_time = time.time()

def compare_touples4(t1, t2):
    counter = 0
    for i in range(0,4):
        if(t1[i] == 1 and t1[i] == t2[i]):
            counter += 1
    return counter

def sum_tuple(t):
    c = 0
    for item in t:
        c += item
    return c
interactions_all = dict()

interactions_all["NA19238"] = loadInteractions("/mnt/raid/ctcf_prediction_anal/trios_new_ctcf/ctcf_named/output/NA19238.bedpe") # 1
interactions_all["NA19239"] = loadInteractions("/mnt/raid/ctcf_prediction_anal/trios_new_ctcf/ctcf_named/output/NA19239.bedpe") # 2
interactions_all["NA19240"] = loadInteractions("/mnt/raid/ctcf_prediction_anal/trios_new_ctcf/ctcf_named/output/NA19240.bedpe") # 3
interactions_all["GM12878"] = loadInteractions("/mnt/raid/ctcf_prediction_anal/trios_new_ctcf/ctcf_named/output/GM12878.bedpe") # 4
print("--- Loaded in %s seconds ---" % (time.time() - start_time))
start_time = time.time()

all_combinations = list(itertools.product([0, 1], repeat=4))
all_combinations = sorted(all_combinations, key=lambda x: sum(list(x)))
my_dict = dict()

for combination in all_combinations:
    interactions_list = list()
    if(combination[0] == 1):
        interactions_list.append(interactions_all["NA19238"])
    if(combination[1] == 1):
        interactions_list.append(interactions_all["NA19239"])
    if(combination[2] == 1):
        interactions_list.append(interactions_all["NA19240"])
    if(combination[3] == 1):
        interactions_list.append(interactions_all["GM12878"])
    my_dict[combination] = len(getOverlappingMultiple(interactions_list))
print("--- Combinations calculated in %s seconds ---" % (time.time() - start_time))
start_time = time.time()

dictlist = list()
for key, value in my_dict.items():
    temp = (key,value)
    dictlist.append(temp)
dictlist = sorted(dictlist, key=lambda x: sum(list(x[0])))
labels = dict()
for item in dictlist:
    to_remove = sum(x[1] for x in dictlist if sum_tuple(item[0]) == compare_touples4(item[0], x[0]) and item[0] != x[0])
    item = (item[0], item[1]-to_remove)
    labels["".join(map(str, item[0]))] = item[1]

fig, ax = venn.venn4(labels, names=['NA19238', 'NA19239', 'NA19240', 'GM12878'])
fig.savefig('venn4.png', bbox_inches='tight')
plt.close()

print("--- The rest done in %s seconds ---" % (time.time() - start_time))


#venn(interactions_all)
#plt.show()