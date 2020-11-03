import argparse
from sortedcontainers import SortedList
import progressbar

def checkOverlap(interaction1, interaction2):
        if(interaction1.chr1 != interaction2.chr1 or interaction1.chr2 != interaction2.chr2):
            return False

        if(max(interaction1.pos1, interaction2.pos1) <= min(interaction1.end1, interaction2.end1)):
            if(max(interaction1.pos2, interaction2.pos2) <= min(interaction1.end2, interaction2.end2)):
                return True
        return False

def loadInteractions(fileName):
    interactions = SortedList([], key=lambda x: (x.chr1, x.pos1, x.end1, x.chr2, x.pos2, x.end2))
    with open(fileName, 'r') as f: #open the file
        print("Parsing file: " + fileName+"...")
        lines = f.readlines()
        bar = progressbar.ProgressBar(maxval=len(lines), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
        i = 0
        bar.start()
        for line in lines:
            bar.update(i+1)
            i = i + 1
            values = line.split("\t")
            interaction = Interaction(values[0], int(values[1]), int(values[2]), values[3], int(values[4]), int(values[5]), int(values[6]))
            interactions.add(interaction)
        bar.finish()
        return interactions
class Interaction:
    def __init__(self, chr1, pos1, end1, chr2, pos2, end2, pet):
        self.chr1 = chr1
        self.chr2 = chr2
        self.pos1 = pos1
        self.pos2 = pos2
        self.end1 = end1
        self.end2 = end2
        self.pet = pet
    def __eq__(self, other):
        """Overrides the default implementation"""
        if isinstance(other, Interaction):
            return self.chr1 == other.chr1 and self.pos1 == other.pos1 and self.end1 == other.end1 and self.chr2 == other.chr2 and self.pos2 == other.pos2 and self.end2 == other.end2
        return False
    def generateLine(self):
        return self.chr1+"\t"+str(self.pos1)+"\t"+str(self.end1)+"\t"+self.chr2+"\t"+str(self.pos2)+"\t"+str(self.end2)+"\t"+str(self.pet)+"\n"

parser = argparse.ArgumentParser()
parser.add_argument('interactions_files', metavar='file', help="Files containing ", nargs='+')

args = parser.parse_args()

if(len(args.interactions_files) < 2):
    print("Merging one file makes no sense, does it?")
    exit(1)

interactionsReplicates = list()

for file in args.interactions_files:
    interactions = loadInteractions(file)
    interactionsReplicates.append(interactions)

allInteractions = interactionsReplicates[0]
fileId = 2
for interactionsReplicate in interactionsReplicates[1:]:
    print("Merging interactions from file #" + str(fileId))
    bar = progressbar.ProgressBar(maxval=len(interactionsReplicate), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    bar.start()
    i = 0
    for interaction in interactionsReplicate:
        bar.update(i + 1)
        i = i + 1
        beforeInteraction = allInteractions.bisect_left(Interaction(interaction.chr1, interaction.pos1-5000, interaction.pos1, interaction.chr1, interaction.pos1, interaction.pos1, interaction.pet))
        if(beforeInteraction > len(allInteractions)-1): # is after all, just to be sure
            allInteractions.add(interaction)
            continue
        afterInteraction = allInteractions.bisect_right(Interaction(interaction.chr1, interaction.end1+5000, interaction.end1, interaction.chr1, interaction.end1, interaction.end1, interaction.pet))
        
        found_interaction = False
        for interactionId in range(beforeInteraction, afterInteraction+1):
            if(checkOverlap(allInteractions[interactionId], interaction)): # check if overlaps on both sides
                allInteractions[interactionId].pet += interaction.pet
                # merge anchors
                edited += 1
                allInteractions[interactionId].pos1 = min(allInteractions[interactionId].pos1, interaction.pos1)
                allInteractions[interactionId].end1 = max(allInteractions[interactionId].end1, interaction.end1)

                allInteractions[interactionId].pos2 = min(allInteractions[interactionId].pos2, interaction.pos2)
                allInteractions[interactionId].end2 = max(allInteractions[interactionId].end2, interaction.end2)

                found_interaction = True
                break
        if not found_interaction: # it does not overlap with any, lets add it as independent one
            allInteractions.add(interaction)
    bar.finish()
with open('result.bedpe', 'w') as f:
    for interaction in allInteractions:
        f.write(interaction.generateLine())

print("Finished parsing!")