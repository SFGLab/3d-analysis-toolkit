def prepare_file_anchors(fileName):
    anchors = list()
    pet3 = list()
    pet4 = list()
    sample_name = fileName.split(".")[0]
    with open(fileName, 'r') as f: #open the file
        lines = f.readlines()
        for line in lines:
            values = line.split("\t")
            if(int(values[6]) < 4):
                pet3.append(values[0]+"\t"+values[1]+"\t"+values[2]+"\t"+values[3]+"\t"+values[4]+"\t"+values[5]+"\t1\n")
            else:
                anchors.append(values[0]+"\t"+values[1]+"\t"+values[2]+"\tR\n")
                anchors.append(values[3]+"\t"+values[4]+"\t"+values[5]+"\tL\n")
                pet4.append(line)
    with open(sample_name+"_anchors.bedpe", 'w') as f:
        for line in anchors:
            f.write(line)
    with open(sample_name+"_singletons.bedpe", 'w') as f:
        for line in pet3:
            f.write(line)
    with open(sample_name+"_clusters.bedpe", 'w') as f:
        for line in pet4:
            f.write(line)

def prepare_file_segments(fileName):
    segments = list()
    sample_name = fileName.split(".")[0]
    with open(fileName, 'r') as f: #open the file
        lines = f.readlines()
        for line in lines:
            values = line.split("\t")
            segments.append(values[0]+"\t"+values[1]+"\t"+values[1]+"\n")
            val2 = values[2].split("\n")[0]
            segments.append(values[0]+"\t"+val2+"\t"+val2+"\n")
    with open(sample_name+"_segments.bed", 'w') as f:
        for line in segments:
            f.write(line)

def main():
    prepare_file_anchors("WTC11.bedpe")
    prepare_file_segments("WTC11.bed")

if __name__ == "__main__":
    main()
