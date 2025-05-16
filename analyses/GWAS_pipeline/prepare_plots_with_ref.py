fileSW = open("AutosomalPCA_VariantLoading.pcl")
fileDict = {}
for i in range(1,51):
    fileDict[f"PC{i}"] = open(f"PC{i}ToMan.txt", "w")
    fileDict[f"PC{i}"].write(f"SNP\tChromosome\tPosition\tPC{i}\n")
header = True
for line in fileSW:
    if header:
        header = False
    else:
        split = line.strip().split("\t")
        chrom, pos, A1, A2 = split[0].split(":")
        for i in range(1,51):
            if float(split[i+3]) < 0:
                split[i+3] = float(split[i+3])*-1
            fileDict[f"PC{i}"].write(f"{split[0]}\t{chrom}\t{pos}\t{split[i+3]}\n")
for i in range(1,51):
    fileDict[f"PC{i}"].close()
