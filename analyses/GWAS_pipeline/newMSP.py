import argparse
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run ADMIX-KIT')
    requiredGeneral = parser.add_argument_group("Required arguments")
    requiredGeneral.add_argument('-p', '--psam', help='PSAM file', required=True)
    requiredGeneral.add_argument('-m', '--msp', help='MSP file', required=True)
    requiredGeneral.add_argument('-o', '--output', help='Name of output folder', required=True)
    args = parser.parse_args()
    file = open(f"{args.psam}")
    header = True
    sampleList = []
    for line in file:
        if header:
            header = False
        else:
            split = line.strip().split()
            sampleList.append(split[0])
    file.close()
    fileIn = open(args.msp)
    fileOut = open(args.output, "w")
    dictIDs = {}
    header = True
    for line in fileIn:
        if header:
            if "#chm" in line:
                IDs = line.strip().split("\t")
                for i in range(6, len(IDs),2):
                    IDComponent = IDs[i].split(".")
                    IDNoHap = IDComponent[0]
                    for j in range(1, len(IDComponent)-1):
                        IDNoHap = IDNoHap+f".{IDComponent[j]}"
#                    print(f"{IDs[i]} - {IDNoHap}")
                    dictIDs[IDNoHap] = i
                #Header of new MSP
                fileOut.write(f"{IDs[0]}")
                for i in range(1,6):
                    fileOut.write(f"\t{IDs[i]}")
                for sample in sampleList:
                    fileOut.write(f"\t{sample}.0\t{sample}.1")
                fileOut.write("\n")
                header = False
            else:
                fileOut.write(line)
        else:
            data = line.strip().split("\t")
            fileOut.write(f"{data[0].replace('chr', '')}")
            for i in range(1,6):
                fileOut.write(f"\t{data[i]}")
            for sample in sampleList:
                index = dictIDs[sample]
                fileOut.write(f"\t{data[index]}\t{data[index+1]}")
            fileOut.write("\n")
    fileOut.close()
for ID in dictIDs:
        print(f"{ID}: {dictIDs[ID]}")

