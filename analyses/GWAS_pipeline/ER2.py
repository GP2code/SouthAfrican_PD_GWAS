import gzip
outputFile = open("AllER2.txt", "w")
outputFile.write("CHR\tEmpRsq\n")
for chrom in range(1,23):
    info = gzip.open(f"chr{chrom}.info.gz")
    header = True
    for line in info:
        line = line.decode("utf-8")
        if header:
            if "#CHROM" in line:
                header = False
        else:
                if "ER2" in line:
                    infos = line.split()[-1].split(";")
                    for info in infos:
                        if "ER2" in info:
                            outputFile.write(f"{chrom}\t{info.replace('ER2=', '')}\n")
outputFile.close()
