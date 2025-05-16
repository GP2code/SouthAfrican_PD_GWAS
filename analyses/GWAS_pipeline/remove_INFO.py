import os
os.system("mkdir ./fixed")
for i in range(1,23):
	fileOld = open(f"chr{i}.pvar")
	fileNew = open(f"./fixed/chr{i}.pvar", "w")
	
	os.system(f"cp chr{i}.pgen ./fixed")
	os.system(f"cp chr{i}.psam ./fixed")
	
	header = True
	for line in fileOld:
		if header:
			if "#CHROM" in line:
				header = False
				split = line.strip().split()
				fileNew.write(f"{split[0]}")
				for j in range(1, 5):
					fileNew.write(f"\t{split[j]}")
				fileNew.write(f"\n")
		else:
			split = line.strip().split()
			fileNew.write(f"{split[0]}")
			for j in range(1, 5):
				fileNew.write(f"\t{split[j]}")
			fileNew.write(f"\n")
	fileNew.close()
	fileOld.close()
