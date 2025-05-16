# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 11:13:13 2024
@author: PEIXOTT
"""
fileSW = open("AutosomalPCA_VariantLoading.pcl")
fileDict = {}
for i in range(1,51):
    fileDict[f"PC{i}"] = open(f"PC{i}ToMan.txt", "w")
    fileDict[f"PC{i}"].write(f"SNP\tCHR\tBP\tPC\n")
header = True
count = 1
for line in fileSW:
    if header:
        header = False
    else:
        split = line.strip().split("\t")
        chrom, pos, A1, A2 = split[0].split(":")
        count = count+1
        for i in range(1,51):
            if float(split[i+3]) < 0:
                split[i+3] = float(split[i+3])*-1
            fileDict[f"PC{i}"].write(f"{split[0]}\t{chrom}\t{count}\t{split[i+3]}\n")
for i in range(1,51):
    fileDict[f"PC{i}"].close()
