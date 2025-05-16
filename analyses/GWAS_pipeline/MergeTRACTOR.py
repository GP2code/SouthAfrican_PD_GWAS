import scipy.stats
import argparse
def calculateP(beta, se):
    if (beta != "NA") and (se != "NA") and (float(se) != 0.0):
        z = float(beta)/float(se)
        return 1 - scipy.stats.norm.cdf(z), beta, se
    return -1, 0, 0
def addInDict(dictPval, idSNP, p, dataset):
    if idSNP not in dictPval:
        dictPval[idSNP] = {}
    dictPval[idSNP][dataset] = p
    return dictPval
#parser = argparse.ArgumentParser(description="P-value split")
#required = parser.add_argument_group("Mandatory arguments")
#required.add_argument('-c', '--chrom', help = "Chromomsome", required = True)
#args = parser.parse_args()
for j in ["Stepwise"]:
        dictPval = {}
        dictPval["TRACTOR-AFR"] = open(f"TRACTOR-AFR_hg38_P_{j}", "w")
        dictPval["TRACTOR-EUR"] = open(f"TRACTOR-EUR_hg38_P_{j}", "w")
        dictPval["TRACTOR-MALAY"] = open(f"TRACTOR-MALAY_hg38_P_{j}", "w")
        dictPval["TRACTOR-NAMA"] = open(f"TRACTOR-NAMA_hg38_P_{j}", "w")
        dictPval["TRACTOR-SAS"] = open(f"TRACTOR-SAS_hg38_P_{j}", "w")
        dictPval["TRACTOR-ALL"] = open(f"TRACTOR-ALL_hg38_P_{j}", "w")
        header = "Chrom\tPosition\tSNP\tRef\tAlt\tN\tBETA\tSE\tP"
        dictPval["TRACTOR-AFR"].write(f"{header}\n")
        dictPval["TRACTOR-EUR"].write(f"{header}\n")
        dictPval["TRACTOR-MALAY"].write(f"{header}\n")
        dictPval["TRACTOR-NAMA"].write(f"{header}\n")
        dictPval["TRACTOR-SAS"].write(f"{header}\n")
        dictPval["TRACTOR-ALL"].write(f"{header}\n")
        for i in range(1,23):
                print(f"{j} - {i}")
                IDPos = {}
                file = open(f"./PD_GWAS_SouthAfrica/admix-kit/PGEN_PC/stellen_{i}.pvar")
                header = True
                for line in file:
                        if header:
                                if line.startswith("#CHROM"):
                                        header = False
                        else:
                                CHROM, POS, ID, REF, ALT, INFO = line.strip().split()
                                IDPos[ID]= f"{CHROM}\t{POS}\t{ID}\t{REF}\t{ALT}"
                print(f"Open stellen_{j}_TRACTOR_{i}.TRACTOR.assoc")
                file = open(f"stellen_{j}_TRACTOR_{i}.TRACTOR.assoc")
                header = True
                for line in file:
                        if header:
                                header = False
                        else:
                                split = line.strip().split()
                                idSNP = split[0]
                                p = split[-1]
                                if p != "NA":
                                        #info = idSNP.replace("chr", "").split(":")
                                        #chrom = info[0]
                                        #pos38 = info[1]
                                        #A1 = info[2]
                                        #A2 = info[3]
                                        basicInfo = f"{IDPos[idSNP]}\t{split[-2]}"
                                        beta = split[1]
                                        se = split[2]
                                        if (beta != "NA") and (se != "NA") and (float(se) != 0.0):
                                                p, beta, se = calculateP(beta, se)
                                                dictPval["TRACTOR-AFR"].write(f"{basicInfo}\t{beta}\t{se}\t{p}\n")
                                        beta = split[3]
                                        se = split[4]
                                        if (beta != "NA") and (se != "NA") and (float(se) != 0.0):
                                                p, beta, se = calculateP(beta, se)
                                                dictPval["TRACTOR-EUR"].write(f"{basicInfo}\t{beta}\t{se}\t{p}\n")
                                        beta = split[5]
                                        se = split[6]
                                        if (beta != "NA") and (se != "NA") and (float(se) != 0.0):
                                                p, beta, se = calculateP(beta, se)
                                                dictPval["TRACTOR-MALAY"].write(f"{basicInfo}\t{beta}\t{se}\t{p}\n")
                                        if (beta != "NA") and (se != "NA") and (float(se) != 0.0):
                                                p, beta, se = calculateP(beta, se)
                                        
                                        dictPval["TRACTOR-NAMA"].write(f"{basicInfo}\t{beta}\t{se}\t{p}\n")
                                        if (beta != "NA") and (se != "NA") and (float(se) != 0.0):
                                                p, beta, se = calculateP(beta, se)
                                        dictPval["TRACTOR-SAS"].write(f"{basicInfo}\t{beta}\t{se}\t{p}\n")
                                        if (beta != "NA") and (se != "NA") and (float(se) != 0.0):
                                                p, beta, se = calculateP(beta, se)
                                                dictPval["TRACTOR-ALL"].write(f"{basicInfo}\t{split[5]}\t{split[6]}\t{split[-1]}\n")
                file.close()
        print("TRACTOR Done")
        dictPval["TRACTOR-AFR"].close()
        dictPval["TRACTOR-EUR"].close()
        dictPval["TRACTOR-MALAY"].close()
        dictPval["TRACTOR-NAMA"].close()
        dictPval["TRACTOR-SAS"].close()
        dictPval["TRACTOR-ALL"].close()

