import os
import gzip

#changing directory to parent directory of all files
os.chdir("./PD_GWAS_SouthAfrica/")

#file paths for external summary stats files
fileKim = "./PD_GWAS_SumStats/Kim2023/supplementaryTable3.txt"
fileLARGE = "./PD_GWAS_SumStats/Loesch2021/large.assoc_single.ALL.04302020_hg38.txt"
fileRizig = "./PD_GWAS_SumStats/Rizig2023/Rizig_et_al_2023_AFR_AAC_metaGWAS_no23andMe_hg38.txt"
fileNalls = "./PD_GWAS_SumStats/Nalls2019/supplementaryTable2.txt"

#file path for Step SAPDSC GWAS summary stats
dictStepConventional = {}
fileStepConventional = "./PD_GWAS_SumStats/Step2025/sorted_new_results_including_malay_nama_hg38.txt"

#reading Step SAPDSC GWAS summary stats file, creating a dictionary, and extracting beta values and p-vals
header = True
with open(fileStepConventional,'r') as fileStep:
    for line in fileStep:
        #skipping all header lines
        if header:
            if line.startswith("CHR"):
                header = False
        else:
            split = line.strip().split()
            chrom = split[0].replace("chr","")
            bp = split[1]
            other_allele = split[3].upper()
            effect_allele = split[4].upper()
            beta = split[5]
            p_val = split[7]
            if f"{chrom}:{bp}" not in dictStepConventional:
                #defining dictionary based on chromosome and position
                dictStepConventional[f"{chrom}:{bp}"] = {}
                #including the ref and other alleles in the dictionary to compare to those from other studies
                dictStepConventional[f"{chrom}:{bp}"]['full'] = f"{chrom}:{bp}:{other_allele}:{effect_allele}"
                #extract beta and p values
                dictStepConventional[f"{chrom}:{bp}"]['beta'] = float(beta)
                dictStepConventional[f"{chrom}:{bp}"]['p-val'] = float(p_val)
print("Done reading conventional GWAS Step file.")

#file path for SAPDSC GWAS summary stats
fileStepATT = "./PD_GWAS_SumStats/Step2025/ATT/merged_summary_stats_with_chr_BP_REF_ALT.txt.gz"
dictStepATT = {}
dictStepATTSuggestive = {}

#reading Step SA GWAS summary stats file, creating a dictionary, and extracting beta values and p-vals
with gzip.open(fileStepATT, 'r') as openfileStep:
    header = True
    for line in openfileStep:
        line = line.decode('utf-8')
        # skipping all header lines
        if header:
            if line.startswith("Chrom"):
                header = False
        else:
            chr, rsID, other_allele, effect_allele, bp, beta, se, n, p_val = line.strip().split()
            chrom = chr.replace("chr","")
            effect = effect_allele.upper()
            other = other_allele.upper()
            if f"{chrom}:{bp}" not in dictStepATT:
                #defining dictionary based on chromosome and position
                dictStepATT[f"{chrom}:{bp}"] = {}
                #including the ref and other alleles in the dictionary to compare to those from other studies
                dictStepATT[f"{chrom}:{bp}"]['full'] = f"{chrom}:{bp}:{other}:{effect}"
                #extract beta and p values
                dictStepATT[f"{chrom}:{bp}"]['beta'] = float(beta)
                dictStepATT[f"{chrom}:{bp}"]['p-val'] = float(p_val)
print("Done reading LA-GWAS with ATT Step file.")

#defineLeadSNPs parameters:
#   leadSNPs: file name of the list of lead SNPs
#   chromsome: index of chromosome column
#   basepair: index of position column
def defineLeadSNPs(leadSNPs, chromosome, basepair):
    topVars = open(leadSNPs)
    dictionary = {}
    header = True
    for line in topVars:
        if header:
            header = False
        else:
            split = line.strip().split("\t")
            chrom = split[chromosome].replace("chr","")
            bp = split[basepair]
            if f"{chrom}:{bp}" not in dictionary:
                #writing dictionary of lead SNPs
                dictionary[f"{chrom}:{bp}"] = {}
    print(f"{leadSNPs}: {len(dictionary)}")
    return dictionary

#creating dictionaries of all lead SNPs for external studies and SA GWAS
dictStepConLead = defineLeadSNPs("./PD_GWAS_SumStats/Step2025/leadSNPs_hg38", 0, 1)
dictStepATTLead = defineLeadSNPs("./PD_GWAS_SumStats/Step2025/ATT/leadSNPs.txt", 4, 5)
dictKim = defineLeadSNPs("./PD_GWAS_SumStats/Kim2023/leadSNPs_fromJeff_hg38.txt", 1, 2)
dictLoesch = defineLeadSNPs("./PD_GWAS_SumStats/Loesch2021/leadSNPs_hg38.txt", 4, 5)
dictNalls = defineLeadSNPs("./PD_GWAS_SumStats/Nalls2019/Nalls_topHits.txt", 1, 2)
dictRizig = defineLeadSNPs("./PD_GWAS_SumStats/Rizig2023/leadSNPs_hg38.txt", 1, 2)
print("Done extracting lead SNPs.")

#defining output file path and header for combined summary stats for conventional GWAS lead SNPs
fileOutConventional = open("./PD_GWAS_SumStats/Step2025/outputSummaryStats/Combined_SummaryStats_LeadSNPs.tsv","w")
fileOutConventional.write("varID\tbeta_ext\tbeta_Step\tp-val_ext\tp-val_Step\tleadSource\tstudy\tleadInStep?\n")

#defining output file for locations that had different other and/or effect alleles from those in Step summary stats
fileInconsistentConventional = open("./PD_GWAS_SumStats/Step2025/outputSummaryStats/inconsistentLeadSNPs.tsv","w")

#defining output file path and header for combined summary stats for ATT results
fileOutATT = open("./PD_GWAS_SumStats/Step2025/outputSummaryStats/ATT_Combined_SummaryStats.tsv","w")
fileOutATT.write("varID\tbeta_ext\tbeta_Step\tp-val_ext\tp-val_Step\tleadSource\tstudy\tleadInStep?\n")

#defining output file for locations that had different reference and/or other alleles from those in Step summary stats for ATT
fileInconsistentATT = open("./PD_GWAS_SumStats/Step2025/outputSummaryStats/ATT_inconsistentRefOther.tsv","w")

#checkPositionsAddBetaPval parameter guide
#   other_allele: other allele
#   effect_allele: effect allele
#   beta: beta
#   p_val: p-value
#   chr_bp: position in format chr:pos
#   study: name of study to be diplayed later in plots ex. "Loesch et al. 2021"
#   leadSource: "yes"-lead SNP from SAPDSC, "no"-lead SNP in previous study
#   count: count betas inverted, some duplicates may occur if SAPDSC and external studies have overlapping lead SNPs
#   dictType: type of GWAS used (ex. dictStepConventional or dictStepATT)
#   outFile: output file
#   mismatchFile: file for inconsistencies
def checkPositionsAddBetaPval(other_allele, effect_allele, beta, p_val, chr_bp, study, leadSource, count, dictType, outFile, mismatchFile):
    #checking if the position is in the summary stats for Step SAPDSC GWAS
    if chr_bp in dictType:
        #checking if the alleles align between Step SAPDSC and external summary stats
        if f"{chr_bp}:{other_allele}:{effect_allele}" == dictType[chr_bp]['full']:
            outFile.write(f"{chr_bp}:{other_allele}:{effect_allele}\t{float(beta)}\t{dictType[chr_bp]['beta']}\t{float(p_val)}\t{dictType[chr_bp]['p-val']}\t{leadSource}\t{study}\t{leadSource}\n")
        #checking if the alleles are reversed between Step SAPDSC and external summary stats and inverting the beta
        elif f"{chr_bp}:{effect_allele}:{other_allele}" == dictType[chr_bp]['full']:
            outFile.write(f"{chr_bp}:{effect_allele}:{other_allele}\t{float(beta)*(-1)}\t{dictType[chr_bp]['beta']}\t{float(p_val)}\t{dictType[chr_bp]['p-val']}\t{leadSource}\t{study}\t{leadSource}\n")
            count += 1
        #writing any positions that were in both but had conflicting alleles reported into fileInconsistent
        else:
            mismatchFile.write(f"Not adding {chr_bp}:{other_allele}:{effect_allele} from {study} to file because effect_allele and ref_allele are different from Kate's ({dictType[chr_bp]['full']})\n")
    return count

#extractBetaPVal parameter guide
#   file: file name of summary stats file
#   chromosome: column number of the chromosome
#   basepair: column number of the bp
#   other: column number of the other allele
#   effect: column number of the effect allele
#   betaval: column number of the beta
#   pval: column number of the p-val
#   dictionary: dictionary of leadSNPs
#   study: name of study to be diplayed later in plots ex. "Loesch et al. 2021"
#   dictType: type of GWAS used (ex. dictStepConventional or dictStepATT)
#   outFile: output file
#   mismatchFile: file for inconsistencies
#   dictStepLead: dictionary name for lead SNPs (ex. dictStepConLead or dictStepATTLead)
def extractBetaPVal(file, chromosome, basepair, other, effect, betaval, pval, dictionary, study, dictType, outFile, mismatchFile, dictStepLead):
    print(study)
    count = 0
    with open(file,'r') as openfile:    
        header = True
        for line in openfile:
            if header:
                header=False
            else:
                split = line.strip().split('\t')
                chrom = split[chromosome].replace("chr","")
                bp = split[basepair]
                other_allele = split[other].upper()
                effect_allele = split[effect].upper()
                beta = split[betaval]
                p_val = split[pval]
                chr_bp = f"{chrom}:{bp}"
                #adding variants that were lead SNPs in external studies
                #checking if the SNP in the external study is in the dictionary of lead SNPs
                if dictStepLead == dictStepConLead:
                    if chr_bp in dictionary:
                        count = checkPositionsAddBetaPval(other_allele, effect_allele, beta, p_val, chr_bp, study, "no", count, dictType, outFile, mismatchFile)
                #checking if the SNP in the external study is in the dictionary of lead SNPs
                if dictStepLead == dictStepATTLead:
                    if chr_bp in dictionary:
                        count = checkPositionsAddBetaPval(other_allele, effect_allele, beta, p_val, chr_bp, study, "no", count, dictType, outFile, mismatchFile)
                #adding variants that were lead SNPs in Step SAPDSC GWAS
                if chr_bp in dictStepLead:
                    count = checkPositionsAddBetaPval(other_allele, effect_allele, beta, p_val, chr_bp, study, "yes", count, dictType, outFile, mismatchFile)
    return count

# Conventional GWAS
#extracting beta and p-vals for external studies
countKim = extractBetaPVal(fileKim, 3, 4, 6, 5, 12, 14, dictKim, "Kim et al. 2024", dictStepConventional, fileOutConventional, fileInconsistentConventional, dictStepConLead)
countLARGE = extractBetaPVal(fileLARGE, 2, 3, 4, 5, 15, 14, dictLoesch, "Loesch et al. 2021", dictStepConventional, fileOutConventional, fileInconsistentConventional, dictStepConLead)
countRizig = extractBetaPVal(fileRizig, 0, 1, 2, 3, 4, 7, dictRizig, "Rizig et al. 2023", dictStepConventional, fileOutConventional, fileInconsistentConventional, dictStepConLead)
countNalls = extractBetaPVal(fileNalls, 1, 2, 6, 5, 8, 10, dictNalls, "Nalls et al. 2019", dictStepConventional, fileOutConventional, fileInconsistentConventional, dictStepConLead)
         
fileOutConventional.close()
fileInconsistentConventional.close()
#printing number of positions in which the other and effect were inverted (and therefore the beta value was inverted manually)
print(f"Number of positions that betas were inverted for conventional GWAS lead SNPs, Kim: {countKim}, Loesch: {countLARGE}, Rizig: {countRizig}, Nalls: {countNalls}")

# LA-GWAS with ATT
#extracting beta and p-vals for external studies
countLARGE = extractBetaPVal(fileLARGE, 2, 3, 5, 4, 15, 14, dictLoesch, "Loesch et al. 2021", dictStepATT, fileOutATT, fileInconsistentATT, dictStepATTLead)
countNalls = extractBetaPVal(fileNalls, 1, 2, 6, 5, 8, 10, dictNalls, "Nalls et al. 2019", dictStepATT, fileOutATT, fileInconsistentATT, dictStepATTLead)
countKim = extractBetaPVal(fileKim, 3, 4, 6, 5, 12, 14, dictKim, "Kim et al. 2024", dictStepATT, fileOutATT, fileInconsistentATT, dictStepATTLead)
         
fileOutATT.close()
fileInconsistentATT.close()
#printing number of positions in which the other and effect were inverted (and therefore the beta value was inverted manually)
print(f"Number of positions that betas were inverted for LA-GWAS with ATT lead SNPs, Loesch: {countLARGE}, and for non lead SNPs, Kim: {countKim}, Nalls: {countNalls}")