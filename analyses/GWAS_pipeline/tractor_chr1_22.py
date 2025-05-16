import pandas
import os
for i in range(1, 23):
    for method in ["TRACTOR"]:
        for PC in ["Stepwise"]:
            pgenFolder = f"./PD_GWAS_SouthAfrica/admix-kit/PGEN_PC"
            infoFolder = f"./PD_GWAS_SouthAfrica/covars"
            outFolder = f"./PD_GWAS_SouthAfrica/admix-kit/AK_results"

            os.system(f"admix assoc --pfile {pgenFolder}/stellen_{i} --family binary --pheno {infoFolder}/Covars_{PC}_recode.txt --method {method} --quantile-normalize True --out {outFolder}/stellen_{PC}_{method}_{i}")
