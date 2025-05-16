import os
for i in range(22, 0, -1):
    for method in ["TRACTOR"]:
        for PC in ["PC1-10", "Stepwise"]:
            pgenFolder = f"./PD_GWAS_SouthAfrica/admix-kit/PGEN_PC"
            infoFolder = f"./PD_GWAS_SouthAfrica/covars"
            outFolder = f"./PD_GWAS_SouthAfrica/admix-kit/AK_results"

            os.system(f"admix append-snp-info --pfile {pgenFolder}/stellen_{i} --out {outFolder}/stellen_{PC}_{method}_{i}")
