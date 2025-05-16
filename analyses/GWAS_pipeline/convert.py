import os

# Define folder paths
pgenFolder = "./PD_GWAS_SouthAfrica/admix-kit/PGEN_PC"  # Folder to store PGEN files
imputationFolder = "./PD_GWAS_SouthAfrica/imputation/post_imputation/related/imputed_files"  # Imputed VCF data
relationship = "./PD_GWAS_SouthAfrica/nba_qced_output_toRemove.txt"  # Individuals to remove

# Iterate over chromosomes 
for i in list(range(1, 23)):
    mspFolder = f"./PD_GWAS_SouthAfrica/STELLEN_{i}_NR/query_results.msp"  # Local ancestry results

    # Construct commands
    plink_cmd = (
        f"plink2 --vcf {imputationFolder}/chr{i}.dose.vcf.gz "
        f"--remove {relationship} --make-pgen --out {pgenFolder}/stellen_{i} "
        f"--maf 0.005 --max-alleles 2 --rm-dup exclude-all --snps-only "
        f"--set-missing-var-ids @:#:\\$r:\\$a"
    )

    newmsp_cmd = (
        f"python newMSP.py -p {pgenFolder}/stellen_{i}.psam "
        f"-m {mspFolder} -o {pgenFolder}/MSP_{i}.msp"
    )

    admix_cmd = (
        f"admix lanc-convert {pgenFolder}/stellen_{i} "
        f"--rfmix {pgenFolder}/MSP_{i}.msp --out {pgenFolder}/stellen_{i}.lanc"
    )

    # Run the commands
    os.system(plink_cmd)
    os.system(newmsp_cmd)
    os.system(admix_cmd)
