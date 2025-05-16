import os
# Define input binary files
dataIn = "nba_qced_FI"
# Convert binary files to VCF format with PLINK
os.system(f"plink --bfile ./{dataIn} --recode vcf-iid bgz --out ./nba_FI --output-chr chrMT")

# Index the resulting VCF file with BCFtools
os.system(f"bcftools index ./nba_FI.vcf.gz")

# Define the path to the reference FASTA file
fasta = f"./PD_GWAS_SouthAfrica/reference_files/Homo_sapiens_assembly38.fasta"

# Fix the reference alignment and create a new VCF file
os.system(f"bcftools +fixref ./nba_FI.vcf.gz -Oz -o nba_FI_1.vcf.gz -- -f {fasta} -m flip -d")

# Re-index the normalized VCF file
os.system(f"bcftools index ./nba_FI_1.vcf.gz -f")

# Split the normalized VCF file by chromosome and output each as a separate file
for i in range(1, 23):
    os.system(f"bcftools view --regions chr{i} -Oz -o NBA_FI_chr{i}.vcf.gz ./nba_FI_1.vcf.gz")
