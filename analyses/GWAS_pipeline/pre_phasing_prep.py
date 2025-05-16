import os
# Define input binary files
dataIn = "merge_1000G_MALAY_NAMA"
# Convert binary files to VCF format with PLINK
os.system(f"plink --bfile ./{dataIn} --recode vcf-iid bgz --out ./reference --output-chr chrMT")
# Index the resulting VCF file with BCFtools
os.system(f"bcftools index ./reference.vcf.gz")
# Define the path to the reference FASTA file
fasta = f"./PD_GWAS_SouthAfrica/reference_files/Homo_sapiens_assembly38.fasta"
# Fix the reference alignment and create a new VCF file
os.system(f"bcftools +fixref ./reference.vcf.gz -Oz -o reference_1.vcf.gz -- -f {fasta} -m flip -d")
# Re-index the normalized VCF file
os.system(f"bcftools index ./reference_1.vcf.gz -f")
# Split the normalized VCF file by chromosome and output each as a separate file
for i in range(1, 23):
 os.system(f"bcftools view --regions chr{i} -Oz -o reference_chr{i}.vcf.gz ./reference_1.vcf.gz")
