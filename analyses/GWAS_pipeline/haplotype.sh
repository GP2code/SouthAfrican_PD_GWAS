# Copy chr 14 to haplotype folder in ./PD_GWAS_SouthAfrica/haplotype
cp ./PD_GWAS_SouthAfrica/imputation/post_imputation/related/imputed_files/chr14.* ./PD_GWAS_SouthAfrica/haplotype

cp ./PD_GWAS_SouthAfrica/imputation/post_imputation/related/imputed_files/chr8.* ./PD_GWAS_SouthAfrica/haplotype

# extract the variant information from the inputed .dose.vcf.gz file for chr 14
zgrep -E '#CHROM|rs17098735' chr14.dose.vcf.gz > chr14_snp.tsv

zgrep -E '#CHROM|rs116582124' chr8.dose.vcf.gz > chr8_snp.tsv

# Need to look at the msp for the window the variant falls into
# Copy the information into excel sheets

=IF(VLOOKUP(C14,[chr14_msp.xlsx]Sheet2!$B:$C,2,FALSE)=0,"AFR",IF(VLOOKUP(C14,[chr14_msp.xlsx]Sheet2!$B:$C,2,FALSE)=1,"EUR",IF(VLOOKUP(C14,[chr14_msp.xlsx]Sheet2!$B:$C,2,FALSE)=2,"MALAY",IF(VLOOKUP(C14,[chr14_msp.xlsx]Sheet2!$B:$C,2,FALSE)=3,"NAMA",IF(VLOOKUP(C14,[chr14_msp.xlsx]Sheet2!$B:$C,2,FALSE)=4,"SAS","")))))

=IF(VLOOKUP(C898,Sheet2!$A:$B,2,FALSE)=0,"AFR",IF(VLOOKUP(C898,Sheet2!$A:$B,2,FALSE)=1,"EUR",IF(VLOOKUP(C898,Sheet2!$A:$B,2,FALSE)=2,"MALAY",IF(VLOOKUP(C898,Sheet2!$A:$B,2,FALSE)=3,"NAMA",IF(VLOOKUP(C898,Sheet2!$A:$B,2,FALSE)=4,"SAS","")))))

# Example with rs562618836 on chr 10:109299759
# Copy file to working directory
cp ./PD_GWAS_SouthAfrica/imputation/post_imputation/related/imputed_files/chr10.dose.vcf.gz ./PD_GWAS_SouthAfrica/haplotype

# extract the variant information from the inputed .dose.vcf.gz file for chr 10
zgrep -E '#CHROM|rs565473607' chr10.dose.vcf.gz > chr10_snp3.tsv

zgrep -E '#CHROM|rs190769295' chr2.dose.vcf.gz > chr2_snp.tsv

zgrep -E '#CHROM|rs8046041' chr16.dose.vcf.gz > chr16_snp.tsv

zgrep -E '#CHROM|rs75842871' chr17.dose.vcf.gz > chr17_snp1.tsv

zgrep -E '#CHROM|rs140681053' chr17.dose.vcf.gz > chr17_snp2.tsv

zgrep -E '#CHROM|rs6490582' chr13.dose.vcf.gz > chr13_snp.tsv

zgrep -E '#CHROM|rs150739158' chr6.dose.vcf.gz > chr6_snp.tsv


# Download msp file for chr10 to local device so you can use excel
# Open msp file in text editor and then copy and paste to excel
# Look for row that have your snp in the window. spos is start position and epos is end position. 
# Copy row 1 and 2 from excel file and paste as transposed in new sheet. 
# Add the ancestry codes to top of sheet (AFR=0	EUR=1	MALAY=2	NAMA=3	SAS=4) 
# In following columns add the headers in line 10:
# c left name 
# d left genotype
# e right name
# f right genotype
# g has variant
# h ancestry left
# i ancestry right

# For column C use formula: =A11&".0"
# For column E use formula: =A11&".1"
# For column D use formula: =VALUE(LEFT(B11,1))
# For column F use formula: =VALUE(RIGHT(LEFT(B11,3),1))
# Then highlight row 10 and add filter to headers 
# Need to make sure column D and F are numbers 
# For column G use formular: =IF(OR(D11=1,F11=1),"Yes","No")

# Copy row 2 and the row with the correct window from the msp file and paste transposed into new sheet. 
# Filter column D for has variant and use VLOOKUP to put in the ancestries
# Filter then for column F for has variant and repeat VLOOKUP. 
# Add VLOOKUP to column H and I to get the ancestries from the msp file
# You can then filter column G for "Yes" and use COUNTIF to get the totals

 =IF(VLOOKUP(C331,Sheet2!$A:$B,2,FALSE)=0,"AFR",IF(VLOOKUP(C331,Sheet2!$A:$B,2,FALSE)=1,"EUR",IF(VLOOKUP(C331,Sheet2!$A:$B,2,FALSE)=2,"MALAY",IF(VLOOKUP(C331,Sheet2!$A:$B,2,FALSE)=3,"NAMA",IF(VLOOKUP(C331,Sheet2!$A:$B,2,FALSE)=4,"SAS","")))))
 
 
=IF(VLOOKUP(E216,Sheet2!$A:$B,2,FALSE)=0,"AFR",IF(VLOOKUP(E216,Sheet2!$A:$B,2,FALSE)=1,"EUR",IF(VLOOKUP(E216,Sheet2!$A:$B,2,FALSE)=2,"MALAY",IF(VLOOKUP(E216,Sheet2!$A:$B,2,FALSE)=3,"NAMA",IF(VLOOKUP(E216,Sheet2!$A:$B,2,FALSE)=4,"SAS","")))))
# Use =COUNTIF to get the number of occurrences of each ancestry and then get proportion
