# Two QCs were run. One including related individuals and one that will exclude the related individuals. 

########################################################################################
# QC (https://github.com/MataLabCCF/GWASQC)
########################################################################################
# Relatedness uses NAToRA (https://github.com/ldgh/NAToRA_Public) so you need the NAToRA_Pulic.py file from the GitHub page. Also need the NetworkX software (https://networkx.org) 

# Create QC steps file with the thresholds (nano qc_steps.txt). Needs to contain the thresholds you want as follows:
sex-check
ATCG
mind	0.05
geno	0.05
duplicate
heterozygosity	3
HWE control	1e-6
HWE case	1e-10
relationship	0.0884

# Create info file with the ID, sex, and status
awk 'BEGIN{OFS=","; print "ID","SEX","STATUS"} {sex=($5==1?"M":($5==2?"F":"Unknown")); status=($6==2?"Case":($6==1?"Control":"Unknown")); print $2,sex,status}' flipped_GP2_merge_STELLENBOS.fam > info_file.csv
sed 's/,/\t/g' info_file.csv > info_file.tsv

# To run main.py script
python GWASQC/main.py -i flipped_GP2_merge_STELLENBOS -o nba_qced -O qced_data -I info_file.tsv -s qc_steps.txt --NAToRA NAToRA_Public/NAToRA_Public.py

# qced_data has the logs. QC.log will have a summary with the chr1:XXX being the number of remaining variants. Final QC data is in FinalData folder.
# Convert final output to bfiles
plink2 --pfile nba_qced_QCed --make-bed --out nba_qced

# You can plot the relatedness using (plotted on local computer):
# nba_qced_input.txt is input needed for createGML_python.py and plot with yEd (https://www.yworks.com/products/yed/download#download).
python ../NAToRA_Public/createGML_python.pl -i nba_qced_input.txt -o GML_plot
# The output GML_plot.gml file is what will be needed to plot on local computer. 

########################################################################################
# QC families included(https://github.com/MataLabCCF/GWASQC)
########################################################################################
# Relatedness uses NAToRA (https://github.com/ldgh/NAToRA_Public) so you need the NAToRA_Pulic.py file from the GitHub page. Also need the NetworkX software (https://networkx.org) 
module load chpc/BIOMODULES
module load python 
pip list | grep networkx 

# Create QC steps file with the thresholds (nano qc_steps_families_included.txt)
sex-check
ATCG
mind	0.05
geno	0.05
duplicate
heterozygosity	3
HWE control	1e-6
HWE case	1e-10

# To run main.py script
python ../GWASQC/main.py -i flipped_GP2_merge_STELLENBOS -o nba_qced_families_included -O qced_data_FI -I info_file.tsv -s qc_steps_families_included.txt --NAToRA NAToRA_Public/NAToRA_Public.py

# qced_data has the logs. QC.log will have a summary with the chr1: XXX being the number of remaining variants. Final QC data is in FinalData (./PD_GWAS_SouthAfrica/quality_control_fixed_files/qced_data_families_included/qced_data_FI/FinalData)

# convert to binary files
plink2 --pfile nba_qced_families_included_QCed --make-bed --out nba_qced_FI

########################################################################################
# Imputation
########################################################################################
# To resolve strand issues pre Imputation. Used files including related individuals
# Pre-processing
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta

samtools faidx Homo_sapiens_assembly38.fasta
bwa index Homo_sapiens_assembly38.fasta

# run the following: fix_files.py to prep files for imputation
# Use this output for imputation then download the imputed results
# Download imputed results (replace Password with the one received from TOPMed)
7za x zippedfile
for chr in $(seq 1 22)
do
	7za x "Password" chr_$chr.zip
done

# To plot the ER2 values from Imputed data run the ER2.py

########################################################################################
# Creation of reference files
########################################################################################
# Need to run a file QC on the Malay and Nama ancestry files. This is done using the Mata Lab GWASQC pipeline in the same manner as the NeuroBooster Array data was QCed. Make sure files are in binary format. 

# Download the 1000 Genomes reference files using download_1000G.py. 
# You may need to remove the INFO column for the VCFs to get the files to merge. Note: You do not need the X chromosome included in generating the reference files. You can remove the INFO column using remove_INFO.py. This will also put the fixed files in a new directory. 

# The Malay and Nama files need to be run through TOPMed for imputation. Run as imputation method above. Once imputation is complete, you can proceed with creating the reference files. 

# To make reference files, you need to get the intersecting SNPs for the reference files and the NBA data to submit for Phasing on TOPMED. For this you need to merge all the reference files together to create one big reference file. 

# To merge Malay, 1000G_AFR, and recoded_nama_grch38
# For intersecting SNPs
# Step 1: Extract SNP IDs from each dataset and sort them
awk '{print $3}' recoded_1000G_AFR_EUR_SAS_EAS.pvar | sort > recoded_1000G_AFR_EUR_SAS_EAS_sorted.txt
awk '{print $3}' recoded_malay.pvar | sort > malay_sorted.txt
awk '{print $3}' recoded_nama_grch38.pvar | sort > nama_sorted.txt

# Step 2: Find and extract the common SNP IDs between all three datasets
comm -12 recoded_1000G_AFR_EUR_SAS_EAS_sorted.txt malay_sorted.txt | comm -12 - nama_sorted.txt > intersecting_snps.txt

# Step 3: Check the number of intersecting SNPs
wc -l intersecting_snps.txt
# The number of common SNPs will be displayed (942325)

# Step 4: Extract intersecting SNPs from each dataset
plink2 --pfile recoded_malay --extract intersecting_snps.txt --make-bed --out malay_intersect_snps
plink2 --pfile recoded_1000G_AFR_EUR_SAS_EAS --extract intersecting_snps.txt --make-bed --out recoded_1000G_AFR_EUR_SAS_EAS_intersect_snps
plink2 --pfile recoded_nama_grch38 --extract intersecting_snps.txt --make-bed --out recoded_nama_grch38_intersect_snps

# Step 5: Merge the files
# First, merge Malay and 1000G_AFR_EUR_SAS_EAS
plink --bfile malay_intersect_snps --bmerge recoded_1000G_AFR_EUR_SAS_EAS_intersect_snps --make-bed --out merge_1000G_MALAY

# Then, merge the result with recoded_nama_grch38
plink --bfile merge_1000G_MALAY --bmerge recoded_nama_grch38_intersect_snps --make-bed --out merge_1000G_MALAY_NAMA

# Run pre_phasing_prep.py 
# Upload to TOPMED imputation server. No rsq filter needed as not imputing. 

# To download TOPMed results use the following, replacing "Password" with the password received from TOPMed. 
for chr in $(seq 1 22)
do
	7za x "Password" chr_$chr.zip
done

########################################################################################
# Project PCA and model
########################################################################################
# Run using Mata Lab scripts (https://github.com/MataLabCCF/ProjectedPCAAndModelSelection)
# Need to run this twice. Once with -R and once without. The projected PCA uses the referene file and this will be the input for ancestry inference. The non-projected PCA is used in the GWAS. 

python /ProjectedPCAAndModelSelection/covarProjectedPCA.py \
 -A ./PD_GWAS_SouthAfrica/projected_pca/autosomes_nba_qced \
 -t ./PD_GWAS_SouthAfrica/projected_pca/covariates.tsv -R ./PD_GWAS_SouthAfrica/projected_pca/merge_1000G_MALAY_NAMA \
 -n covarProjected_including_nama_12Dec -f ./PD_GWAS_SouthAfrica/projected_pca/updated_pipeline_dec24/with_reference \
 --gcta gcta64 \
 --selectModel ./PD_GWAS_SouthAfrica/projected_pca/updated_pipeline_dec24/ProjectedPCAAndModelSelection/selectModel.R --plink1 plink --plink2 plink2
 
# recommended covariates are written to covarProjected_including_nama_12Dec_variables.tsv 

# Repeat this while not including the -R
 python /ProjectedPCAAndModelSelection/covarProjectedPCA.py \
 -A ./PD_GWAS_SouthAfrica/projected_pca/autosomes_nba_qced \
 -t ./PD_GWAS_SouthAfrica/projected_pca/covariates.tsv \
 -n covarProjected_13Dec_no_ref -f ./PD_GWAS_SouthAfrica/projected_pca/updated_pipeline_dec24/no_reference \
 --gcta gcta64 \
 --selectModel ./PD_GWAS_SouthAfrica/projected_pca/updated_pipeline_dec24/ProjectedPCAAndModelSelection/selectModel.R --plink1 plink --plink2 plink2

## The following works for plotting with reference, but not without.
# Plot the manhattans for the PCs. Input file is AutosomalPCA_VariantLoading.pcl (Use prepare_plots_with_ref.py and plot with R script plot_pcs.R). 

## For the plotting where there is no reference included in the PCA. Use the following:
# For preparation to plot use preparation_no_ref_pca.py

# To plot in R:
library(dplyr)
library(gtools)
library(ggplot2)
library(tidyverse)
cols =  c("1" = "black", "2" = "grey", "3" = "black", "4" = "grey", "5" = "black", "6" = "grey",
          "7" = "black", "8" = "grey", "9" = "black", "10" = "grey", "11" = "black", "12" = "grey",
          "13" = "black", "14" = "grey", "15" = "black", "16" = "grey", "17" = "black", "18" = "grey",
          "19" = "black", "20" = "grey", "21" = "black", "22" = "grey")
manhattan_plots <- list()
for (i in 1:50) {
  file_name <- paste0("PC", i, "ToMan.txt")
  data <- read.table(file_name, header = TRUE)  # Assuming tab-delimited
  manhattan_plots[[i]] <- ggplot(data, aes(x = BP, y = PC)) + geom_point(aes(colour = factor(CHR)))+
    scale_color_manual(values =cols)+ theme_minimal() +
    theme(legend.position = "none",axis.text.x=element_blank(), axis.ticks.x=element_blank())+
    labs(x= "CHR", y = "SNP loading", title = paste("PC", i))
}
png("Loading.png", dpi = 450)
grob = grid.arrange(grobs = manhattan_plots, ncol = 5)
dev.off()

########################################################################################
# GWAS with SAIGE
########################################################################################
# Run using Mata Lab scripts (https://github.com/MataLabCCF/)
# Notes: Need to run step 1 and 2 on genotyping data (not imputed data). PCAs calculated with genotyping data. Need to recode sex and phenotypes to 0 and 1 not 1 and 2. Need to index the vcf file. Don't do QC on imputed data (whatever TOPMed didnt do, SAIGE will do). Use imputed vcf files because it has dosage information already. (website: https://saigegit.github.io/SAIGE-doc/docs/single_step1.html and GitHub: https://github.com/saigegit/SAIGE). 
## To run SAIGE
# Create sparse Genetic Relationship Matrix (GRM) (All chromosomes)
# createSpareseGRM.R is available at https://github.com/saigegit/SAIGE/tree/main/extdata
# SAIGE needs to be run in a singularity container on the CHPC

# Need to only run this on autosomes so remove the X chr
plink --bfile nba_qced_FI --autosome --make-bed --out nba_qced_autosomes

singularity shell -B /mnt:/mnt $SAIGE

Rscript createSparseGRM.R \
--plinkFile nba_qced_autosomes \
--outputPrefix=saige_step1 \
--relatednessCutoff=0.088 \
--nThreads=8 --numRandomMarkerforSparseKin 50000 

# Step 1. Need the same plink file as above (QCED genotype file), the generated PC file but the cases need to be recoded to 0 and 1 (this was done for the _toModel.tsv file). 
# Make Results directory
mkdir ./PD_GWAS_SouthAfrica/Results

# Run step1 of SAIGE in singularity shell. Change the PCs if you want a different combination ie PC1-PC10.
Rscript step1_fitNULLGLMM.R --plinkFile=nba_qced_autosomes --phenoFile=./covarProjected_including_nama_toModel.tsv --phenoCol=DISEASE --covarColList=AGE,SEX,Autosomal_WR_PC1,Autosomal_WR_PC3,Autosomal_WR_PC7,Autosomal_WR_PC10,Autosomal_WR_PC12,Autosomal_WR_PC15,Autosomal_WR_PC16,Autosomal_WR_PC19,Autosomal_WR_PC29,Autosomal_WR_PC37,Autosomal_WR_PC42,Autosomal_WR_PC45,Autosomal_WR_PC46,Autosomal_WR_PC49,Autosomal_WR_PC50 --qCovarColList=SEX --sampleIDColinphenoFile=IID --isCateVarianceRatio=TRUE --traitType=binary --outputPrefix=./Results/Model_NBA --nThreads=8 --IsOverwriteVarianceRatioFile=TRUE --useSparseGRMtoFitNULL=TRUE --sparseGRMFile saige_step1_relatednessCutoff_0.088_50000_randomMarkersUsed.sparseGRM.mtx --sparseGRMSampleIDFile saige_step1_relatednessCutoff_0.088_50000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt

#First: index the imputed data in csi.
for chr in $(seq 1 22)
do
	bcftools index -c -f chr${chr}.dose.vcf.gz
done

# Need to loop for step2. Takes a while so run in chunks. Results are written to ./Results/Results_nba_Chrom22"
for chr in $(seq 19 22); do
Rscript step2_SPAtests.R --vcfFile=./PD_GWAS_SouthAfrica/imputation/post_imputation/related/imputed_files/chr${chr}.dose.vcf.gz --vcfFileIndex=./PD_GWAS_SouthAfrica/imputation/post_imputation/related/imputed_files/chr${chr}.dose.vcf.gz.csi --vcfField=DS --SAIGEOutputFile=./PD_GWAS_SouthAfrica/Results/Result_nba_Chrom${chr} --minMAF=0 --minMAC=20 --GMMATmodelFile=./Results/Model_NBA.rda --varianceRatioFile=./Results/Model_NBA.varianceRatio.txt --sparseGRMFile=./saige_step1_relatednessCutoff_0.088_50000_randomMarkersUsed.sparseGRM.mtx --sparseGRMSampleIDFile=./saige_step1_relatednessCutoff_0.088_50000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt --LOCO=FALSE --is_Firth_beta=TRUE --pCutoffforFirth=0.1 --is_output_moreDetails=TRUE --AlleleOrder=ref-first --chrom=chr${chr}
done

# Concatenate files to get results. Note: In SAIGE allele 2 is the effect allele 
for chr in {1..22}; do cat Result_nba_Chrom${chr} >> concatenated_results.txt; done

# If you want to check the local ancestry window of top hits, first run the ancestry inference below then follow the instructions in haplotype.sh

########################################################################################
# Ancestry inference 
########################################################################################
# Run using Mata Lab scripts (https://github.com/MataLabCCF/)
# Need to remove the EAS samples from the phased reference files (you may need to remove ancestries based on how admixed your population is. For the South African data we needed to remove EAS before proceeding).
for chr in {1..22}
do input_file="./PD_GWAS_SouthAfrica/ancestry_inference/create_reference_files/topmed_phased_files/chr${chr}.phased.vcf.gz"
    output_file="chr${chr}.filtered.phased.vcf.gz"
    
    # Remove samples listed in EAS.txt
    bcftools view --samples-file ^EAS.txt -Oz -o $output_file $input_file
    
    echo "Processed $input_file -> $output_file"
done

# For the input file. Need to extract the INFO/TYPED=1 from the imputed files. Dont need to remove related individuals for this. 
for chr in {1..22}
do
bcftools view -i "INFO/TYPED=1" ./PD_GWAS_SouthAfrica/imputation/post_imputation/related/imputed_files/chr${chr}.dose.vcf.gz -Oz -o CHR${chr}_OnlyTyped.vcf.gz
done

# Run prepareDataToLA.py. If you want the full list of input flags use
python prepareDataToLA.py -h

# Run as follows but will need to loop through each chromosome
for chr in {1..22}; do
    python prepareDataToLA.py \
    -i input_files_typed1/CHR${chr}_OnlyTyped.vcf.gz \
    -o LA_GWAS -Y python -O ./LAToEpistasis -t 10 \
    -c CorrespondenceFiles/sample_with_population.txt  \
    -l reference_phased/chr${chr}.filtered.phased.vcf.gz \
    -p reference_phased/chr${chr}.filtered.phased.vcf.gz \
    -b ${chr} -e ${chr} \
    -g ./PD_GWAS_SouthAfrica/reference_files/geneticMap/geneticMapShapeit/chr${chr}.b38.gmap \
    -f ./PD_GWAS_SouthAfrica/reference_files/Homo_sapiens_assembly38.fasta -n 10
done

# Run Gnomix for local ancestry
# False means not rephasing. gnomix.py is downloaded from https://github.com/MataLabCCF/LA_LARGE-PD. Need to use the GeneticMap files as the phasing for gnomix differs from the SHAPEIT phasing. This was run using a singularity shell on the HPC. 
for chr in {1..22}; do
singularity run -B /mnt:/mnt $GNOMIX LAToEpistasis/Phased/TargetInCommon_chr${chr}.vcf.gz \
STELLEN_${chr}_NR \
${chr} \
False \
./PD_GWAS_SouthAfrica/reference_files/geneticMap/geneticMapShapeit/G_chr${chr}.b38.gmap \
LAToEpistasis/Phased/RefInCommon_chr${chr}.vcf.gz \
CorrespondenceFiles/sample_with_population.txt config.yaml; done

# Edit number of nodes on configBest.yaml or when loading session make sure you've requested enough. The inference needs to indicate best in the configBest.yaml file. 

# To plot gnomix results
 python $pip list | grep click click 7.1.2 rich-click 1.6.1 [globus ~ 09:21:45 ]$python -c 'import click'
 
########################################################################################
# Local ancestry GWAS
########################################################################################
# Run using Mata Lab scripts (https://github.com/MataLabCCF/)

# To run LA-GWAS using admix-kit:
# Using imputed files but need to convert files and include the --set-missing-var-ids flag or code will not run because of missing rsIDs from imputation and include MAF of 0.005 to remove monomorphic variants. Input files are your imputed VCFs per chromosome. The msp files that were the output from gnomix. The NAToRA file that lists the related individuals. Then run convert.py and newMSP.py. The msp files generated from gnomix include related individuals so you need to edit the msp files so that it contains the same number of the pgen files (after removing related in the previous step)

# To run TRACTOR: First need a covar file that has IID STATUS AGE SEX PCS. The STATUS needs to be coded as 0/1 and not 1/2. the psam file here is the post QC file not including related. The eigenvec file is created during SAIGE processing. BuildCovar file: Add sex to psam file using add_age_covar.py

# Tractor is run through admix-kit and in a singularity shell on CHPC. Input file is the imputed files. The flag --remove allows you to remove related individuals or individuals according to the PCA. 

# Run the appended version for chr_ammend_snp_info.py, run tractor using tractor_chr1_22.py, and run ATT using att_chr1_22.py 
# The you going to need to merge tractor results using the following: MergeTRACTOR.py

########################################################################################
# Effect directions analysis 
########################################################################################
# To check the effect directions from previous studies, use the EffectDirections_combined.py script and for plotting use EffectDirectionsPlotting.R

########################################################################################
# Runs of homozygosity analysis 
########################################################################################
# Using the imputed files 
for chr in {1..22}; do
  # Filter by MAF >= 0.0005 and only biallelic SNPs
  bcftools view --min-af 0.0005 -Oz -o tmp_chr${chr}.biallelic.vcf.gz ./PD_GWAS_SouthAfrica/imputation/post_imputation/related/imputed_files/chr${chr}.dose.vcf.gz

  # Convert to plink format with missing variant IDs set
  plink2 --vcf tmp_chr${chr}.biallelic.vcf.gz dosage=DS \
    --set-missing-var-ids @:#:\$r:\$a --new-id-max-allele-len 100 \
    --make-bed --out chr${chr}_plink
done

# Create a merge list
ls chr*_plink.bed | sed 's/.bed//' | grep -v "chr1_plink" > mergelist.txt

# Start with chr1_plink
plink --bfile chr1_plink --merge-list mergelist.txt --make-bed --out topmed_allchr_merged

# Remove related
plink --bfile topmed_allchr_merged --remove ../Related_individuals_toRemove.txt --make-bed --out all_chr_noRelated

# Remove LD
plink --bfile all_chr_noRelated --indep-pairwise 50 5 0.5 \
--out all_chr_noRelated_noLD --maf 0.05

plink --bfile all_chr_noRelated --make-bed \
--out all_chr_noRelated_noLD \
--extract all_chr_noRelated_noLD.prune.in

# Add covars to fam file
# For this need covariate file, but with IID, SEX, DISEASE so use edit_covariate_file.py
mv chrAll_QC_covar.fam ALL_OnlyTyped_noRelated_noLD.fam

# Remove outliers due to PCs (if applicable)
# Generate outliers with LARGE_Phase2_Outliers.py
for chr in {1..22}; do /home/waldoe/beegfs/Programs/plink2 --pfile \
./RemovedRelated_LARGEPhase1/chr${chr}_LPD_phase1_filtered \
--remove OutliersToRemove.txt --make-pfile --out \
/home/waldoe/beegfs/Projects/ROH/RemovedRelated_LARGEPhase1/chr${chr}_LPD_phase1_filtered_noOutliers; done

# To calculate ROH 
plink --bfile all_chr_noRelated_noLD  \
        --homozyg group --homozyg-density 50 --homozyg-gap 1000 --homozyg-kb 1500 --homozyg-snp 100 \
        --homozyg-window-het 1 --homozyg-window-missing 5 --homozyg-window-snp 50  --homozyg-window-threshold 0.05 \
        --out stellen_ROH
        
# To get the total number of ROHs
awk 'NR > 1 {sum += $4} END {print "Total number of ROH in the dataset:", sum}' stellen_ROH.hom.indiv

# Need to summarise ROH for four ROH parameters
## To calculate FROH, run the following in R
%%R
library(data.table)
library(magrittr)
library(stringr)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)
library(RColorBrewer)

# Read the individual-level ROH data file
data <- fread("stellen_ROH.hom.indiv", header = TRUE)

# Calculate FROH: total length of ROH segments (KB) divided by the length of the human autosomal genome (GRCh38)
data$FROH <- round(data$KB / 2886428, digits = 5)  # Autosomal genome length: 2,886,428 KB

# Write the updated data with FROH to a new file
fwrite(data, "stellen_ROH_withFROH.hom.indiv")


## Summarise the four ROH parameters using R:
# Load individual ROH data with FROH from the stellen file
data <- fread("stellen_ROH_withFROH.hom.indiv", header = TRUE)

# Recode Cases and Controls: PHE == -9 is missing data, so we filter those out
dat <- data[data$PHE != -9, ]
dat$PHE <- dat$PHE - 1  # Recode phenotype: 1 -> control (0), 2 -> case (1)

# Summarize data by phenotype (case/control) using group_by and summarise from dplyr
data_tidy <- dat %>%
    group_by(PHE) %>%
    summarise(
        NSEG_mean = mean(NSEG, na.rm = TRUE),
        KB_mean = mean(KB, na.rm = TRUE),
        KBAVG_mean = mean(KBAVG, na.rm = TRUE),
        FROH_mean = mean(FROH, na.rm = TRUE),
        NSEG_se = sd(NSEG, na.rm = TRUE) / sqrt(n()),
        KB_se = sd(KB, na.rm = TRUE) / sqrt(n()),
        KBAVG_se = sd(KBAVG, na.rm = TRUE) / sqrt(n()),
        FROH_se = sd(FROH, na.rm = TRUE) / sqrt(n())
    )

# Summarize data for the entire dataset
data_all <- dat %>%
    summarise(
        NSEG_mean = mean(NSEG, na.rm = TRUE),
        KB_mean = mean(KB, na.rm = TRUE),
        KBAVG_mean = mean(KBAVG, na.rm = TRUE),
        FROH_mean = mean(FROH, na.rm = TRUE),
        NSEG_se = sd(NSEG, na.rm = TRUE) / sqrt(n()),
        KB_se = sd(KB, na.rm = TRUE) / sqrt(n()),
        KBAVG_se = sd(KBAVG, na.rm = TRUE) / sqrt(n()),
        FROH_se = sd(FROH, na.rm = TRUE) / sqrt(n())
    )

# Perform t-tests for the four parameters
model1_imp <- t.test(NSEG ~ PHE, data = dat)    # NROH
model2_imp <- t.test(KB ~ PHE, data = dat)      # SROH
model3_imp <- t.test(KBAVG ~ PHE, data = dat)   # AVROH
model4_imp <- t.test(FROH ~ PHE, data = dat)    # FROH

# Calculate mean, SE, and p-value for each parameter
# SROH (Sum of ROH lengths)
SROH_cases <- paste(round(data_tidy$KB_mean[2], 1), "±", round(data_tidy$KB_se[2], 1))
SROH_controls <- paste(round(data_tidy$KB_mean[1], 1), "±", round(data_tidy$KB_se[1], 1))
SROH_all <- paste(round(data_all$KB_mean, 1), "±", round(data_all$KB_se, 1))
SROH_pvalue <- signif(model2_imp$p.value, 3)

# AVROH (Average ROH length)
AVROH_cases <- paste(round(data_tidy$KBAVG_mean[2], 1), "±", round(data_tidy$KBAVG_se[2], 1))
AVROH_controls <- paste(round(data_tidy$KBAVG_mean[1], 1), "±", round(data_tidy$KBAVG_se[1], 1))
AVROH_all <- paste(round(data_all$KBAVG_mean, 1), "±", round(data_all$KBAVG_se, 1))
AVROH_pvalue <- signif(model3_imp$p.value, 3)

# NROH (Number of ROH segments)
NROH_cases <- paste(round(data_tidy$NSEG_mean[2], 1), "±", round(data_tidy$NSEG_se[2], 1))
NROH_controls <- paste(round(data_tidy$NSEG_mean[1], 1), "±", round(data_tidy$NSEG_se[1], 1))
NROH_all <- paste(round(data_all$NSEG_mean, 1), "±", round(data_all$NSEG_se, 1))
NROH_pvalue <- signif(model1_imp$p.value, 3)

# FROH (Genomic inbreeding coefficient)
FROH_cases <- paste(round(data_tidy$FROH_mean[2], 4), "±", round(data_tidy$FROH_se[2], 4))
FROH_controls <- paste(round(data_tidy$FROH_mean[1], 4), "±", round(data_tidy$FROH_se[1], 4))
FROH_all <- paste(round(data_all$FROH_mean, 4), "±", round(data_all$FROH_se, 4))
FROH_pvalue <- signif(model4_imp$p.value, 3)

# Consolidate results
results <- data.frame(
    Metric = c("SROH", "AVROH", "NROH", "FROH"),
    Cases = c(SROH_cases, AVROH_cases, NROH_cases, FROH_cases),
    Controls = c(SROH_controls, AVROH_controls, NROH_controls, FROH_controls),
    All = c(SROH_all, AVROH_all, NROH_all, FROH_all),
    Pvalue = c(SROH_pvalue, AVROH_pvalue, NROH_pvalue, FROH_pvalue)
)

# Save the results to a text file
fwrite(results, "Pvalues_T-test_stellen_cases_controls_all_SE.txt")

# Display the results data frame
print(results)

### Intersect PD genes with ROH
# Need the GeneList_updated.txt and post_plink_ROH_mapping.pl files
# Prep GeneList_updated.txt file
cat GeneList_updated.txt | cut -f1-2 | sed 1d > geneList.txt

cat stellen_ROH.hom.overlap  | grep 'CON' | awk -F' ' '{if($4 != "-9")print}' > stellen_ROH.hom.overlap.CON.txt
cat stellen_ROH.hom.overlap  | grep -v 'CON' | grep -v 'UNION' | awk -F' ' '{if($4 != "-9")print}' | grep -v 'NA' > stellen_ROH.hom.overlap.clean.txt

perl ../post_plink_ROH_mapping.pl stellen_ROH.hom.overlap.clean.txt stellen_ROH.hom.overlap.CON.txt  geneList.txt  stellen_ROH.hom.overlap.CON.PDgenes.txt

# The run the following in R: 
library(data.table)
library(plyr)

ROH <- fread("stellen_ROH.hom.overlap.CON.PDgenes.txt", header = T)
ROH$case <- ldply(strsplit(as.character(ROH$case_to_control), split = ":"))[[1]]
ROH$control <- ldply(strsplit(as.character(ROH$case_to_control), split = ":"))[[2]]
ROH$NcasesWithROHs <- as.numeric(ROH$case)
ROH$NcontrolsWithROHs <- as.numeric(ROH$control)
ROH$caseN <- as.numeric(ROH$case)
ROH$controlN <- as.numeric(ROH$control)
ROH$combinedN <- ROH$caseN + ROH$controlN
ROH$P <- NA
for(i in 1:length(ROH$case))
{
  thisP <- prop.test(x = c(ROH$NcasesWithROHs[i], ROH$NcontrolsWithROHs[i]), n = c(661,  737)) # note the last part should be the total N for cases and then the total N for controls.
  ROH$P[i] <- thisP$p.value
}
ROH$total_ROH_count <- ROH$caseN + ROH$controlN

ROH$propCases <- ROH$caseN/661 # number of cases 
ROH$propControls <- ROH$controlN/737  #number of controls
ROH$caseEnriched <- ifelse(ROH$propCases > ROH$propControls, 1, 0)
ROH_subsetted <- subset(ROH, total_ROH_count >= 1 & caseEnriched == 1)
Ngenes <- length(ROH_subsetted$PD)
ROH_subsetted$passMultiTest <- ifelse(ROH_subsetted$P <= (0.05/Ngenes),1,0)

str(ROH_subsetted)
fwrite(ROH_subsetted,"stellen_ROH_caseEnriched.txt")

# Run the python script overlap.py to get a summary of the overlapping regions and output to new file
python overlap.py > outputROHtotals.txt

# To plot the results run ROH_plotting.R





















