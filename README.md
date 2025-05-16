# Genome-wide association analyses reveal susceptibility variants linked to Parkinson's disease in the South African population using inferred global and local ancestry

`GP2 ❤️ Open Science 😍`

[![DOI](https://zenodo.org/badge/984313052.svg)](https://doi.org/10.5281/zenodo.15442537)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Last Updated:** May 2025

## **Summary**

This repository contains the code, data workflows, and results associated with the manuscript titled “**Genome-wide association analyses reveal susceptibility variants linked to Parkinson's disease in the South African population using inferred global and local ancestry**”. 

This study looks at ancestry inference, a conventional genome-wide association study (GWAS), a GWAS incorporating local ancestry, and a runs of homozygosity analysis (ROH).

**Summary statistics generated from this study will be available for use on the GP2 approved analysis platforms (Terra/Verily Workbench) on GP2 Tier 1 starting release 10 (scheduled for July 1st 2025).**

## **Data Statement**

* The data was obtained from the Global Parkinson’s Genetics Program (GP2) release 7 (GP2 release 7, DOI: 10.5281/zenodo.10962119) and access can be requested through the Accelerating Medicines Partnership in Parkinson’s Disease (AMP-PD) via the online application process ([https://www.amp-pd.org/](https://www.amp-pd.org/)).
* Requests to access the Nama datasets should be directed to Prof. Marlo Moller ([marlom@sun.ac.za](mailto:marlom@sun.ac.za)). 
* The full summary statistics from LARGE-PD are available upon request. Additionally, LARGE-PD summary statistics can be found in the Neurodegenerative Disease Knowledge Portal ([https://ndkp.hugeamp.org/phenotype.html?phenotype=Parkinsons](https://ndkp.hugeamp.org/phenotype.html?phenotype=Parkinsons)).
    * The data analyzed in this study is subject to the following licenses/restrictions: No new genetic data was generated for this study however, summary statistics for the quality and accuracy assessment of the genetic data for the NAMA participants will be made available to researchers who meet the criteria for access after application to the Health Research Ethics Committee of Stellenbosch University. 
---

## **Citation**

If you use this repository or find it helpful for your research, please cite the corresponding manuscript:

> **Genome-wide association analyses reveal susceptibility variants linked to Parkinson's disease in the South African population using inferred global and local ancestry**
> *[Step K, Peixoto Leal T, Waldo E, Madula L, Swart Y, Hernández CF, Kim JJ, Bandres-Ciga S, Global Parkinson's Genetics Program (GP2), Mata IF, Bardien S, 2025]* (DOI: )
---

## **Workflow Diagram**

This figure outlines the key steps in the analysis pipeline, including data preprocessing, ancestry inference, association analysis approaches, and runs of homozygosity overview.

![workflow](https://github.com/GP2code/SouthAfrican_PD_GWAS/blob/main/figures/GWAS_Pipeline.png)


---

## **Repository Orientation**

- The `analyses/` directory includes all the analyses discussed in the manuscript or links to the respective GitHubs where analyses were followed.
- The `figures/` directory includes all the figures and supplementary figures referenced in the manuscript (*pending publication*).
- The `tables/` directory includes all the tables and supplementary tables referenced in the manuscript (*pending publication*).
- The `data/` directory includes the output from TRACTOR (the local ancestry GWAS approach)
    - This will be available on GDrive due to size restrictions 

```
├── analyses
│   ├── GWAS_manuscript.sh
│   └── GWAS_pipeline
│       ├── add_age_covar.py
│       ├── att_chr1_22.py
│       ├── chr_ammend_snp_info.py
│       ├── convert.py
│       ├── download_1000G.py
│       ├── EffectDirections_combined.py
│       ├── EffectDirectionsPlotting.R
│       ├── ER2.py
│       ├── fix_files.py
│       ├── geneList.txt
│       ├── haplotype.sh
│       ├── MergeTRACTOR.py
│       ├── newMSP.py
│       ├── overlap.py
│       ├── post_plink_ROH_mapping.pl
│       ├── pre_phasing_prep.py
│       ├── preparation_no_ref_pca.py
│       ├── prepare_plots_with_ref.py
│       ├── remove_INFO.py
│       ├── ROH_plotting.R
│       └── tractor_chr1_22.py
├── figures
│   └── GWAS_Pipeline.png
├── LICENSE
└── TRACTOR_output [available on GDrive due to size]
    ├── TRACTOR-AFR_hg38_P_Stepwise
    ├── TRACTOR-ALL_hg38_P_Stepwise
    ├── TRACTOR-EUR_hg38_P_Stepwise
    ├── TRACTOR-MALAY_hg38_P_Stepwise
    ├── TRACTOR-NAMA_hg38_P_Stepwise
    └── TRACTOR-SAS_hg38_P_Stepwise
```

---

## **Key Analyses**

1. Quality control and imputation preparation
    - The parameters used for the quality control as well as how the files were prepped for imputation is presented in this analysis.
2. Ancestry inference
    - The global and local inferred ancestry for the study participants was investigated to use in downstream analysis.
3. Association analysis using SAIGE
    - A sparse genetic relationship matrix was created and a conventional GWAS was run using projected PCs.
4. Association analysis using local ancestry
    - Local ancestry was used as a covariate (ATT model) and as separate dosage information for ancestry (Tractor) in two approaches incorporating it into a GWAS.
5. Runs of homozygosity
    - ROH regions were identified in the datasets for four ROH parameters and statistical analysis was performed to identify association with disease status. Additional, the ROH regions were intersected with known Parkinson’s disease genes and risk loci.

---

### **Figures and Supplementary Figures**

(*Pending publication*)

### **Tables and Supplementary Tables**

(*Pending publication*)

---

## **Analysis Notebooks**

*Languages: Python, Bash, and R*

| **Notebook** | **Description** |
| --- | --- |
| GWAS_Manuscript.sh | Provides the links to the respective GitHub repositories used in the analysis as well as any additional scripts for quality control, creation of reference files, ancestry inference, GWAS using SAIGE, GWAS using Admix-Kit, effect direction analysis, and runs of homozygosity analysis. |

---

## **Software**

| **Software** | **Version(s)** | **Resource URL** | **RRID** | **Notes** |
| --- | --- | --- | --- | --- |
| Python Programming Language | 3.7.0 | [http://www.python.org/](http://www.python.org/) | RRID:SCR_008394 | *Used for quality control, imputation preparation, ancestry inference* |
| PLINK | 1.9 and 2.0 | [http://www.nitrc.org/projects/plink](http://www.nitrc.org/projects/plink) | RRID:SCR_001757 | *Used for quality control, ancestry inference, and runs of homozygosity analysis.* |
| R Project for Statistical Computing | 4.2.0 | [http://www.r-project.org/](http://www.r-project.org/) | RRID:SCR_001905 | *Used for plotting throughout the analysis, to create the genetic relationship matrix, to run the conventional association analysis,* |
| SAIGE | 1.0.0 | [https://saigegit.github.io/SAIGE-doc/](https://saigegit.github.io/SAIGE-doc/) | - | *Used to run the conventional association analysis.* |
| Gnomix | 1.0 | [https://github.com/AI-sandbox/gnomix](https://github.com/AI-sandbox/gnomix) | - | *Used to infer local and global ancestries.* |
| ADMIXTURE | 1.3.0 | [https://dalexander.github.io/admixture/](https://dalexander.github.io/admixture/) | RRID:SCR_001263 | *Used to investigate population substructure.* |
| Admix-Kit | 0.1.3 | [https://github.com/KangchengHou/admix-kit](https://github.com/KangchengHou/admix-kit) | - | *Used to run the ATT and Tractor models for the local ancestry GWAS.* |
