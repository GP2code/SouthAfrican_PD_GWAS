# Run the following as a python script
import pandas as pd

# Load the files using pandas
roh_overlap_pd_genes = pd.read_csv('stellen_ROH.hom.overlap.CON.PDgenes.txt', sep='\t')
roh_case_enriched = pd.read_csv('stellen_ROH_caseEnriched.txt', sep=',')

# Count the number of ROH overlapping PD genes (excluding the header is automatically handled)
roh_overlap_pd_genes_count = len(roh_overlap_pd_genes)
print(f"Number of ROH overlapping PD genes: {roh_overlap_pd_genes_count}")

# Count the number of ROH enriched in Cases
roh_case_enriched_count = len(roh_case_enriched)
print(f"Number of ROH enriched in Cases: {roh_case_enriched_count}")

# Count the number of ROH that Pass Bonferroni using column name 'passMultiTest'
roh_pass_bonferroni_count = roh_case_enriched[roh_case_enriched['passMultiTest'] == 1].shape[0]
print(f"Number of ROH that Pass Bonferroni: {roh_pass_bonferroni_count}")

