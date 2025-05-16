import pandas as pd

def add_age_to_psam(psam_file, covariates_file, output_file):
    # Load the psam file
    psam = pd.read_csv(psam_file, delim_whitespace=True)
    psam.rename(columns={"#IID": "IID"}, inplace=True)  # Ensure consistent column names

    # Load the covariates file
    covariates = pd.read_csv(covariates_file, sep="\t")
    covariates.rename(columns={"GP2sampleID": "IID"}, inplace=True)  # Ensure matching column names

    # Merge psam with covariates on IID
    merged = psam.merge(covariates[["IID", "AGE"]], on="IID", how="left")

    # Report missing age values
    missing_ages = merged["AGE"].isna().sum()
    print(f"Missing AGE values for {missing_ages} individuals.")

    # Save the updated psam file (keep original formatting)
    merged.to_csv(output_file, sep="\t", index=False)
    print(f"Updated .psam file saved as {output_file}!")

# Example usage
psam_file = "nba_qced_QCed.psam"
covariates_file = "covariates.tsv"
output_file = "nba_qced_QCed_with_age.psam"

add_age_to_psam(psam_file, covariates_file, output_file)

