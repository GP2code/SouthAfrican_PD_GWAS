setwd("./PD_GWAS_SouthAfrica/ROH")
library(ggplot2)

# Define the ancestries to loop through
ancestries <- c("LARGE-PD", "SAPDSC")

# Initialize an empty data frame to store all the results
data_all <- data.frame(NULL)

# Loop through each ancestry to process the ROH data
for (ANC in ancestries) {
  
  # Read the ROH file for the current ancestry
  data <- read.table(paste0("output", ANC, "/", ANC, "_ROH.hom.indiv"), header = TRUE)
  
  # Select the columns of interest: Individual ID, Number of Segments (NSEG), and Length in KB (KB)
  data_summary <- data[c("IID", "NSEG", "KB")]
  
  # Convert KB (length in kilobases) to megabases (MB)
  data_summary$KB <- data_summary$KB / 1000
  
  # Rename the columns for clarity
  colnames(data_summary) <- c("ID", "NROH", "SROH")  # NROH = Number of ROH, SROH = Size of ROH in MB
  
  # Add a column for the ancestry information
  data_summary$ANC <- ANC
  
  # Append the current ancestry's data to the overall data frame
  data_all <- rbind(data_all, data_summary)
}

# Plot size of ROH segments (SROH) > 1.5 Mb against the mean number of ROH segments (NROH)
ggplot(data_all, aes(x = SROH, y = NROH, colour = ANC)) +
  geom_jitter(width = 0.2, height = 0.2) +  # Adds jitter to avoid overplotting
  theme_minimal() +
  labs(colour = "Ancestry", # Label for the legend (colour by ancestry)
       x = "SROH (in Mb)",  # X-axis label for size of ROH segments > 1.5 Mb
       y = "Mean NROH",  # Y-axis label for number of ROH segments
       title = "Comparison of SROH against mean NROH") + # Set title
  scale_colour_brewer(palette = "Set1") +  # Use the "Set1" colour palette for distinct colours
  theme(legend.position = "right")  # Move legend to the right for better visibility
ggsave("./plots/SROH_1.5_vs_NROH.png",dpi=600)