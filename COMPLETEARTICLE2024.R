#Before proceeding with assessing the quality of the current working data frames, 
#several preliminary steps were necessary to successfully prepare the R workplace. 

#Five key packages were installed (and subsequentially loaded into RStudio prior to beginning any analyses).
#Specifically, “httr”, “jsonlite”, “dplyr”, “readr” and “readxl” where installed and used throughout the entirety
#of this study. The Httr package makes http requests in R by providing a wrapper for the curl package.
#The jsonlite package consists of a JSON/parser generator implementing bidirectional mapping between 
#JSON data and R data. The dplyr package provides a set of grammar for working with data manipulation protocols.
#The readr package provides a way to read rectangular data, such as .csv files, into the R studio software. 
#Finally, the readxl package provides a way to read excel files into the working environment. 

# Install necessary packages if not already installed
if (!requireNamespace("httr", quietly = TRUE)) {
  install.packages("httr")
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("readr", quietly = TRUE)) {
  install.packages("readr")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("readxl", quietly = TRUE)) {
  install.packages("readxl")
}
if (!requireNamespace("reshape2", quietly = TRUE)) {
  install.packages("reshape2")
}
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}
if (!requireNamespace("sf", quietly = TRUE)) {
  install.packages("sf")
}
if (!requireNamespace("rnaturalearth", quietly = TRUE)) {
  install.packages("rnaturalearth")
}
if (!requireNamespace("rnaturalearthdata", quietly = TRUE)) {
  install.packages("rnaturalearthdata")
}
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("ropensci/rnaturalearthhires")

# Load packages
library(readr)
library(dplyr)
library(readxl)
library(httr)
library(jsonlite)
library(ggplot2)
library(reshape2)
library(ggrepel)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)

# Set working directory
setwd("/Users/alexboudreau/Desktop/Ete")

#From here, the initial data frame (present as an excel file, hence the implementation of the “readxl” package) 
#was imported into the working environment for initial analysis of data cleanliness and extent of workability.
#Being that the file consisted of a poly-tabular excel file, the spreadsheets of interest 
#(Acadian Peninsula sample data, Restigouche sample data, and Madawaska sample data) were isolated 
#and split into separate mono-tabular excel files to then have all three new data frames be imported 
#into the working environment.

# Import files
Penninsule <- read_csv("Pen.csv")
Restigouche <- read_csv("Res.csv")
Madawaska <- read_csv("Mad.csv")

# Visualize data tables to assess working state
summary(Penninsule)
summary(Restigouche)
summary(Madawaska)
head(Penninsule)
head(Restigouche)
head(Madawaska)

# Removal of non-essential columns in Penninsule dataframe
Penninsule = subset(Penninsule, select = -c(1,3,4,11))

# View column names to assess joining capacity
colnames(Penninsule)
colnames(Restigouche)
colnames(Madawaska)

# Renaming columns to assure full joint capability
Madawaska <- Madawaska %>%
  rename(Gene= `Gene (NM)`)
Restigouche <- Restigouche %>%
  rename(Gene = `Gene (NM)`)
Penninsule <- Penninsule %>%
  rename(Code=`Code GENMB`, proteine = Protein, Gene = `Gène (abbrév) Transcript`)

# Removal of empty and/or incomplete rows
Madawaska <- Madawaska[-(78:94),]
Restigouche <- Restigouche [-(53:62),]
Penninsule <- Penninsule[-c(2,10,17,22,39,42,56,62,65,66,69,76,77,82,86),]

# Joining all three dataframes into one complete dataframe
Results <- Madawaska %>%
  full_join(Restigouche) %>%
  full_join(Penninsule)

# Viewing resultant dataframe
summary(Results)
colnames(Results)

# Checking for duplicates in variable columns
table(Results$Gene)
table(Results$Base)
table(Results$proteine)

# Extract gene names and version numbers
Results$Gene_Name <- gsub("\\s*\\(NM_[0-9]+\\.[0-9]+\\)", "", Results$Gene)
Results$Version <- as.numeric(gsub(".*\\(NM_[0-9]+\\.([0-9]+)\\)", "\\1", Results$Gene))

# Find the latest version for each gene
latest_versions <- aggregate(Version ~ Gene_Name, data = Results, max)

# Merge the latest version back to the original data frame
Results_merged <- merge(Results, latest_versions, by = "Gene_Name", all.x = TRUE)

# Replace older versions with the latest version
Results_merged$Gene <- with(Results_merged, paste(Gene_Name, "(NM_000018.", Version.y, ")", sep = ""))

# Remove unnecessary columns
Results_cleaned <- Results_merged[, c("Gene")]

# Print and verify the cleaned data frame
print(Results_cleaned)
table(Results_merged$Gene)


# Manual insertion of gnomAD results is done in the df Results_Merged at this point
# The second part of this script starts from the point where these data have been integrated into the df

#Being that the goal of this study was to compare allelic frequencies of observed genetic mutations in
#certain key Acadian populations, reference profiles of respective allelic frequencies needed to be 
#outsourced to be able to proceed with statistical comparisons (Fischer’s exact testing, as 
#discussed further in this section). Specifically, genomic data from the gnomAD database was 
#obtained and filtered to integrate only data representing their “European non-Finnish” cohort.
#Before proceeding with the rest of the data manipulation in R, this data was manually incorporated 
#into the “Results” data frame. 

# Import data from CSV file
Results_Merged <- read.table("Results_Merged.csv", sep = ",", header = TRUE)

# Manually change values in allele_count_acadian column for specified rows
rows_to_change <- c(1, 2, 3, 6, 8, 17, 18, 19, 22, 23, 28, 30, 32, 34, 36, 50, 62, 73, 72, 75, 76, 79, 80, 86, 85, 87, 93, 97, 100, 101, 103, 104, 106, 108, 112, 113, 116, 118, 117, 115)
new_values <- c(1, 1, 1, 1, 2, 2, 2, 2, 2, 12, 7, 2, 3, 2, 2, 2, 2, 1, 1, 1, 1, 8, 1, 1, 1, 1, 2, 11, 2, 2, 1, 1, 3, 1, 2, 2, 1, 2, 2, 1)

# Create a vector of NA values with the length of Results_Merged$allele_count_acadian
replacement_values <- rep(NA, nrow(Results_Merged))

# Assign new values to the specified rows
replacement_values[rows_to_change] <- new_values

# Replace the values in allele_count_acadian column
Results_Merged$allele_count_acadian <- replacement_values

# Create a frequency table of genes
gene_counts <- table(Results_Merged$Gene)


#Since every candidate possesses two distinct alleles for every gene tested, gene counts were normalized 
#by using the sum() function to create a new column named Normalized_Counts – the merge() function was 
#subsequentially employed to incorporate this new column into the working data frame.
#By dividing the allele count of a gene by the total allelic count of the sample, the allelic frequency (AF) 
#was able to be obtained for every gene in the data frame. The same principle was applied to obtain the 
#heterozygote frequency (HF) for every gene in the data frame.

#Allele frequency  (Positive alleles)/(Total alleles)
#Heterozygote frequency=  ((Positive alleles-2(Homozygotes)))/((0.5 x Total alleles))

#Being that every gene in the data frame now had nonzero numeric values for the following categories 
#: AF_Acadian, HF_Acadian, AF_gnomAD, HF_gnomAD, AF and HF ratios were obtained and imported in 
#the working data frame.

#AF ratio=  (AF_Acadian)/(AF_gnomAD)
#hF ratio=  (HF_Acadian)/(HF_gnomAD)


# Normalize counts
twice_sum_of_counts <- 2 * sum(gene_counts)
normalized_counts <- gene_counts / twice_sum_of_counts
new_data <- data.frame(Gene = names(gene_counts), `AF acadian` = normalized_counts, stringsAsFactors = FALSE)

# Merge normalized counts with original data
Results_Merged <- merge(Results_Merged, new_data, by = "Gene", all.x = TRUE)

# Calculate AF ratio and HF acadian
Results_Merged$AF_ratio <- Results_Merged$AF.acadian.Freq / Results_Merged$Allele_frequency
Results_Merged$HF_acadian <- Results_Merged$AF.acadian.Freq / 119

# Duplicate AF gnomAD NFE column and rename it to HF gnomAD NFE
Results_Merged$HF_gnomAD_NFE <- Results_Merged$Allele_frequency

# Calculate HF ratio
Results_Merged$HF_ratio <- Results_Merged$HF_acadian / Results_Merged$HF_gnomAD_NFE

# Calculate allele_count_acadian
Results_Merged$allele_count_acadian <- Results_Merged$AF.acadian.Freq * 238

# Add a constant column "acadian_count"
Results_Merged$acadian_count <- 360

# Initialize p-value column
Results_Merged$pval <- NA

# Remove rows with NA values in Allele_count and Allele_number columns
Results_Merged <- Results_Merged[complete.cases(Results_Merged[, c("Allele_count", "Allele.number")]), ]

#Being that the function called to perform the Fischer’s exact two-tailed test operates under a conditional
#loop framework, a column in which the calculated p-values will be able to be stored was created before 
#running the Fischer testing script. Additionally, the complete.cases() function was called before running 
#the Fischer testing script to remove any rows having incomplete data, as this would result in either a null
#or NA p-value output from the Fischer test.
#The matrix() and nrow() functions were used to create a matrix in which a conditional loop would run through 
#every numerical value in the allele_count column and perform Fischer testing. In the event of a value being 
#of the zero or negative format, it was to be skipped via the any() function and the conditional loop would 
#continue along its designated matrix – the creation of this conditional loop and matrix before the calling 
#of the associated Fischer testing function was imperative to the success of such testing.
#The fischer.test() function was performed on every tested gene and its resultant p-value was sent to
#the previously created p-value column in the working data frame. The p.adjust() function was subsequentially
#called with the method = “bonferonni” and alpha = 0.05 conditional inserts to perform Bonferonni corrections 
#on the resultant p-values. Finally, a column called Significant_Rows was created by filtering the p-values
#based on the criteria of being inferior to 0.05. the print() function was once again called to visualize 
#the genes in this study that possessed a significant p-value. The order() function was also called to sort 
#the Gene column via it’s corresponding p-values. This complete data frame was then saved and exported as both
#.xlx and .csv type files.

#p = (((A+C)/A)((B+D)/B))/((N/(A+B)))=  (A+B)!(C+D)!(A+C)!(B+D)!/A!B!C!D!N!
#α_adjusted =  α/n  


# Perform Fisher's exact test and store p-values
for (i in 1:nrow(Results_Merged)) {
  contingency_table <- matrix(c(
    Results_Merged$allele_count_acadian[i], 
    Results_Merged$acadian_count[i],
    Results_Merged$Allele_count[i], 
    Results_Merged$`Allele.number`[i]
  ), ncol = 2, byrow = TRUE)
  if (any(contingency_table < 0)) {
    next  # Skip rows with negative values
  }
  p_val <- fisher.test(contingency_table, simulate.p.value = TRUE)$p.value
  Results_Merged$pval[i] <- p_val
}

# Apply Bonferroni correction
bonferroni_corrected_p_values <- p.adjust(Results_Merged$pval, method = "bonferroni")
Results_Merged$bonferroni_corrected_p_values <- bonferroni_corrected_p_values

# Set the significance threshold (alpha)
alpha <- 0.05

# Filter significant rows based on Bonferroni corrected p-values
significant_rows <- Results_Merged[Results_Merged$bonferroni_corrected_p_values < alpha, ]

# Print significant rows
print("Significant rows:")
print(significant_rows)

# Check for problematic rows causing errors
problematic_rows <- which(Results_Merged$allele_count_acadian < 0 | !is.finite(Results_Merged$allele_count_acadian))
print("Problematic rows:")
print(Results_Merged[problematic_rows, ])

# Check for NA or NaN values in allele_count_acadian column
na_rows <- which(is.na(Results_Merged$allele_count_acadian) | is.nan(Results_Merged$allele_count_acadian))
print("Rows with NA or NaN values in allele_count_acadian column:")
print(Results_Merged[na_rows, ])

# Check for negative values in allele_count_acadian column
negative_rows <- which(Results_Merged$allele_count_acadian < 0)
print("Rows with negative values in allele_count_acadian column:")
print(Results_Merged[negative_rows, ])

# Check data completeness and structure
str(Results_Merged)

# Check column names
colnames(Results_Merged)

# Check if necessary columns exist in the dataframe
required_columns <- c("allele_count_acadian", "acadian_count", "Allele_count", "Allele_number")
all_present <- all(required_columns %in% colnames(Results_Merged))
print("Are all required columns present?")
print(all_present)

# Sort the dataframe by p-value in ascending order
Results_Merged <- Results_Merged[order(Results_Merged$pval), ]

# Print the top most significant rows
print("Top most significant rows:")
print(head(Results_Merged))

# Sort the dataframe by Bonferroni corrected p-values in ascending order
Results_Merged <- Results_Merged[order(Results_Merged$bonferroni_corrected_p_values), ]

# Print the sorted dataframe
print("Sorted dataframe by Bonferroni corrected p-values:")
print(Results_Merged)

print(colnames(Results_Merged))

# Assuming each row is unique per subject, we'll create a subject identifier
Results_Merged$Subject <- paste0("P", seq_len(nrow(Results_Merged)))

# Check the updated column names
print(colnames(Results_Merged))

# Group by Gene and Subject, then summarize Variant_Count
summary_data <- Results_Merged %>%
  group_by(Gene, Subject) %>%
  summarize(Variant_Count = n()) %>%
  ungroup()

# Display the summarized data
print(summary_data)

# Create the plot
plot <- ggplot(summary_data, aes(x = Subject, y = Gene, fill = Variant_Count)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "darkblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        legend.position = "right") +
  labs(title = "Number of Variant Alleles per Subject", x = "Subject", y = "Gene") +
  geom_text(aes(label = Variant_Count), size = 4, color = "black")

# Display the plot
print(plot)

# Inspect the first few rows of the data
head(Results_Merged)

# Check if there are any duplicate Gene-Subject pairs
duplicate_check <- Results_Merged %>%
  group_by(Gene, Subject) %>%
  tally() %>%
  filter(n > 1)

print(duplicate_check)

# Group by Gene and Subject, then summarize Variant_Count by summing Allele_count
summary_data_sum <- Results_Merged %>%
  group_by(Gene, Subject) %>%
  summarize(Variant_Count = sum(Allele_count)) %>%
  ungroup()

# Display the summarized data
print(summary_data_sum)

# Create the plot
plot_sum <- ggplot(summary_data_sum, aes(x = Subject, y = Gene, fill = Variant_Count)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "purple") +
  labs(
    title = "Sum of Variant Alleles per Subject",
    x = "Subject",
    y = "Gene",
    fill = "Variant_Count"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.text.y = element_text(size = 6)
  )

# Print the plot
print(plot_sum)

# Aggregate data by Gene to count the number of unique subjects per gene
allele_summary <- Results_Merged %>%
  group_by(Gene) %>%
  summarize(Total_Variant_Count = n()) %>%
  ungroup()

# Display the summarized data
print(allele_summary)

# Create the heatmap
plot_allele_heatmap <- ggplot(allele_summary, aes(x = Gene, y = 1, fill = Total_Variant_Count)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "purple") +
  labs(
    title = "Total Number of Variant Alleles per Gene",
    x = "Gene",
    y = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# Print the plot
print(plot_allele_heatmap)

# Read the data from the Excel file
file_path <- "~/Desktop/DATA 2.xlsx"
df <- read_excel(file_path, sheet = 1)

# Clean the column names
colnames(df) <- c("Index", "Gene", "Base", "Protein", "Allele_count_gnomAD", "Allele_number_gnomAD",
                  "Number_of_homozygotes", "Allele_frequency_gnomAD", "Mode_of_transmission", 
                  "Zygosity", "Associated_condition", "Frequency_Acadian", "Allele_Frequency")

# Convert relevant columns to numeric, coercing any problematic values to NA
df <- df %>%
  mutate(
    Frequency_Acadian = as.numeric(Frequency_Acadian),
    Allele_count_gnomAD = as.numeric(Allele_count_gnomAD),
    Allele_number_gnomAD = as.numeric(Allele_number_gnomAD),
    Allele_Frequency = as.numeric(Allele_Frequency),
    Allele_frequency_gnomAD = as.numeric(Allele_frequency_gnomAD)
  )

# Remove rows with NA values in relevant columns
df <- df %>% filter(!is.na(Frequency_Acadian) & !is.na(Allele_count_gnomAD) & !is.na(Allele_number_gnomAD))

# Function to perform one-tailed Fisher's exact test (alternative = "greater")
fisher_test <- function(frequency_acadian, allele_count_gnomad, allele_number_gnomad) {
  positive_acadian <- as.numeric(frequency_acadian)
  negative_acadian <- 360 - as.numeric(frequency_acadian)
  positive_gnomad <- as.numeric(allele_count_gnomad)
  negative_gnomad <- as.numeric(allele_number_gnomad) - as.numeric(allele_count_gnomad)
  
  # Ensure all values are non-negative
  if (positive_acadian < 0) positive_acadian <- 0
  if (negative_acadian < 0) negative_acadian <- 0
  if (positive_gnomad < 0) positive_gnomad <- 0
  if (negative_gnomad < 0) negative_gnomad <- 0
  
  table <- matrix(c(positive_acadian, negative_acadian, positive_gnomad, negative_gnomad), nrow = 2)
  test_result <- fisher.test(table, alternative = "greater")
  return(test_result$p.value)
}

# Apply Fisher's test to each row
df$p_value <- mapply(fisher_test, df$Frequency_Acadian, df$Allele_count_gnomAD, df$Allele_number_gnomAD)

# Apply Bonferroni correction
df$bonferroni_corrected_p <- p.adjust(df$p_value, method = "bonferroni")

# Apply FDR correction
df$fdr_corrected_p <- p.adjust(df$p_value, method = "BH")

# Filter significant results for each type of p-value
significant_fisher <- df %>% filter(p_value < 0.05)
significant_bonferroni <- df %>% filter(bonferroni_corrected_p < 0.05)
significant_fdr <- df %>% filter(fdr_corrected_p < 0.05)

# Cartesian plot for Bonferroni corrected p-value without log scale
cartesian_plot_bonferroni <- ggplot(df, aes(x = Allele_Frequency, y = -log10(bonferroni_corrected_p))) +
  geom_point(color = "grey40", alpha = 0.6, size = 2) +
  geom_point(data = significant_bonferroni, aes(x = Allele_Frequency, y = -log10(bonferroni_corrected_p)), color = "red", size = 2) +
  geom_text_repel(data = significant_bonferroni, aes(x = Allele_Frequency, y = -log10(bonferroni_corrected_p), label = Gene), 
                  color = "red", size = 3, max.overlaps = 10) +
  labs(title = "Significant Genes based on Bonferroni Corrected p-value",
       x = "Allele Frequency",
       y = "-log10(Bonferroni Corrected p-value)") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())

# Print the improved Cartesian plot without log scale
print(cartesian_plot_bonferroni)

# Save the improved Cartesian plot as an image
ggsave("cartesian_plot_bonferroni_no_log.png", cartesian_plot_bonferroni)

# Filter out rows with Allele_frequency_gnomAD equal to 0 to avoid division by zero
filtered_df <- df %>%
  filter(Allele_frequency_gnomAD != 0)

# Calculate the allele frequency ratio
filtered_df <- filtered_df %>%
  mutate(allele_frequency_ratio = Allele_Frequency / Allele_frequency_gnomAD)

# Get the top 20 genes with the highest allele frequency ratios
top_genes <- filtered_df %>%
  arrange(desc(allele_frequency_ratio)) %>%
  head(20)

# Visualization of allele frequency ratios
plot_top_genes <- ggplot(top_genes, aes(x = reorder(Gene, -allele_frequency_ratio), y = allele_frequency_ratio, fill = Gene)) +
  geom_bar(stat = "identity", fill = "grey80") +
  coord_flip() +
  labs(title = "Top 20 Allele Frequency Ratios (Allele Frequency / Allele Frequency gnomAD)",
       x = "Gene",
       y = "Allele Frequency Ratio") +
  theme_minimal() +
  theme(legend.position = "none", # Remove legend
        axis.title.y = element_blank(), # Remove y-axis title
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) # Add plot margin

# Save the plot as an image
ggsave("top_20_allele_frequency_ratio_plot.png", plot_top_genes)

# Print the data with allele frequency ratios
print("Top 20 data with allele frequency ratios:")
print(top_genes)

# Print the top 20 genes plot
print(plot_top_genes)

# Calculate statistics for allele frequency ratios
allele_freq_stats <- filtered_df %>%
  summarise(
    mean_ratio = mean(allele_frequency_ratio, na.rm = TRUE),
    median_ratio = median(allele_frequency_ratio, na.rm = TRUE),
    sd_ratio = sd(allele_frequency_ratio, na.rm = TRUE),
    min_ratio = min(allele_frequency_ratio, na.rm = TRUE),
    max_ratio = max(allele_frequency_ratio, na.rm = TRUE)
  )

# Print the statistics
print("Statistics for allele frequency ratios:")
print(allele_freq_stats)

# Calculate the number of genes with a positive allele frequency ratio
positive_ratio_count <- filtered_df %>%
  filter(allele_frequency_ratio > 0) %>%
  summarise(n_distinct(Gene)) %>%
  pull()

# Print the count of genes with a positive ratio
print(paste("Number of genes with a positive allele frequency ratio:", positive_ratio_count))

# Filter the genes where the ratio is > 0 but < 1
ratio_between_0_and_1_genes <- filtered_df %>%
  filter(allele_frequency_ratio > 0 & allele_frequency_ratio < 1)

# Arrange by allele_frequency_ratio
ratio_between_0_and_1_genes <- ratio_between_0_and_1_genes %>%
  arrange(allele_frequency_ratio)

# Get the top 8 rows
ratio_between_0_and_1_genes <- head(ratio_between_0_and_1_genes, 8)

# Print the 8 genes with an allele frequency ratio > 0 but < 1
print("Genes with an allele frequency ratio > 0 but < 1 (top 8):")
print(ratio_between_0_and_1_genes)

# Calculate the number of genes where the ratio is >= 1
ratio_gte_1_count <- filtered_df %>%
  filter(allele_frequency_ratio >= 1) %>%
  summarise(n_distinct(Gene)) %>%
  pull()

# Print the count of genes where the ratio is >= 1
print(paste("Number of genes with an allele frequency ratio >= 1:", ratio_gte_1_count))

# Visualization: Histogram of allele frequency ratios
hist_plot <- ggplot(filtered_df, aes(x = allele_frequency_ratio)) +
  geom_histogram(binwidth = 0.1, fill = "grey20", color = "black") +
  scale_x_log10() + # Log scale for better visualization
  labs(title = "Histogram of Allele Frequency Ratios",
       x = "Allele Frequency Ratio (Log Scale)",
       y = "Count") +
  theme_minimal()

# Save the histogram plot
ggsave("allele_frequency_ratio_histogram.png", hist_plot)

# Prepare data for the bar plot of gene counts in different ratio ranges
counts_df <- data.frame(
  range = c("0 < ratio < 1", "ratio >= 1"),
  count = c(nrow(ratio_between_0_and_1_genes), ratio_gte_1_count)
)

# Visualization: Bar plot of gene counts in different ratio ranges
bar_plot <- ggplot(counts_df, aes(x = range, y = count, fill = range)) +
  geom_bar(stat = "identity", fill = "grey60") +
  labs(title = "Counts of Genes in Different Allele Frequency Ratio Ranges",
       x = "Allele Frequency Ratio Range",
       y = "Gene Count") +
  theme_minimal() +
  theme(legend.position = "none") # Remove legend

# Save the bar plot
ggsave("allele_frequency_ratio_bar_plot.png", bar_plot)

# Show the plots
print(hist_plot)
print(bar_plot)

# Load the map data
canada <- ne_states(country = "Canada", returnclass = "sf")

# Define coordinates for the regions to highlight (more accurate)
highlight_regions <- data.frame(
  region = c("Restigouche", "Madawaska", "Acadian Peninsula"),
  lon = c(-66.67, -68.32, -64.68),
  lat = c(47.95, 47.36, 47.75)
)

# Convert the data frame to an sf object
highlight_sf <- st_as_sf(highlight_regions, coords = c("lon", "lat"), crs = 4326)

# Create the map
ggplot(data = canada) +
  geom_sf(fill = "white", color = "black") +
  geom_sf(data = highlight_sf, fill = "black", color = "black", size = 2, shape = 21) +
  geom_text(data = highlight_regions, aes(x = lon, y = lat, label = region), 
            nudge_y = 0.1, size = 3, color = "black") +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white")
  ) +
  coord_sf(xlim = c(-69, -63), ylim = c(46, 49)) +
  ggtitle("Highlighted Areas in Atlantic Canada") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14)
  )

# Save the map
ggsave("/mnt/data/highlighted_areas_map_with_labels.png", width = 10, height = 8)





# Load necessary libraries
library(readxl)
library(dplyr)
library(ggplot2)

# Read the data from the Excel file
file_path <- "~/Desktop/DATA 2.xlsx"
df <- read_excel(file_path, sheet = 1)

# Clean the column names
colnames(df) <- c("Index", "Gene", "Base", "Protein", "Allele_count_gnomAD", "Allele_number_gnomAD",
                  "Number_of_homozygotes", "Allele_frequency_gnomAD", "Mode_of_transmission", 
                  "Zygosity", "Associated_condition", "Frequency_Acadian", "Allele_Frequency")

# Convert relevant columns to numeric, coercing any problematic values to NA
df <- df %>%
  mutate(
    Frequency_Acadian = as.numeric(Frequency_Acadian),
    Allele_count_gnomAD = as.numeric(Allele_count_gnomAD),
    Allele_number_gnomAD = as.numeric(Allele_number_gnomAD),
    Allele_Frequency = as.numeric(Allele_Frequency),
    Allele_frequency_gnomAD = as.numeric(Allele_frequency_gnomAD)
  )

# Remove rows with NA values in relevant columns
df <- df %>% filter(!is.na(Frequency_Acadian) & !is.na(Allele_count_gnomAD) & !is.na(Allele_number_gnomAD))

# Function to perform one-tailed Fisher's exact test (alternative = "greater")
fisher_test <- function(frequency_acadian, allele_count_gnomad, allele_number_gnomad) {
  positive_acadian <- as.numeric(frequency_acadian)
  negative_acadian <- 360 - as.numeric(frequency_acadian)
  positive_gnomad <- as.numeric(allele_count_gnomad)
  negative_gnomad <- as.numeric(allele_number_gnomad) - as.numeric(allele_count_gnomad)
  
  # Ensure all values are non-negative
  if (positive_acadian < 0) positive_acadian <- 0
  if (negative_acadian < 0) negative_acadian <- 0
  if (positive_gnomad < 0) positive_gnomad <- 0
  if (negative_gnomad < 0) negative_gnomad <- 0
  
  table <- matrix(c(positive_acadian, negative_acadian, positive_gnomad, negative_gnomad), nrow = 2)
  test_result <- fisher.test(table, alternative = "greater")
  return(test_result$p.value)
}

# Apply Fisher's test to each row
df$p_value <- mapply(fisher_test, df$Frequency_Acadian, df$Allele_count_gnomAD, df$Allele_number_gnomAD)

# Apply Bonferroni correction
df$bonferroni_corrected_p <- p.adjust(df$p_value, method = "bonferroni")

# Apply FDR correction
df$fdr_corrected_p <- p.adjust(df$p_value, method = "BH")

# Filter significant results for each type of p-value
significant_fisher <- df %>% filter(p_value < 0.05)
significant_bonferroni <- df %>% filter(bonferroni_corrected_p < 0.05)
significant_fdr <- df %>% filter(fdr_corrected_p < 0.05)

# Plot significant genes based on raw p-value
plot_significant_fisher <- ggplot(significant_fisher, aes(x = reorder(Gene, -p_value), y = -log10(p_value))) +
  geom_bar(stat = "identity", fill = "blue") +
  coord_flip() +
  labs(title = "Significant Genes based on Raw p-value",
       x = "Gene",
       y = "-log10(p-value)") +
  theme_minimal()

# Plot significant genes based on Bonferroni corrected p-value
plot_significant_bonferroni <- ggplot(significant_bonferroni, aes(x = reorder(Gene, -bonferroni_corrected_p), y = -log10(bonferroni_corrected_p))) +
  geom_bar(stat = "identity", fill = "red") +
  coord_flip() +
  labs(title = "Significant Genes based on Bonferroni Corrected p-value",
       x = "Gene",
       y = "-log10(Bonferroni Corrected p-value)") +
  theme_minimal()

# Plot significant genes based on FDR corrected p-value
plot_significant_fdr <- ggplot(significant_fdr, aes(x = reorder(Gene, -fdr_corrected_p), y = -log10(fdr_corrected_p))) +
  geom_bar(stat = "identity", fill = "green") +
  coord_flip() +
  labs(title = "Significant Genes based on FDR Corrected p-value",
       x = "Gene",
       y = "-log10(FDR Corrected p-value)") +
  theme_minimal()

# Print the plots
print(plot_significant_fisher)
print(plot_significant_bonferroni)
print(plot_significant_fdr)

# Filter out rows with Allele_frequency_gnomAD equal to 0 to avoid division by zero
filtered_df <- df %>%
  filter(Allele_frequency_gnomAD != 0)

# Calculate the allele frequency ratio
filtered_df <- filtered_df %>%
  mutate(allele_frequency_ratio = Allele_Frequency / Allele_frequency_gnomAD)

# Get the top 20 genes with the highest allele frequency ratios
top_genes <- filtered_df %>%
  arrange(desc(allele_frequency_ratio)) %>%
  head(20)

# Visualization of allele frequency ratios
plot_top_genes <- ggplot(top_genes, aes(x = reorder(Gene, -allele_frequency_ratio), y = allele_frequency_ratio, fill = Gene)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top 20 Allele Frequency Ratios (Allele Frequency / Allele Frequency gnomAD)",
       x = "Gene",
       y = "Allele Frequency Ratio") +
  theme_minimal() +
  theme(legend.position = "none", # Remove legend
        axis.title.y = element_blank(), # Remove y-axis title
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) # Add plot margin

# Save the plot as an image
ggsave("top_20_allele_frequency_ratio_plot.png", plot_top_genes)

# Print the data with allele frequency ratios
print("Top 20 data with allele frequency ratios:")
print(top_genes)

# Print the top 20 genes plot
print(plot_top_genes)

# Calculate statistics for allele frequency ratios
allele_freq_stats <- filtered_df %>%
  summarise(
    mean_ratio = mean(allele_frequency_ratio, na.rm = TRUE),
    median_ratio = median(allele_frequency_ratio, na.rm = TRUE),
    sd_ratio = sd(allele_frequency_ratio, na.rm = TRUE),
    min_ratio = min(allele_frequency_ratio, na.rm = TRUE),
    max_ratio = max(allele_frequency_ratio, na.rm = TRUE)
  )

# Print the statistics
print("Statistics for allele frequency ratios:")
print(allele_freq_stats)

# Calculate the number of genes with a positive allele frequency ratio
positive_ratio_count <- filtered_df %>%
  filter(allele_frequency_ratio > 0) %>%
  summarise(n_distinct(Gene)) %>%
  pull()

# Print the count of genes with a positive ratio
print(paste("Number of genes with a positive allele frequency ratio:", positive_ratio_count))

# Calculate the number of genes where the ratio is > 0 but < 1
ratio_between_0_and_1_genes <- filtered_df %>%
  filter(allele_frequency_ratio > 0 & allele_frequency_ratio < 1) %>%
  select(Gene, allele_frequency_ratio) %>%
  distinct() %>%
  arrange(allele_frequency_ratio) %>%
  head(8)

# Print the 8 genes with an allele frequency ratio > 0 but < 1
print("Genes with an allele frequency ratio > 0 but < 1 (top 8):")
print(ratio_between_0_and_1_genes)

# Calculate the number of genes where the ratio is >= 1
ratio_gte_1_count <- filtered_df %>%
  filter(allele_frequency_ratio >= 1) %>%
  summarise(n_distinct(Gene)) %>%
  pull()

# Print the count of genes where the ratio is >= 1
print(paste("Number of genes with an allele frequency ratio >= 1:", ratio_gte_1_count))

# Visualization: Histogram of allele frequency ratios
hist_plot <- ggplot(filtered_df, aes(x = allele_frequency_ratio)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black") +
  scale_x_log10() + # Log scale for better visualization
  labs(title = "Histogram of Allele Frequency Ratios",
       x = "Allele Frequency Ratio (Log Scale)",
       y = "Count") +
  theme_minimal()

# Save the histogram plot
ggsave("allele_frequency_ratio_histogram.png", hist_plot)

# Prepare data for the bar plot of gene counts in different ratio ranges
counts_df <- data.frame(
  range = c("0 < ratio < 1", "ratio >= 1"),
  count = c(nrow(ratio_between_0_and_1_genes), ratio_gte_1_count)
)

# Visualization: Bar plot of gene counts in different ratio ranges
bar_plot <- ggplot(counts_df, aes(x = range, y = count, fill = range)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = rainbow(2)) +
  labs(title = "Counts of Genes in Different Allele Frequency Ratio Ranges",
       x = "Allele Frequency Ratio Range",
       y = "Gene Count") +
  theme_minimal() +
  theme(legend.position = "none") # Remove legend

# Save the bar plot
ggsave("allele_frequency_ratio_bar_plot.png", bar_plot)

# Show the plots
print(hist_plot)
print(bar_plot)

