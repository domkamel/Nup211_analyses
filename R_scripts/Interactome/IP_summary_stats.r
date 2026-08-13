# ############################################################################
# Nup211 IP-MS Data Filtering and Analysis Pipeline
# ############################################################################
# Purpose: Filter and identify high-confidence protein interaction partners
#          Generate summary statistics of IP data
# ############################################################################
# rm(list = ls())

# Load required libraries
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(ggrepel)

# ############################################################################
# Data Import
# ############################################################################

# Read your data file
# Adjust the file path and delimiter as needed
data <- read.delim("~/Documents/Nup211_IP_MS_Results.txt", 
                   header = TRUE, 
                   stringsAsFactors = FALSE)

# Preview the data
head(data); str(data)

# Read in identifiers
df_identifiers <- read.delim("~/Downloads/gene_IDs_names_products (1).tsv") %>%
  mutate_if(is.character, list(~na_if(.,""))) %>%
  mutate(
    gene_name = coalesce(gene_name, gene_systematic_id)
  ) %>%
  select(c(1,3))

# ############################################################################
# Data Preparation
# ############################################################################

# Create a clean experiment type column
data <- data %>%
  mutate(
    exp_type = case_when(
      grepl("CONTROL", experiment, ignore.case = TRUE) ~ "Control",
      grepl("IP", experiment, ignore.case = TRUE) ~ "IP",
      grepl("NUP211", experiment, ignore.case = TRUE) ~ "IP",
      TRUE ~ "Unknown"
    ),
    exp_short = case_when(
      experiment == "ZHONGHL_1-CONTROL_053124" ~ "Control_1",
      experiment == "ZHONGHL_CONTROL_070524" ~ "Control_2",
      experiment == "ZHONGHL_2-IP-NUPZH_053124" ~ "IP_1",
      experiment == "ZHONGHL_IP-SAMPLE_CONC_070524" ~ "IP_2",
      experiment == "ZHONGHL_NUP211-E1_110424" ~ "IP_3",
      experiment == "ZHONGHL_NUP211-E2_110424" ~ "IP_4_E2",
      TRUE ~ experiment
    ) 
  )

# ############################################################################
# Filter out non-informative ZHONGHL_NUP211-E2_110424
#
# ZHONGHL_NUP211-E2_110424 is Elution 2 (different from other samples - EXCLUDE)
# ############################################################################

cat("\n=== SAMPLE QUALITY CHECK ===\n")
cat("Proteins detected per experiment:\n")
sample_counts <- data %>%
  group_by(experiment, exp_type) %>%
  summarise(n_proteins = n(), .groups = "drop")
print(sample_counts)

# ZHONGHL_CONTROL_070524 only has 1 protein but it's also in Control_1, so it's fine to keep or exclude

# Create filtered dataset excluding E2 ----------------------------------------

data_filtered <- data %>%
  filter(experiment != "ZHONGHL_NUP211-E2_110424")

cat("Original dataset:", nrow(data), "rows\n")
cat("Filtered dataset:", nrow(data_filtered), "rows\n")

# Use filtered data for all downstream analyses  ------------------------------

data <- data_filtered %>%
  select(c(1, 11, 12,2, 8, 4:6, 3, 7, 10)) %>%   # 9 removed
  left_join(df_identifiers, by = "gene_name") %>%
  mutate(pombase_link = paste0('=HYPERLINK("https://www.pombase.org/gene/',
                               gene_systematic_id, '","', gene_name, '")'))
head(data)
data <- data %>% 
  select(c(1:10,13)) #include 11 if you need sequence header info

# Rename columns for clarity  -------------------------------------------------
names(data) <- c("experiment", "exp_type", "exp_short", "hit", "uniprot_id", "no_of_peptides", "gene_name",
                 "gene_product", "molecular_weight","relative_abundance", "pombase_link")

# Convert abundance to numeric (remove % if present)  -------------------------
data$relative_abundance <- as.numeric(gsub("%", "", data$relative_abundance))

# ############################################################################
# Contaminant / non-informative protein list (Work in progress)
# ############################################################################

# Define common MS contaminants ---------------------------------------------
# Customize this list based on experience / literature

contaminant_list <- list(
  # Ribosomal proteins (40S and 60S subunits)
  ribosomal = c("^rps", "^rpl", "^rpp"),
  
  # Metabolic enzymes (common contaminants)
  metabolic = c("tdh1", "gapdh", "glyceraldehyde", "aldolase", "enolase", 
                "pyruvate kinase", "phosphoglycerate"),
  
  # Cytoskeletal proteins
  cytoskeletal = c("act1", "actin", "profilin", "clathrin"),
  
  # Chaperones and heat shock proteins
  chaperones = c("hsp", "heat shock", "ssa", "ssb", "grp"),
  
  # Translation machinery (often non-specific)
  translation = c("eef", "eif", "elongation factor", "initiation factor"),
  
  # Vesicle trafficking (often sticky)
  trafficking = c("arf", "rab", "sec", "cop"),
  
  # Other common contaminants
  other = c("keratin", "trypsin", "albumin")
)

# Function to classify proteins as contaminants  -------------------------------------
classify_contaminant <- function(gene_name, gene_product, contaminant_list) {
  gene_lower <- tolower(paste(gene_name, gene_product))
  
  for (category in names(contaminant_list)) {
    patterns <- contaminant_list[[category]]
    if (any(sapply(patterns, function(p) grepl(p, gene_lower, ignore.case = TRUE)))) {
      return(category)
    }
  }
  return("non_contaminant")
}

# Apply classification   -------------------------------------------------------------
data <- data %>%
  mutate(
    contaminant_type = mapply(classify_contaminant, 
                               gene_name, 
                               gene_product, 
                               MoreArgs = list(contaminant_list = contaminant_list)),
    is_contaminant = contaminant_type != "non_contaminant",
    is_bait = gene_name == "nup211"
  )

# ############################################################################
# Get Summary Statistics
# ############################################################################

# Calculate protein-level statistics across experiments ------------------------------------------
protein_summary <- data %>%
  group_by(gene_name, gene_product, contaminant_type) %>%
  summarise(
    # Detection frequency
    n_total_detections = n(),
    n_IP_detections = sum(exp_type == "IP"),
    n_control_detections = sum(exp_type == "Control"),
    
    # Abundance statistics
    mean_abundance_IP = ifelse(sum(exp_type == "IP") > 0,
                               mean(relative_abundance[exp_type == "IP"], na.rm = TRUE),
                               0),
    mean_abundance_control = ifelse(sum(exp_type == "Control") > 0,
                                    mean(relative_abundance[exp_type == "Control"], na.rm = TRUE),
                                    0),
    max_abundance_IP = if (any(exp_type == "IP"))
      max(relative_abundance[exp_type == "IP"], na.rm = TRUE) else 0,
    
    # Peptide statistics
    mean_peptides_IP = ifelse(sum(exp_type == "IP") > 0,
                              mean(no_of_peptides[exp_type == "IP"], na.rm = TRUE),
                              0),
    max_peptides_IP = ifelse(sum(exp_type == "IP") > 0,
                             max(no_of_peptides[exp_type == "IP"], na.rm = TRUE),
                             0),
    
    # Protein characteristics
    protein_mass = dplyr::first(molecular_weight),
    uniprot_id = dplyr::first(uniprot_id),
    
    .groups = "drop"
  ) %>%
  mutate(
    # Calculate fold enrichment (IP vs Control)
    fold_enrichment = ifelse(mean_abundance_control > 0,
                             mean_abundance_IP / mean_abundance_control,
                             ifelse(mean_abundance_IP > 0, Inf, 0)),
    
    # Calculate specificity score
    specificity = n_IP_detections / (n_IP_detections + n_control_detections + 0.1) 
  ) %>%
  
  # Sort by mean relative absorbance in IP fractions
  arrange(desc(mean_peptides_IP),desc(mean_abundance_IP))

# ############################################################################
# Collapse Replicates (Optional - for data presentation mostly)
# ############################################################################

# This is APPROPRIATE when:
# - wanting to simplify presentation
# - looking at presence/absence rather than quantitative differences
# - comparing "IP" vs "Control" as categories
# - technical variation between replicates is not dramatic

# Create collapsed dataset ---------------------------------------------------
collapsed_data <- data %>%
  group_by(gene_name, gene_product, exp_type) %>%
  summarise(
    # For each protein in each condition (IP or Control)
    n_replicates_detected = n(),
    
    # Abundance metrics
    mean_abundance = mean(relative_abundance, na.rm = TRUE),
    median_abundance = median(relative_abundance, na.rm = TRUE),
    sd_abundance = sd(relative_abundance, na.rm = TRUE),
    max_abundance = max(relative_abundance, na.rm = TRUE),
    
    # Peptide metrics  
    mean_peptides = mean(no_of_peptides, na.rm = TRUE),
    median_peptides = median(no_of_peptides, na.rm = TRUE),
    max_peptides = max(no_of_peptides, na.rm = TRUE),
    sum_peptides = sum(no_of_peptides, na.rm = TRUE),
    
    # Other info
    molecular_weight = dplyr::first(molecular_weight),
    uniprot_id = dplyr::first(uniprot_id),
    contaminant_type = dplyr::first(contaminant_type),
    is_contaminant = dplyr::first(is_contaminant),
    is_bait = dplyr::first(is_bait),
    
    .groups = "drop"
  ) %>%
  # Pivot to have IP and Control as columns for easy comparison
  pivot_wider(
    id_cols = c(gene_name, gene_product, molecular_weight, uniprot_id, 
                contaminant_type, is_contaminant, is_bait),
    names_from = exp_type,
    values_from = c(n_replicates_detected, mean_abundance, median_abundance, 
                    sd_abundance, max_abundance, mean_peptides, median_peptides,
                    max_peptides, sum_peptides),
    values_fill = 0
  ) %>%
  # Calculate enrichment and specificity
  mutate(
    # IP-specific detection (detected in IP but not Control)
    IP_specific = (n_replicates_detected_IP > 0) & (n_replicates_detected_Control == 0),
    Control_specific = (n_replicates_detected_IP == 0) & (n_replicates_detected_Control > 0),
    # Fold enrichment (using mean abundance)
    fold_enrichment = ifelse(mean_abundance_Control > 0,
                             mean_abundance_IP / mean_abundance_Control,
                             ifelse(mean_abundance_IP > 0, Inf, 0)),
    
    # Reproducibility in IPs (out of 3 IP experiments - after excluding E2)
    reproducibility = n_replicates_detected_IP / 3,
    
    # Flag high-confidence hits
    high_confidence = # (n_replicates_detected_IP >= 2) &  # In at least 2 IPs
      (mean_abundance_Control < 5) &      # Low in controls
      (max_peptides_IP >= 2) &            # Good peptide evidence
      (!is_contaminant) | gene_name %in% c("sum3","zuo1","mts4","erf1","fib1")    # Not a contaminant OR gene in custom list
  ) %>%
  
  # Sort by mean relative absorbance in IP fractions
  arrange(desc(mean_peptides_IP), desc(mean_abundance_IP))


# Export collapsed dataset ---------------------------------------------------

# write.csv(collapsed_data, 
#           "~/Documents/Nup211_IP Mass Spec Analysis/Nup211_IP_vs_Control_summary.csv", 
#           row.names = FALSE)


# Create a simplified comparison table for high-confidence hits ------------------------
high_confidence_table <- collapsed_data %>%
  filter(high_confidence) %>%
  select(
    Gene = gene_name,
    Product = gene_product,
    `IP Detections` = n_replicates_detected_IP,
    `Control Detections` = n_replicates_detected_Control,
    `Mean Abundance (IP)` = mean_abundance_IP,
    `Mean Abundance (Control)` = mean_abundance_Control,
    `Fold Enrichment` = fold_enrichment,
    `Mean Peptides (IP)` = mean_peptides_IP,
    `Max Peptides (IP)` = max_peptides_IP,
    `IP Specific` = IP_specific,
    `Control Specific` = Control_specific
  ) %>%
  
  # Sort by mean relative abundance in IP fractions
  arrange(desc(`Mean Peptides (IP)`),desc(`Mean Abundance (IP)`))

# write.csv(high_confidence_table,
#           "~/Documents/Nup211_IP Mass Spec Analysis/Nup211_high_confidence_hits_table.csv",
#           row.names = FALSE)

cat("\n### COLLAPSED REPLICATE ANALYSIS ###\n")
cat("Total unique proteins:", nrow(collapsed_data), "\n")
cat("Detected in control :", sum(collapsed_data$n_replicates_detected_Control > 0), "\n")
cat("Detected in IP :", sum(collapsed_data$n_replicates_detected_IP > 0), "\n")
cat("Detected in Control only :", sum(collapsed_data$Control_specific, na.rm = TRUE), "\n")
cat("Detected in IP only:", sum(collapsed_data$IP_specific, na.rm = TRUE), "\n")

cat("High-confidence hits:", sum(collapsed_data$high_confidence & !collapsed_data$is_bait, na.rm = TRUE), "\n\n")

cat("### HIGH-CONFIDENCE INTERACTION PARTNERS (Collapsed) ###\n")
print(high_confidence_table %>% head(10), row.names = FALSE)


# ############################################################################
# Create "publication-ready" table
# ############################################################################

# Generate a manuscript-ready supplementary table ----------------------------
publication_table <- collapsed_data %>%
  #filter(!is_bait) %>%
  arrange(desc(mean_peptides_IP),(mean_abundance_IP)) %>%
  transmute(
    `Protein` = gene_name,
    `Description` = gene_product,
    `Molecular Weight (kDa)` = (molecular_weight/1000),
    `UniProt ID` = uniprot_id,
    `IP Detections (n=3)` = n_replicates_detected_IP,
    `IP Specific` = IP_specific,
    `Control Detections (n=2)` = n_replicates_detected_Control,
    `Mean Peptides IP` = mean_peptides_IP,
    `Max Peptides IP` = max_peptides_IP,
    `Sum Peptides IP` = sum_peptides_IP,
    `Mean Abundance IP (%)` = round(mean_abundance_IP, 2),
    `Median Abundance IP (%)` = round(median_abundance_IP, 2),
    `Max Abundance IP (%)` = round(max_abundance_IP, 2),
    `SD Abundance IP` = round(sd_abundance_IP, 2),
    `Mean Peptides Control` = mean_peptides_Control,
    `Median Peptides Control` = round(median_peptides_Control, 2),
    `Max Peptides Control` = max_peptides_Control,
    `Sum Peptides Control` = sum_peptides_Control,
    `Mean Abundance Control (%)` = round(mean_abundance_Control, 2),
    `Median Abundance Control (%)` = round(median_abundance_Control, 2),
    `Max Abundance Control (%)` = round(max_abundance_Control, 2),
    `High Confidence` = high_confidence,
    `Flags` = contaminant_type
  )


publication_table <- publication_table %>%
  left_join(df_identifiers, by = c("Protein" = "gene_name")) %>%
  mutate(Link = paste0('=HYPERLINK("https://www.pombase.org/gene/', 
                               gene_systematic_id, '","', 
                               Protein, '")')) 
publication_table <- publication_table %>% select(c(25,1:24))

library(openxlsx2)
# write_xlsx(publication_table,
#            file = "~/Documents/Nup211_IP Mass Spec Analysis/Nup211_supplementary_table.xlsx",
#           row_names = FALSE)

cat("\n=== Publication-ready table created ===\n")
cat("File: Nup211_supplementary_table.csv\n\n")


# ############################################################################
# Calculating % sequence coverage (protein)
# ############################################################################
library(readxl)
library(readr)
library(tidyr)
library(dplyr)
library(stringr)
library(writexl)
library("openxlsx2")


# Function to read and parse FASTA file -------------------------------------
parse_fasta <- function(fasta_file) {
  # Read the file
  lines <- readLines(fasta_file)
  
  # Initialize vectors
  headers <- c()
  sequences <- c()
  current_seq <- ""
  
  # Parse line by line
  for (line in lines) {
    if (grepl("^>", line)) {
      # Header line
      if (current_seq != "") {
        sequences <- c(sequences, current_seq)
        current_seq <- ""
      }
      headers <- c(headers, line)
    } else {
      # Sequence line
      current_seq <- paste0(current_seq, line)
    }
  }
  # Add last sequence
  if (current_seq != "") {
    sequences <- c(sequences, current_seq)
  }
  
  # Create dataframe
  df <- data.frame(
    header = headers,
    sequence = sequences,
    stringsAsFactors = FALSE
  )
  
  return(df)
}

# Read the FASTA file
fasta_df <- parse_fasta("~/Rstudio/peptide.fa")

# Extract gene name from header
# Format: >SPAC1002.01.1:pep mrx11|mitochondrial expression...
fasta_df <- fasta_df %>%
  mutate(
    # Extract systematic ID and remove the .1 at the end
    gene_systematic_id = str_extract(header, "(?<=>)[^:]+") %>%
      str_replace("\\.\\d+$", ""),  # Remove .1, .2, etc. at the end
    
    # Extract gene name (after 'pep ' and before '|')
    protein = str_extract(header, "(?<=pep\\s)[^|]+") %>%
      str_trim(),
    
    # Remove stop codon (*) from sequence
    sequence_clean = str_replace(sequence, "\\*$", ""),
    
    # Calculate sequence length (without stop codon)
    sequence_length = nchar(sequence_clean)
  ) %>%
  select(gene_systematic_id, protein, protein_sequence = sequence_clean, sequence_length) %>%
  mutate_if(is.character, list(~na_if(.,""))) %>%
  mutate(
    protein = coalesce(protein, gene_systematic_id))

# View results
head(fasta_df)

fasta_df <- fasta_df %>%
  mutate(
    protein = case_when(
      protein == "hom3" ~ "SPBC1703.13c",
      protein == "mir1" ~ "SPBC19F5.04",
      TRUE ~ protein
    )
  )



IP_peptides <- read_excel("~/Documents/Nup211_IP Mass Spec Analysis/Nup211_IP_MS_Results.xlsx", sheet = "Peptides")

IP_peptides <- IP_peptides %>%
  mutate(
    exp_type = case_when(
      grepl("CONTROL", experiment, ignore.case = TRUE) ~ "Control",
      grepl("IP", experiment, ignore.case = TRUE) ~ "IP",
      grepl("NUP211", experiment, ignore.case = TRUE) ~ "IP",
      TRUE ~ "Unknown"
    ),
    exp_short = case_when(
      experiment == "ZHONGHL_1-CONTROL_053124" ~ "Control_1",
      experiment == "ZHONGHL_CONTROL_070524" ~ "Control_2",
      experiment == "ZHONGHL_2-IP-NUPZH_053124" ~ "IP_1",
      experiment == "ZHONGHL_IP-SAMPLE_CONC_070524" ~ "IP_2",
      experiment == "ZHONGHL_NUP211-E1_110424" ~ "IP_3",
      experiment == "ZHONGHL_NUP211-E2_110424" ~ "IP_4_E2",
      TRUE ~ experiment
    ) 
  ) %>%
  filter(experiment != "ZHONGHL_NUP211-E2_110424")

IP_peptides <- IP_peptides %>% select(c(9,10,5,4))  %>%
  distinct()

IP_peptides_average <- IP_peptides %>% select(-c(2))  %>%
  distinct()




# Rename columns for clarity
names(IP_peptides) <- c("exp_type", "experiment","protein","peptide_sequence")
names(IP_peptides_average) <- c("exp_type","protein","peptide_sequence")

# ---------- Quick diagnostics: which proteins don't match ----------
unmatched_proteins <- setdiff(unique(IP_peptides$protein), unique(fasta_df$protein))
cat("Number of unique TurboID proteins:", length(unique(IP_peptides$protein)), "\n")
cat("Number not found in FASTA:", length(unmatched_proteins), "\n")
if(length(unmatched_proteins) > 0) {
  cat("Examples of unmatched proteins (up to 20):\n")
  print(head(unmatched_proteins, 20))
}

# ---------- 4) Function to merge intervals (robust) ----------
merge_intervals <- function(starts, ends) {
  starts <- as.integer(starts)
  ends <- as.integer(ends)
  # drop NA or empty
  keep <- !is.na(starts) & !is.na(ends) & length(starts) > 0
  if(length(starts) == 0 || sum(!is.na(starts)) == 0) return(0L)
  # sort
  ord <- order(starts)
  s <- starts[ord]; e <- ends[ord]
  merged_s <- c(); merged_e <- c()
  cur_s <- s[1]; cur_e <- e[1]
  for(i in 2:length(s)) {
    if(s[i] <= cur_e + 1) {
      cur_e <- max(cur_e, e[i])
    } else {
      merged_s <- c(merged_s, cur_s)
      merged_e <- c(merged_e, cur_e)
      cur_s <- s[i]; cur_e <- e[i]
    }
  }
  merged_s <- c(merged_s, cur_s)
  merged_e <- c(merged_e, cur_e)
  sum(as.integer(merged_e - merged_s + 1))
}

# ---------- 5) Locate peptides inside protein sequences ----------
# helper to find all start/end locations for a given peptide in a sequence
find_all_positions <- function(peptide, protein_seq) {
  if(is.na(protein_seq) || is.na(peptide) || peptide == "") return(data.frame(start = integer(0), end = integer(0)))
  # use str_locate_all with fixed matching
  locs <- str_locate_all(protein_seq, fixed(peptide))[[1]]
  if(nrow(locs) == 0) return(data.frame(start = integer(0), end = integer(0)))
  data.frame(start = locs[,1], end = locs[,2])
}

# ---------- 6) Create dataframe with positions for each peptide hit ----------
# join peptides to fasta sequences for matching proteins
joined <- IP_peptides %>%
  left_join(fasta_df %>% select(protein, protein_sequence, sequence_length), by = "protein")

joined_avg <- IP_peptides_average %>%
  left_join(fasta_df %>% select(protein, protein_sequence, sequence_length), by = "protein")


# check how many peptide rows lack a protein sequence after join (i.e. unmatched)
no_seq_rows <- joined %>% filter(is.na(protein_sequence))
cat("Peptide rows with NO protein sequence after join:", nrow(no_seq_rows), "\n")
if(nrow(no_seq_rows) > 0) {
  cat("Example rows without sequence:\n")
  print(head(no_seq_rows, 10))
}

no_seq_rows_avg <- joined_avg %>% filter(is.na(protein_sequence))

library(purrr)
# For rows that have sequences, find positions
positions_df <- joined %>%
  filter(!is.na(protein_sequence)) %>%
  mutate(pos_df = map2(peptide_sequence, protein_sequence, find_all_positions)) %>%
  # unnest the pos_df so each occurrence becomes a row
  unnest(pos_df)

positions_df_avg <- joined_avg %>%
  filter(!is.na(protein_sequence)) %>%
  mutate(pos_df = map2(peptide_sequence, protein_sequence, find_all_positions)) %>%
  # unnest the pos_df so each occurrence becomes a row
  unnest(pos_df)

# If positions_df ends up empty, that means none of the peptides were found as cleaned.
cat("Number of peptide occurrences located in sequences:", nrow(positions_df), "\n")
nrow(joined)


missing_peptides <- anti_join(joined, positions_df, 
                              by = c("protein", "peptide_sequence"))

#write.csv(positions_df, "~/Documents/TurboID Mass Spec Analysis/merged_peptide_positions.csv", row.names = FALSE)
#write.csv(missing_peptides, "~/Documents/TurboID Mass Spec Analysis/missing_peptides.csv", row.names = FALSE)


# ---------- 7) Summarize per experiment/protein ----------

# coverage_df <- positions_df %>%
#   # keep only valid matches
#   filter(!is.na(start), !is.na(end)) %>%
#   
#   # For each peptide, convert interval → vector of residues
#   mutate(residues = purrr::map2(start, end, ~ seq(.x, .y))) %>%
#   
#   # Expand
#   unnest(residues) %>%
#   
#   # Deduplicate overlapping residues within experiment × protein
#   distinct(exp_type, experiment, protein, residues, sequence_length) %>%
#   
#   # Count unique covered residues
#   group_by(exp_type, experiment, protein, sequence_length) %>%
#   summarise(
#     covered_residues = n(),
#     coverage_percent = 100 * covered_residues / sequence_length,
#     .groups = "drop"
#   )
coverage_df <- positions_df %>%
  filter(!is.na(start), !is.na(end)) %>%
  mutate(residues = purrr::map2(start, end, ~ seq(.x, .y))) %>%
  unnest(residues) %>%
  distinct(exp_type, experiment, protein, residues, sequence_length) %>%
  group_by(exp_type, experiment, protein, sequence_length) %>%
  summarise(covered_residues = dplyr::n(), .groups = "drop") %>%
  mutate(coverage_percent = 100 * covered_residues / sequence_length) %>%
  arrange(desc(coverage_percent))

coverage_df <- coverage_df %>% distinct() %>% arrange(desc(coverage_percent)) 


# write.csv(coverage_df, "~/Documents/TurboID Mass Spec Analysis/protein_coverage_results.csv", row.names = FALSE)


# Now do the same for average

coverage_df_avg <- positions_df_avg %>%
  # keep only valid matches
  filter(!is.na(start), !is.na(end)) %>%
  
  # For each peptide, convert interval → vector of residues
  mutate(residues = purrr::map2(start, end, ~ seq(.x, .y))) %>%
  
  # Expand
  unnest(residues) %>%
  
  # Deduplicate overlapping residues within exp_type × protein
  distinct(exp_type, protein, residues, sequence_length) %>%
  
  # Count unique covered residues
  group_by(exp_type, protein, sequence_length) %>%
  summarise(covered_residues = dplyr::n(), .groups = "drop") %>%
  
  mutate(coverage_percent = 100 * covered_residues / sequence_length) %>%
  arrange(desc(coverage_percent))

coverage_df_avg <- coverage_df_avg %>% distinct() %>% arrange(desc(coverage_percent)) 



coverage_summary <- coverage_df %>%
  group_by(protein, sequence_length, experiment) %>%
  # Pivot to wide format for easy comparison
  pivot_wider(
    id_cols = c(protein, sequence_length),
    names_from = experiment,
    values_from = c(covered_residues, coverage_percent),
    values_fill = list(covered_residues = 0, coverage_percent = 0, detected = FALSE)
  )

# # If you want the actual mean across replicates, aggregate coverage_df instead:
# 
# coverage_df_avg <- coverage_df %>%
#   group_by(exp_type, protein, sequence_length) %>%
#   summarise(
#     n_experiments   = dplyr::n(),
#     mean_coverage   = mean(coverage_percent),
#     max_coverage    = max(coverage_percent),
#     .groups = "drop"
#   )

avg_coverage_summary <- coverage_df_avg %>%
  group_by(protein, sequence_length, exp_type) %>%
  # Pivot to wide format for easy comparison
  pivot_wider(
    id_cols = c(protein, sequence_length),
    names_from = exp_type,
    values_from = c(covered_residues, coverage_percent),
    values_fill = list(covered_residues = 0, coverage_percent = 0, detected = FALSE)
  )




merge_intervals_safe <- function(df) {
  df <- df %>% filter(!is.na(start), !is.na(end))
  
  if (nrow(df) == 0) return(list())  # no intervals
  
  df <- df[order(df$start), ]
  
  merged <- list()
  s <- df$start[1]
  e <- df$end[1]
  
  for (i in 2:nrow(df)) {
    s_i <- df$start[i]
    e_i <- df$end[i]
    
    # skip any unexpected NA
    if (is.na(s_i) || is.na(e_i)) next
    
    if (s_i <= e + 1) {
      e <- max(e, e_i)
    } else {
      merged <- append(merged, list(c(s, e)))
      s <- s_i
      e <- e_i
    }
  }
  
  merged <- append(merged, list(c(s, e)))
  merged
}



positions_df_avg

library(dplyr)
library(purrr)
library(tidyr)

collapsed_avg_df <- positions_df_avg %>%
  group_by(exp_type, protein) %>%
  summarise(
    intervals = list(merge_intervals_safe(cur_data())),
    peptides = list(unique(peptide_sequence)),
    .groups = "drop"
  ) %>%
  mutate(
    # format intervals like "5-20; 30-45"
    intervals = map_chr(intervals, ~ 
      paste(purrr::map_chr(.x, ~ paste(.x[1], .x[2], sep = "-")), collapse = "; ")
    ),
    
    # collapse peptides into comma-separated string
    peptides = map_chr(peptides, ~ paste(sort(.x), collapse = ", ")),
    
    # number of intervals
    n_intervals = map_int(intervals, length),
    
    # number of unique peptides (if useful)
    n_peptides = map_int(peptides, length)
  ) %>%
  select(
    exp_type, protein,
    intervals, peptides,
    n_intervals, n_peptides
  )


library(dplyr)
library(purrr)

merge_intervals_safe <- function(df) {
  df2 <- df %>% filter(!is.na(start), !is.na(end))
  if (nrow(df2) == 0) return(list())
  df2 <- df2[order(df2$start), ]
  
  merged <- list()
  s <- df2$start[1]
  e <- df2$end[1]
  
  for (i in 2:nrow(df2)) {
    s_i <- df2$start[i]
    e_i <- df2$end[i]
    if (is.na(s_i) || is.na(e_i)) next
    
    if (s_i <= e + 1) {
      e <- max(e, e_i)
    } else {
      merged <- append(merged, list(c(s, e)))
      s <- s_i; e <- e_i
    }
  }
  
  merged <- append(merged, list(c(s, e)))
  merged
}

rows <- positions_df_avg %>%
  # ensure input columns exist and clean NA peptide entries for counting
  mutate(peptide_sequence = ifelse(is.na(peptide_sequence), NA_character_, peptide_sequence)) %>%
  group_by(exp_type, protein) %>%
  group_map(.f = function(df, key) {
    # df: subset for this experiment+protein
    # key: tibble with experiment & protein
    ex <- key$exp_type; prot <- key$protein
    
    # unique peptides for this protein (drop NA)
    unique_peps <- unique(na.omit(df$peptide_sequence))
    
    # merged intervals from rows that have start/end
    intervals <- merge_intervals_safe(df)
    
    # for each interval, find unique peptides overlapping it
    if (length(intervals) == 0) {
      interval_str <- ""
      peptides_per_interval <- ""
      n_intervals <- 0L
    } else {
      interval_strs <- map_chr(intervals, ~ paste(.x[1], .x[2], sep = "-"))
      # for each interval, collect peptides that overlap: peptide start <= end_interval & peptide end >= start_interval
      pep_lists <- map(intervals, function(iv) {
        s_iv <- iv[1]; e_iv <- iv[2]
        overlapping <- df %>%
          filter(!is.na(start) & !is.na(end) & start <= e_iv & end >= s_iv) %>%
          pull(peptide_sequence)
        overlapping_unique <- unique(na.omit(overlapping))
        # collapse to single string (comma separated); empty -> ""
        if (length(overlapping_unique) == 0) "" else paste(overlapping_unique, collapse = ",")
      })
      interval_str <- paste(interval_strs, collapse = "; ")
      peptides_per_interval <- paste(pep_lists, collapse = " ; ")  # semicolon separates intervals, comma separates peptides
      n_intervals <- length(intervals)
    }
    
    tibble(
      exp_type = ex,
      protein = prot,
      interval_str = interval_str,
      peptides_per_interval = peptides_per_interval,
      n_intervals = n_intervals,
      n_unique_peptides = length(unique_peps),
      peptides_all = if (length(unique_peps) == 0) "" else paste(sort(unique_peps), collapse = ",")
    )
  }, .keep = TRUE) %>%
  bind_rows()

# rows is the final collapsed dataframe
collapsed_avg_df <- rows

# quick sanity checks:
# 1) see top rows
print(head(collapsed_avg_df, 20))
# 2) show any proteins where n_intervals == 0 but have peptides (diagnostic)
diag <- collapsed_avg_df %>% filter(n_intervals == 0 & n_unique_peptides > 0)
if (nrow(diag) > 0) {
  message("Warning: these proteins have peptides but no valid start/end intervals (check missing starts/ends):")
  print(diag)
}

collapsed_avg_df

coverage_df_avg
avg_coverage_summary


wide_peptide_summary <- collapsed_avg_df %>%
  group_by(protein, exp_type) %>%
  # Pivot to wide format for easy comparison
  pivot_wider(
    id_cols = c(protein),
    names_from = exp_type,
    values_from = c(n_unique_peptides, peptides_all, n_intervals, peptides_per_interval),
    values_fill = list(n_unique_peptides = 0, peptides_all = "none" ,n_unique_peptides = 0, peptides_per_interval = "none",
                       n_intervals = 0)
  )

# write.csv(coverage_df_avg, "~/Documents/Nup211_IP Mass Spec Analysis/IP_avg_protein_coverage.csv", row.names = FALSE)
# write.csv(collapsed_avg_df, "~/Documents/Nup211_IP Mass Spec Analysis/IP_avg_peptide_regions.csv", row.names = FALSE)

# write.csv(avg_coverage_summary, "~/Documents/Nup211_IP Mass Spec Analysis/IP_avg_peptide_coverage.csv", row.names = FALSE)
# write.csv(wide_peptide_summary, "~/Documents/Nup211_IP Mass Spec Analysis/IP_avg_peptide_regions_wide.csv", row.names = FALSE)

IP_peptide_analysis <- avg_coverage_summary %>%
  left_join(wide_peptide_summary, by = c("protein")) %>%
  arrange(desc(n_unique_peptides_IP))

# write.csv(IP_peptide_analysis, "~/Documents/Nup211_IP Mass Spec Analysis/IP_avg_peptide_analysis.csv", row.names = FALSE)


IP_data_plotting <- (as.data.frame(publication_table) %>% select(c(2,7,23,24,12,9,10))) %>%
  left_join(IP_peptide_analysis, by = c("Protein"="protein")) %>%
  arrange(desc(`Mean Abundance IP (%)`))

#then filter out contaminants and those with 0 in percent abundance in IP

IP_filtered_plotting <- IP_data_plotting %>%
  filter(
    # CRITERION 1: Must be detected in IP samples
    `Mean Abundance IP (%)` > 0,
    
    # CRITERION 2:  Exclude classified contaminants (except bait)
    Flags == "non_contaminant" | Protein == "nup211",
    
    # CRITERION 3: Exclude the bait itself from interaction partners
    !Protein %in% c("nup211", "zuo1", "erf1", "pma1", "shm1")  
    )

# install.packages("ggbreak")
library(ggbreak)
library(patchwork)
library(ggrepel)


IP_filtered_plotting_2 <- IP_filtered_plotting %>%
  mutate(label = ifelse(Protein %in% c("spo15", "ubi2", "cdc8"), Protein, ""))

# Plot with visible annotation
test_plot <- ggplot(IP_filtered_plotting_2, aes(x = coverage_percent_IP, y = `Mean Abundance IP (%)`, label = Protein)) +
  geom_point(shape = 21, size = 5, color = "black", fill = "steelblue", alpha = 0.8, stroke = 1) +
  geom_text_repel(aes(label = label),
                  size = 7,
                  fontface = "plain",
                  box.padding = 0.8,
                  max.overlaps = 20) +
  labs(
    title = "IP: Abundance vs. peptide #",
    x = "Coverage (%)",
    y = "Mean Abundance IP (%)"
  ) +
  # annotate("text", x = -2, y = 3,
  #          label = paste0("r = ", cor_value, "\n", "p = ", p_value),
  #          size = 5, fontface = "bold", hjust = 0) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 16, colour = "black"),
    axis.text.y = element_text(size = 16, colour = "black"),
    axis.title.x = element_text(size = 18, face = "bold", colour = "black"),
    axis.title.y = element_text(size = 18, face = "bold", colour = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 2)
   )  # +
      # scale_y_break(breaks = c(5, 10)) + # Breaks the y-axis between 10 and 45 
      # scale_x_break(breaks = c(10, 45)) # Breaks the x-axis between 10 and 45 

test_plot

IP_filtered_inset  <- IP_filtered_plotting %>%
  filter(
    # CRITERION 1: Must be detected in IP samples
    `Mean Abundance IP (%)` < 5 & coverage_percent_IP < 10, 
    
    # CRITERION 2:  Exclude classified contaminants (except bait)
    Flags == "non_contaminant" | Protein == "nup211",
    
    # CRITERION 3: Exclude the bait itself from interaction partners
    !Protein %in% c("nup211", "zuo1", "erf1", "pma1", "shm1")  
  )

inset_plot <- ggplot(IP_filtered_inset, aes(x = coverage_percent_IP, y = `Mean Abundance IP (%)`, label = Protein)) +
  geom_point() +
  geom_point(shape = 21, size = 5, color = "black", fill = "steelblue", alpha = 0.8, stroke = 1) +
  geom_text_repel(aes(label = Protein),
                  size = 7,
                  fontface = "plain",
                  box.padding = 0.8,
                  max.overlaps = 20) +
  coord_cartesian(xlim = c(0, 5), ylim = c(0, 3)) + # Zoom into specific region
  labs(title = "Inset (Zoomed)") +
  labs(
    title = "Inset",
    x = "Coverage (%)",
    y = "Mean Abundance IP (%)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = NULL, #element_text(size = 14, face = "bold"),
    axis.text.x =  element_text(size = 16, colour = "black"),
    axis.text.y =  element_text(size = 16, colour = "black"),
    axis.title.x =  element_text(size = 18, face = "bold", colour = "black"),
    axis.title.y = element_text(size = 18, face = "bold", colour = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 2)
  ) +
  theme(legend.position = "none") # Remove legend from inset for cleaner look


print(test_plot)
print(inset_plot)

library(grid)
# Add inset
vps <- viewport(width = 0.57, height = 0.59, x = 0.35, y = 0.35)
print(inset_plot, vp = vps)

library(cowplot)
library(imager)

ggdraw() +
  draw_plot(test_plot) +
  draw_plot(inset_plot, x = 0.08, y = 0.1, width = 0.45, height = 0.45)


border <- rectGrob(
  gp = gpar(col = "black", fill = NA, lwd = 2, lty = "dotted")  # lwd = border thickness
)

ggdraw() +
  draw_plot(test_plot) +
  draw_plot(inset_plot, x = 0.08, y = 0.1, width = 0.45, height = 0.45) +
  draw_grob(border, x = 0.08, y = 0.1, width = 0.45, height = 0.45)

p1 <- ggdraw() +
  draw_plot(test_plot) +
  draw_grob(border, x = 0.11, y = 0.11, width = 0.2, height = 0.2)

p1+inset_plot

plot_grid(p1, inset_plot, ncol = 2, rel_widths = c(2, 2))



Summary_Full <- read_excel("~/Documents/Nup211_IP Mass Spec Analysis/Nup211_supplementary_table_2.xlsx", sheet = "Sheet2")

names(Summary_Full) <- c("protein",
                         "description",
                         "nup211_mean_abundance",
                         "nup211_peptides",
                         "nup211_coverage",
                         "igg_mean_abundance",
                         "igg_peptides",
                         "igg_coverage")


library(dplyr)
library(flextable)

ft <- flextable(Summary_Full)

ft <- colformat_num(
  ft,
  j = c("nup211_mean_abundance", "igg_mean_abundance", 
        "nup211_coverage", "igg_coverage"),
  digits = 2
)


# Rename visible column labels
ft <- set_header_labels(
  ft,
  protein = "Protein",
  description = "Description",
  nup211_mean_abundance = "Mean relative abundance (%)",
  nup211_peptides = "# peptides",
  nup211_coverage = "Protein coverage (%)",
  igg_mean_abundance = "Mean relative abundance (%)",
  igg_peptides = "# peptides",
  igg_coverage = "Protein coverage (%)"
)

# Add grouped headers
ft <- add_header_row(
  ft,
  values = c("", "", "Nup211 AP", "IgG AP"),
  colwidths = c(1, 1, 3, 3)
)

# Formatting tweaks
ft <- align(ft, align = "center", part = "all")
ft <- align(ft, j = c("protein", "description"), align = "left")
ft <- autofit(ft)

# Optional: nicer borders
ft <- theme_booktabs(ft)

ft <- align(
  ft,
  i = 1,
  align = "center",
  part = "header"
)


ft <- bg(
  ft,
  i = band_rows,
  j = -c("protein"),   # exclude protein column
  bg = "#F2F2F2",
  part = "body"
)

ft <- bg(
  ft,
  i = seq(1, nrow(Summary_Full), by = 2),
  bg = "#F2F2F2",
  part = "body"
)

library(officer)
vline <- fp_border(color = "black", width = 0.75)

ft <- border(
  ft,
  j = c("description","nup211_mean_abundance","igg_mean_abundance"),
  border.left = vline,
  part = "all"
)

ft <- colformat_num(
  ft,
  j = grep("Protein coverage (%)", names(df), value = TRUE),
  digits = 2
)

str(Summary_Full$nup211_coverage)

# save_as_docx(ft, path = "Nup211_AP_LCMS_results.docx")

