###########################################
# GTF Processing Script
# Purpose: Extract and process gene mapping information from GTF file
###########################################

# Load required libraries
library(dplyr)
library(tidyverse)
library(tidyr)
library(stringr)

###########################################
# File Import and Initial Processing
###########################################
# Set file path
gtf_file <- "~/pombe-downloads/GCF_000002945.2_ASM294v3_genomic.gtf"

# Read GTF file
gtf_data <- read.delim(gtf_file, 
                      header = FALSE, 
                      comment.char = "#", 
                      sep = "\t", 
                      stringsAsFactors = FALSE)

###########################################
# Gene Mapping Creation
###########################################
# Extract relevant information from GTF
gene_mapping <- gtf_data %>%
  filter(V3 == "transcript") %>%
  mutate(
    # Extract and clean identifiers
    pombase_id = str_extract(V9, "PomBase:[^;]+") %>%
      str_remove("PomBase:") %>%
      toupper(),
    transcript_id = str_extract(V9, "GenBank:[^;]+") %>%
      str_remove("GenBank:"),
    gene = str_extract(V9, "gene\\s+[^;]+") %>%
      str_remove("gene\\s+")
  )

# For use with chromoPlot
gene_mapping_coord <- gtf_data %>%
  filter(V3 == "gene") %>%
  mutate(
    # Extract and clean identifiers
    pombase_id = str_extract(V9, "PomBase:[^;]+") %>%
      str_remove("PomBase:") %>%
      toupper(),
    gene_id = str_extract(V9, "gene_id [^;]+") %>%
      str_remove("gene_id SPOM_"),
    gene = str_extract(V9, "gene\\s+[^;]+") %>%
      str_remove("gene\\s+")
  )

# Add column names
new_colnames <- c("chromosome", "source", "type", "start", "stop", 
                           "V6", "strand", "V8", "annotations",
                           "pombase_id", "transcript_id", "gene")

colnames(gene_mapping) <- c("chromosome", "source", "type", "start", "stop", 
                                  "V6", "strand", "V8", "annotations",
                                  "pombase_id", "transcript_id", "gene")

colnames(gene_mapping_coord) <-c("chr", "source", "type", "start", "stop", 
                                 "V6", "strand", "V8", "annotations",
                                 "pombase_id", "gene_id", "gene")



###########################################
# Data Cleaning and Formatting
###########################################
# Convert chromosome names
gene_mapping <- gene_mapping %>%
  mutate(chromosome = case_when(
    chromosome == "NC_003424.3" ~ "I",
    chromosome == "NC_003423.3" ~ "II",
    chromosome == "NC_003421.2" ~ "III",
    chromosome == "NC_088682.1" ~ "MT1",
    TRUE ~ chromosome
  ))

gene_mapping_coord <- gene_mapping_coord %>%
  mutate(chr = case_when(
    chr == "NC_003424.3" ~ "I",
    chr == "NC_003423.3" ~ "II",
    chr == "NC_003421.2" ~ "III",
    chr == "NC_088682.1" ~ "MT1",
    TRUE ~ chr
  ))

# Select and clean relevant columns
gene_mapping <- gene_mapping %>%
  dplyr::select(transcript_id, pombase_id, gene, chromosome, strand) %>%
  filter(!is.na(transcript_id)) %>%
  mutate(
    gene = ifelse(is.na(gene), pombase_id, gene),
    gene = str_trim(gene),
    pombase_id = str_trim(pombase_id),
    transcript_id = str_trim(transcript_id)
  ) %>%
  distinct()

# For coordinate file: 
gene_mapping_coord <- gene_mapping_coord %>%
  filter(chr != 'MT1') %>%
  dplyr::select(chr, start, stop, strand, gene_id, gene) %>%
  filter(!is.na(gene_id)) %>%
  mutate(
    gene = ifelse(is.na(gene), gene_id, gene),
    gene = str_trim(gene),
    gene_id = str_trim(gene_id)
  ) %>%
  distinct()



###########################################
# Quality Control
###########################################
# Check for duplicates (commented out for future troubleshooting)
# duplicate_check <- list(
#   genes = gene_mapping[duplicated(gene_mapping$gene), ],
#   pombase = gene_mapping[duplicated(gene_mapping$pombase_id), ],
#   transcript = gene_mapping[duplicated(gene_mapping$transcript_id), ]
# )

 duplicate_check_coord <- list(
   genes = gene_mapping_coord[duplicated(gene_mapping_coord$gene), ],
   pombase = gene_mapping_coord[duplicated(gene_mapping_coord$gene_id), ]
 )

gene_mapping_coord <- gene_mapping_coord[-6907, ]

gene_mapping_coord <- gene_mapping_coord %>%
  mutate(gene_id = case_when(
    gene_id == "SPNCRNA.1715.1" ~ "SPNCRNA.1715",
    TRUE ~ gene_id),
    pombase_id = if_else(
      substr(gene_id, nchar(gene_id), nchar(gene_id)) == "C",  # Check if the last letter is "c"
      paste0(substr(gene_id, 1, nchar(gene_id) - 1), "c"),       # Replace the last letter with "C"
      gene_id                                                      # Keep the string unchanged otherwise
    ))

rownames(gene_mapping_coord) <- gene_mapping_coord$pombase_id

###########################################
# Create Final Output
###########################################
# Create simplified version for tximport
gene_mapping_txi <- gene_mapping %>% 
  dplyr::select(transcript_id, pombase_id, gene)

# Clean up environment
rm(gtf_data, duplicate_check_coord,new_colnames)
