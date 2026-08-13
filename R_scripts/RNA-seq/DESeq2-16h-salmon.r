#################################################
# RNA-seq Analysis Pipeline (16h) - Script 1
# Purpose: Analyze RNA-seq data for 16h timepoint
#################################################

# Clear workspace
# rm(list=ls())

###########################################
# Setup and Library Loading
###########################################
# Load required libraries
libraries <- c(
  "readr", "dplyr", "tidyr", "tidyverse", "tximport",
  "GenomicFeatures", "txdbmaker", "AnnotationDbi", "AnnotationForge",
  "GenomicAlignments", "BiocParallel", "DESeq2",
  "apeglm", "ggplot2", "plotly", "gplots",
  "pheatmap", "RColorBrewer", "viridis",
  "EnhancedVolcano", "PoiClaClu", "ggpubr"
)

for (lib in libraries) {
  if (!require(lib, character.only = TRUE)) {
    BiocManager::install(lib)
    library(lib, character.only = TRUE)
  }
}

# Set the prefix for each output file name
outputFolder <- "~/Rstudio/expression_analysis_16h/qc_plots/"

###########################################
# Data Import
###########################################
# Set paths
gtf_file <- "~/pombe-downloads/GCF_000002945.2_ASM294v3_genomic.gtf"
quant_dir <- "~/salmon_output"

# Create transcript database
txdb <- makeTxDbFromGFF(gtf_file)
tx2gene <- transcripts(txdb, columns = c("GENEID", "TXNAME"))

# Define samples and files
samples <- paste0("Sample_", c(1:3,10:12))
files <- file.path(quant_dir, samples, "quant.sf")
names(files) <- samples

#################################################
# Import count data (RUN gtf.processing.r script)
#################################################

# Remember to run gtf.processing.r (dont clear workspace)

# Import count data 
txi <- tximport(files, type = "salmon", tx2gene = gene_mapping_txi)

# Count data (save this and/or normalized counts from vst or rld). Round or dont round?

txi$counts    ## txi has 12380 rows...

# To match my original analyses using the full batch of samples analyzed by RNA-seq 

counts_16h <- round(txi$counts)


############################################################################################
# QC check on what was imported (12380 rows but ~7,000 annotated genes in Pombase) 🤔
############################################################################################

# How many rows are actually protein-coding systematic IDs?
  sum(grepl("^SP[AB]C|^SP[AB]P|^SPBC|^SPCC", rownames(txi$counts)))

# What do the non-standard rownames look like?
  head(rownames(txi$counts)[!grepl("^SP", rownames(txi$counts))], 20)

# Any duplicate gene IDs, which would mean collapsing failed?
  sum(duplicated(rownames(txi$counts)))

# How many rows are empty across all samples?
  sum(rowSums(txi$counts) == 0)

# And the mapping itself
  length(unique(gene_mapping_txi$GENEID))
  nrow(gene_mapping_txi)

# What are the remaining types of transcripts? 
  table(sub("^(SP[A-Z]+).*", "\\1", rownames(txi$counts)))


###########################################
# Sample Metadata Setup
###########################################
# Create metadata table
samples_df <- data.frame(
  row.names = colnames(counts_16h),
  condition = c(
    rep("wild_type", 3),
    rep("nup211so_16h", 3)
  )
)

samples_df$condition <- factor(samples_df$condition,
                                 levels = c("wild_type", "nup211so_16h"))
###########################################
# Initial DESeq2 Analysis
###########################################
# Create DESeq dataset

dds_16h <- DESeqDataSetFromMatrix(
  countData = counts_16h,
  colData   = samples_df,
  design    = ~ condition
)

#################################################
# Export counts data as .csv file
#################################################

# Rename columns to make a bit more sense
colnames(counts_16h) <- c("WT_16h_rep1", "WT_16h_rep2", "WT_16h_rep3",
                          "nup211so_16h_rep1", "nup211so_16h_rep2", "nup211so_16h_rep3")

# Fix the last "C" and make it "c" to match PomBase conventions
rownames(counts_16h) <- sub("C$", "c", rownames(counts_16h))

# Move gene_id's to a real column and drop rownames
counts_out <- data.frame(gene_id = rownames(counts_16h),
                         counts_16h,
                         check.names = FALSE,
                         row.names = NULL)

# Save counts as a .csv file
## write.csv(counts_out, "~/Rstudio/expression_analysis_16h/nup211so_16h_counts.csv", row.names = FALSE)

###########################################
# Quality Control and Visualization (Total)
###########################################
# VST transformation ----
vsd <- vst(dds_16h)

# Create PCA plot ----
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + theme_bw() +  
  labs(title = bquote(bold("PCA Plot:" ~ bolditalic("nup211-so") ~ "(16h) vs. wild-type"))) +
  scale_color_manual(name = "Condition",
                     values = c("wild_type_total" = "#52d8da", "16h_total" = "#F8766D"), 
                     labels = c("wild-type", bquote(italic("nup211-so") ~ "(16h)"))) +
  theme(plot.title = element_text(size = 14, face = "bold"),
        strip.text = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 12, face= "plain",colour = 'black'),    # Axis labels
        axis.text.y = element_text(size = 12, face= "plain",colour = 'black'),
        axis.title.x = element_text(size = 12, face= "bold",colour = 'black'),   # Axis titles
        axis.title.y = element_text(size = 12, face="bold", colour = 'black'),
        legend.title = element_text(size = 12, face = "bold"),   # Adjust legend title font
        legend.text = element_text(size = 12)                  # Adjust legend labels font
  )

#ggsave(pca_plot,file=paste0(outputFolder, "pcaplot_total.tiff"))

# Perform MDS on manually calculated distances, using the R function cmdscale: ----
# Starting from the VST counts, using Euclidean distance.
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix( sampleDists )
mdsData <- data.frame(cmdscale(sampleDistMatrix))
mds <- cbind(mdsData, as.data.frame(colData(vsd)))
ggplot(mds, aes(X1,X2,color=condition)) + geom_point(size=3) + theme_bw() +  
  labs(title = bquote(bold("MDS Plot:" ~ bolditalic("nup211-so") ~ "(16h) vs. wild-type")),
       subtitle = bquote(bolditalic("Euclidean distance"))) +
  scale_color_manual(name = "Condition",
                     values = c("wild_type" = "#52d8da", "nup211so_16h" = "#F8766D"),
                     labels = c("wild-type", bquote(italic("nup211-so") ~ "(16h)"))) +
  theme(plot.title = element_text(size = 14, face = "bold"),
        strip.text = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 12, face= "plain",colour = 'black'),    # Axis labels
        axis.text.y = element_text(size = 12, face= "plain",colour = 'black'),
        axis.title.x = element_text(size = 12, face= "bold",colour = 'black'),   # Axis titles
        axis.title.y = element_text(size = 12, face="bold", colour = 'black'),
        legend.title = element_text(size = 12, face = "bold"),   # Adjust legend title font
        legend.text = element_text(size = 12)                  # Adjust legend labels font
  )
#ggsave(file=paste0(outputFolder, "mds_euclidean_total.tiff"))


# Same plot with Poisson Distance:
poisd <- PoissonDistance(t(counts(dds_16h)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
mdsPoisData <- data.frame(cmdscale(samplePoisDistMatrix))
mdsPois <- cbind(mdsPoisData, as.data.frame(colData(dds_16h)))
ggplot(mdsPois, aes(X1,X2,color=condition)) + geom_point(size=3) + theme_bw() +  
  labs(title = bquote(bold("MDS Plot:" ~ bolditalic("nup211-so") ~ "(16h) vs. wild-type")),
       subtitle = bquote(bolditalic("Poisson distance"))) +
  scale_color_manual(name = "Condition",
                     values = c("wild_type" = "#52d8da", "nup211so_16h" = "#F8766D"), 
                     labels = c("wild-type", bquote(italic("nup211-so") ~ "(16h)"))) +
  theme(plot.title = element_text(size = 14, face = "bold"),
        strip.text = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 12, face= "plain",colour = 'black'),    # Axis labels
        axis.text.y = element_text(size = 12, face= "plain",colour = 'black'),
        axis.title.x = element_text(size = 12, face= "bold",colour = 'black'),   # Axis titles
        axis.title.y = element_text(size = 12, face="bold", colour = 'black'),
        legend.title = element_text(size = 12, face = "bold"),   # Adjust legend title font
        legend.text = element_text(size = 12)                  # Adjust legend labels font
  )
#ggsave(file=paste0(outputFolder, "mds_poisson_total.tiff"))

# Transform raw counts into normalized values by rlog transformation ----
rld <- rlogTransformation(dds_16h, blind=T)
# save normalized values
rld.df <- as.data.frame(assay(rld))
vsd.df <- as.data.frame(assay(vsd))

# Look at sample distances using gplots ----
sampleDists <- dist(t(assay(rld))) 
sampleDistMatrix <- as.matrix(sampleDists) 

rownames(sampleDistMatrix) <- paste(rld$condition)
colnames(sampleDistMatrix) <- NULL

colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) 
labels <- c("wild-type 1", "wild-type 2", "wild-type 3",
            expression(italic("nup211-so") ~ "1"), 
            expression(italic("nup211-so") ~ "2"), 
            expression(italic("nup211-so") ~ "3"))
# Plot heatmap
distmap <- pheatmap(sampleDistMatrix, 
                    color = colours, 
                    labels_row = labels, 
                    labels_col = labels,
                    main = "Sample Distance Matrix \n",
                    annotation_legend = TRUE,
                    legend_labels = "Distance")

#ggsave(distmap,file=paste0(outputFolder, "sample_distances_total.tiff"))



###########################################
# Differential Expression Analysis (Total)
###########################################
# Filter low counts and run DESeq ----
dds_16h <- dds_16h[rowSums(counts(dds_16h)) > 10, ]
dds_16h <- DESeq(dds_16h)

# Get results with different shrinkage methods
resNorm <- lfcShrink(dds_16h, 
                     coef = "condition_nup211so_16h_vs_wild_type", 
                     type = "normal",
                     lfcThreshold = 0) # or lfcThreshold = 1

###########################################
# Gene Annotation and Result Processing
###########################################
# Import gene descriptions and clean up
geneIDs_descriptions <- read.delim("~/Downloads/gene_IDs_names_products (1).tsv") %>%
  mutate_if(is.character, list(~na_if(.,""))) %>%
  mutate(
    gene_name = coalesce(gene_name, gene_systematic_id),
    gene_systematic_id = toupper(gene_systematic_id),
    gene_systematic_id_with_prefix = gsub("PomBase:", "", gene_systematic_id_with_prefix)
  )

# Merge results with annotations
results_annotated <- merge(
  as.data.frame(resNorm) %>% 
    mutate(Gene_ID = rownames(.)),
  geneIDs_descriptions,
  by.x = "Gene_ID",
  by.y = "gene_systematic_id"
) %>%
  distinct() %>%
  column_to_rownames("gene_systematic_id_with_prefix") %>%
  mutate(
    pvalue = case_when(
      (is.na(pvalue) | is.infinite(pvalue)) & abs(log2FoldChange) >= 1 ~ 1.663234e-35,
      (is.na(pvalue) | is.infinite(pvalue)) & abs(log2FoldChange) < 1 ~ 0.051,
      TRUE ~ pvalue
    ),
    padj = case_when(
      (is.na(padj) | is.infinite(padj)) & abs(log2FoldChange) >= 1 ~ 1.663234e-35,
      (is.na(padj) | is.infinite(padj)) & abs(log2FoldChange) < 1 ~ 0.051,
      TRUE ~ padj
    )
  )

results_annotated$PomBase_Link <- paste0(
  '<a href="https://www.pombase.org/gene/', rownames(results_annotated), '" target="_blank">', results_annotated$Gene_ID, '</a>'
)
results_annotated$tooltip <- paste0(
  "Gene:", results_annotated$gene_name, "<br>",
  "Log2FC: ", results_annotated$log2FoldChange, "<br>",
  "Padj: ", results_annotated$padj
)

results_annotated <- merge(results_annotated, gene_mapping_coord,
                           by = 0) %>%
  distinct() %>% 
  column_to_rownames("pombase_id") %>% # Set rownames to pombase_id but keep it as a column
  dplyr::select(c(2,9,3:8,11,13,17:20,15,16))


filtered_annotated <- results_annotated %>%
  # Keep row if either:
  # 1. Gene_ID doesn't contain SPNCRNA
  # 2. Gene_ID contains SPNCRNA but gene_name doesn't (i.e., it has a special name)
  filter(
    !str_detect(Gene_ID, "SPNCRNA") | 
      (str_detect(Gene_ID, "SPNCRNA") & !str_detect(gene_name, "SPNCRNA"))
  )

filtered_annotated <- results_annotated %>% filter(gene_type == "protein coding gene")

outputFolder <- "~/Rstudio/expression_analysis_16h/results/DE_spreadsheets/"

#write.csv(results_annotated,file=paste0(outputFolder, "DESeq2_results_all_transcripts.csv"))
#write.csv(filtered_annotated,file=paste0(outputFolder, "DESeq2_results_protein_coding.csv"))


###########################################
# Result Processing and Annotation
###########################################
# Process results
process_results <- function(res) {
  as.data.frame(res) %>%
    rownames_to_column("Gene_ID") %>%
    mutate(Gene_ID = sub("C$", "c", Gene_ID)) %>%
    arrange(log2FoldChange) %>%
    column_to_rownames("Gene_ID")
}

res_normal_df <- process_results(resNorm)

###########################################
# Data Summary
###########################################

# padj cutoff = 0.05
cat("all transcripts (padj cutoff = 0.05):", "\n",
    "2-Fold Upregulated:", sum(results_annotated$padj <= 0.05 & results_annotated$log2FoldChange >= 1), "\n",
    "2-Fold Downregulated:", sum(results_annotated$padj <= 0.05 & results_annotated$log2FoldChange <= -1), "\n",
    "4-Fold Upregulated:", sum(results_annotated$padj <= 0.05 & results_annotated$log2FoldChange >= 2), "\n",
    "4-Fold Downregulated:", sum(results_annotated$padj <= 0.05 & results_annotated$log2FoldChange <= -2), "\n",
    "2-Fold DE:", sum(results_annotated$padj <= 0.05 & abs(results_annotated$log2FoldChange) >=1), "\n",
    "4-Fold DE:", sum(results_annotated$padj <= 0.05 & abs(results_annotated$log2FoldChange) >=2), "\n") 

cat("coding genes only (padj cutoff = 0.05):", "\n",
    "2-Fold Upregulated:", sum(filtered_annotated$padj <= 0.05 & filtered_annotated$log2FoldChange >= 1), "\n",
    "2-Fold Downregulated:", sum(filtered_annotated$padj <= 0.05 & filtered_annotated$log2FoldChange <= -1), "\n",
    "4-Fold Upregulated:", sum(filtered_annotated$padj <= 0.05 & filtered_annotated$log2FoldChange >= 2), "\n",
    "4-Fold Downregulated:", sum(filtered_annotated$padj <= 0.05 & filtered_annotated$log2FoldChange <= -2), "\n",
    "2-Fold DE:", sum(filtered_annotated$padj <= 0.05 & abs(filtered_annotated$log2FoldChange) >=1), "\n",
    "4-Fold DE:", sum(filtered_annotated$padj <= 0.05 & abs(filtered_annotated$log2FoldChange) >=2), "\n") 


# padj cutoff = 10e-6
cat("all transcripts (padj cutoff = 10e-6):", "\n",
    "2-Fold Upregulated:", sum(results_annotated$padj <= 10e-6 & results_annotated$log2FoldChange >= 1), "\n",
    "2-Fold Downregulated:", sum(results_annotated$padj <= 10e-6 & results_annotated$log2FoldChange <= -1), "\n",
    "4-Fold Upregulated:", sum(results_annotated$padj <= 10e-6 & results_annotated$log2FoldChange >= 2), "\n",
    "4-Fold Downregulated:", sum(results_annotated$padj < 10e-6 & results_annotated$log2FoldChange <= -2), "\n",
    "2-Fold DE:", sum(results_annotated$padj <= 10e-6 & abs(results_annotated$log2FoldChange) >=1), "\n",
    "4-Fold DE:", sum(results_annotated$padj <= 10e-6 & abs(results_annotated$log2FoldChange) >=2), "\n") 

cat("coding genes only (padj cutoff = 10e-6):", "\n",
    "2-Fold Upregulated:", sum(filtered_annotated$padj <= 10e-6 & filtered_annotated$log2FoldChange >= 1), "\n",
    "2-Fold Downregulated:", sum(filtered_annotated$padj <= 10e-6 & filtered_annotated$log2FoldChange <= -1), "\n",
    "4-Fold Upregulated:", sum(filtered_annotated$padj <= 10e-6 & filtered_annotated$log2FoldChange >= 2), "\n",
    "4-Fold Downregulated:", sum(filtered_annotated$padj <= 10e-6 & filtered_annotated$log2FoldChange <= -2), "\n",
    "2-Fold DE:", sum(filtered_annotated$padj <= 10e-6 & abs(filtered_annotated$log2FoldChange) >=1), "\n",
    "4-Fold DE:", sum(filtered_annotated$padj <= 10e-6 & abs(filtered_annotated$log2FoldChange) >=2), "\n") 


###########################################
# Results Summary Table 
###########################################

library(knitr)  # For kable()
library(flextable)  # Optional for better table formatting
library(officer)

# Function to summarize DE results into a dataframe ----
summarize_DE <- function(data, label, padj_cutoff) {
  data.frame(
    Dataset = label,
    padj_Cutoff = padj_cutoff,
    `2-Fold Up` = sum(data$padj <= padj_cutoff & data$log2FoldChange >= 1),
    `2-Fold Down` = sum(data$padj <= padj_cutoff & data$log2FoldChange <= -1),
    `4-Fold Up` = sum(data$padj <= padj_cutoff & data$log2FoldChange >= 2),
    `4-Fold Down` = sum(data$padj <= padj_cutoff & data$log2FoldChange <= -2),
    `2-Fold DE` = sum(data$padj <= padj_cutoff & abs(data$log2FoldChange) >= 1),
    `4-Fold DE` = sum(data$padj <= padj_cutoff & abs(data$log2FoldChange) >= 2)
  )
}

# Create summary tables for each dataset and padj cutoff
summary_table <- bind_rows(
  summarize_DE(results_annotated, "All Transcripts", 0.05),
  summarize_DE(results_annotated, "All Transcripts", 10e-6),
  summarize_DE(filtered_annotated, "Coding Genes Only", 0.05),
  summarize_DE(filtered_annotated, "Coding Genes Only", 10e-6)
)

summary_table <- summary_table %>%
  rename( `p.adjusted cutoff` = padj_Cutoff,
          `2-Fold Up` = `X2.Fold.Up`,
          `2-Fold Down` = `X2.Fold.Down`,
          `4-Fold Up` = `X4.Fold.Up`,
          `4-Fold Down` = `X4.Fold.Down`,
          `2-Fold DE` = `X2.Fold.DE`,
          `4-Fold DE` = `X4.Fold.DE`
  )

# Display as a nicely formatted table in RMarkdown
knitr::kable(summary_table, caption = bquote("Differential Expression Summary: *nup211-so* (16h) vs. wild-type"))

# Optional: Use flextable for enhanced formatting
ft <- flextable::flextable(summary_table) %>%
  flextable::theme_box() %>%
  flextable::set_caption("S1 Table: Differential Expression Summary: nup211-so (16h) vs. wild-type",
                         style = fp_text(font.size = 12, font.family = "Arial", italic = FALSE),
                         align_with_table = FALSE) %>%
  flextable::font(part = "header", fontname = "Arial") %>%
  flextable::font(part = "body", fontname = "Arial") %>%
  flextable::bg(bg = "#f0f0f0", part = "header") # or this for greenish color: #d4edda  or #f0f0f0 for classic gray

doc <- read_docx() %>%
  body_add_flextable(ft) %>%
  body_add_par("Generated using R.", style = "Normal")

# Save the document
# print(doc, target = "~/Rstudio/expression_analysis_16h/results/DE_spreadsheets/summary_table.docx")


