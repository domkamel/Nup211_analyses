#################################################
# RNA-seq Analysis (48h) - nup211-so vs. wild-type
#
# Processed data are loaded from GEO GSE216106 (Sookdeo, Kamel et al.), so no
# data reported here is re-derived. The upstream alignment and counting
# steps that produced the deposited data are preserved, commented out, in the
# provenance section below.
#################################################

library(tidyverse)
library(DESeq2)
library(GEOquery)
library(pheatmap)
library(RColorBrewer)
library(flextable)
library(officer)

###########################################
# Preserved Code - how the deposited data were generated ---------------------------------
###########################################
# Retained verbatim as a record. Not run; requires the BAM files and
# sample_table.csv on the original OneDrive path or from GEO. Package installs and library
# calls are omitted here since they are loaded live above.
#
# # Set the working directory  -----------------------------------------------------------
# dir <- "/Users/domkamel/Library/CloudStorage/OneDrive-CUNY/Bam_files"
# list.files(dir)
#
# # Set the prefix for each output file name
# outputPrefix <- "shutoff_48h"
#
# csvfile <- file.path(dir,"sample_table.csv")
# sampleTable <- read.csv(csvfile,row.names=1)
# sampleTable
#
# filenames <- file.path(dir, paste0(sampleTable$Run, "Aligned.sortedByCoord.out.bam"))
# file.exists(filenames)
#
# library("Rsamtools")
# bamfiles <- BamFileList(filenames, yieldSize=2000000)
#
# seqinfo(bamfiles[1])
#
# # Make TxDb TxDbFromBiomart to generate features  -----------------------------------------------------------
#
# txdb <-  makeTxDbFromBiomart(biomart ="fungi_mart" ,dataset="spombe_eg_gene" ,host="https://fungi.ensembl.org")
# ebg <- exonsBy(txdb, by="gene")
# ebg
#
# # Summarize Overlaps -------------------------------------------------------------------------------------
#
# library("GenomicAlignments")
# library("BiocParallel")
#
# register(MulticoreParam(2))
#
# # Run once:
#   # se <-summarizeOverlaps(features=ebg,
#   #                        reads=bamfiles,
#   #                        mode="Union",
#   #                        singleEnd=FALSE,
#   #                        ignore.strand=FALSE,
#   #                        fragments=TRUE,
#   #                        preprocess.reads = invertStrand)
#
#
# # Read in old RDS of saved SE ---------------------------------------------
#
# se <- readRDS("/Users/domkamel/Library/CloudStorage/OneDrive-CUNY/original_summarizedExperiment.rds")
#
# # Quality assesment of the SE object --------------------------------------
# se
# head(assay(se))
# dim(se)
# nrow(se)
# ncol(se)
# rowRanges(se)
# length(rowRanges(se))
# rowRanges(se)[[1]]
# str(metadata(rowRanges(se)))
# colData(se)
# colnames(se)
# bamfiles
# sampleTable$Run
#
# colData(se) <- DataFrame(sampleTable)
# colData(se)
#
# # Create DESeq dataset ---------------------------------------------------------------------------------
#
# colData(se)
#
# # Note: it will be helpful if the first level of a factor is the reference level (e.g. control, or untreated samples).
#
# se$thiamine <- relevel(se$thiamine, "untrt")
# se$thiamine
#
# # DESeqDataSet starting from SummarizedExperiment -------------------------------------------------------
#
# # Pre-filtering rows with very small counts
#
# se <- se[ rowSums(assay(se)) >= 5, ]
#
# se$thiamine <- relevel(se$thiamine, "untrt")
# se$thiamine
#
# library("DESeq2")
# dds <- DESeqDataSet(se, design = ~ thiamine)
#
#
# # Getting count matricies --------------------------------------------------------------------------------
#
# countdata <- assay(se)
# head(countdata, 3)
#
# countdata.df <- as.data.frame(countdata)
# coldata <- colData(se)
#
# # dds <- DESeq(dds)
# # res <- results(dds)
# #
# # Results deposited as GSE216106_DEseq2Results.csv.gz
# # Counts deposited as  GSE216106_counts.csv.gz


###########################################
# Download processed data from GEO
###########################################
# Create the geo directory on your machine (example below)
geo_dir <- path.expand("~/Rstudio/expression_analysis_48h/geo")
dir.create(geo_dir, recursive = TRUE, showWarnings = FALSE)

library(GEOquery)
# Downloads into <geo_dir>/GSE216106/. Skipped if already present.
if (!dir.exists(file.path(geo_dir, "GSE216106"))) {
  getGEOSuppFiles("GSE216106", baseDir = geo_dir)
}
list.files(file.path(geo_dir, "GSE216106"))

library(dplyr)
counts_48h <- read.csv(gzfile(file.path(geo_dir, "GSE216106/GSE216106_counts.csv.gz")),
                       row.names = 1) %>% as.matrix()

res_48h <- read.csv(gzfile(file.path(geo_dir, "GSE216106/GSE216106_DEseq2Results.csv.gz")))

colnames(counts_48h)   # samples 1-3 = control, 4-6 = thiamine treated
colnames(res_48h)      # already annotated - confirm the gene ID column name
head(res_48h, 3)

# First column holds the gene ID; rename for consistency with the 16h script
names(res_48h)[1] <- "Gene_ID"

###########################################
# Optional: refresh annotation against a current PomBase export
###########################################
# The deposited table was annotated at submission (Oct 2022). This swaps those
# columns for current gene names/products and leaves every statistic untouched.
# Set FALSE to keep the table exactly as deposited.

UPDATE_ANNOTATION <- TRUE

if (UPDATE_ANNOTATION) {
  geneIDs_descriptions <- read.delim("~/Downloads/gene_IDs_names_products (1).tsv") %>%
    mutate_if(is.character, list(~na_if(., ""))) %>%
    mutate(gene_name = coalesce(gene_name, gene_systematic_id))

  stale <- setdiff(intersect(names(res_48h), names(geneIDs_descriptions)), "Gene_ID")
  stale   # annotation columns being replaced

  res_48h <- res_48h %>%
    select(-all_of(stale)) %>%
    left_join(geneIDs_descriptions, by = c("Gene_ID" = "gene_systematic_id"))

  sum(is.na(res_48h$gene_name))   # unmatched IDs; should be 0 or close to it
}

res_48h <- res_48h %>% arrange(log2FoldChange)

# write.csv(res_48h, "~/Rstudio/expression_analysis_48h/DESeq2_results_all_transcripts.csv", row.names = FALSE)

###########################################
# Sample metadata
###########################################
coldata <- data.frame(
  row.names = colnames(counts_48h),
  condition = factor(rep(c("wild_type", "nup211so_48h"), each = 3),
                     levels = c("wild_type", "nup211so_48h"))
)

# Used only for the QC transformations below, not for any DE statistic
library(DESeq2)
dds_48h <- DESeqDataSetFromMatrix(counts_48h, coldata, ~ condition)
dds_48h <- estimateSizeFactors(dds_48h[rowSums(counts(dds_48h)) > 10, ])

###########################################
# QC plots
###########################################
vsd <- vst(dds_48h)

cols <- c("wild_type" = "#52d8da", "nup211so_48h" = "#F8766D")
labs <- c("wild-type", expression(italic("nup211-so") ~ "(48h)"))

# PCA
pca <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pv  <- round(100 * attr(pca, "percentVar"))

library(ggplot2)
ggplot(pca, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) + theme_bw() +
  labs(x = paste0("PC1: ", pv[1], "% variance"),
       y = paste0("PC2: ", pv[2], "% variance"),
       title = bquote(bold("PCA:" ~ bolditalic("nup211-so") ~ "(48h) vs. wild-type"))) +
  scale_color_manual(name = "Condition", values = cols, labels = labs)

# Sample distances
sampleDistMatrix <- as.matrix(dist(t(assay(vsd))))
rownames(sampleDistMatrix) <- paste(vsd$condition)
colnames(sampleDistMatrix) <- NULL

library(pheatmap)
library(RColorBrewer)
pheatmap(sampleDistMatrix,
         color = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         main = "Sample Distance Matrix\n")

# Normalized counts for a gene
plot_gene <- function(id, label = id) {
  ggplot(plotCounts(dds_48h, id, "condition", returnData = TRUE),
         aes(condition, count, color = condition)) +
    geom_point(position = position_jitter(w = 0.1), size = 3) +
    scale_y_log10() + scale_color_manual(values = cols) + theme_bw() +
    labs(title = label, x = NULL, y = "normalized count") +
    theme(legend.position = "none", plot.title = element_text(face = "bold.italic"))
}

plot_gene("SPCC162.08c", "nup211")
plot_gene("SPNCRNA.1206", "nup211-antisense-1")

# MA plot
res_48h %>%
  filter(!is.na(padj)) %>%
  mutate(sig = padj <= 0.05 & abs(log2FoldChange) >= 1) %>%
  ggplot(aes(baseMean, log2FoldChange, color = sig)) +
  geom_point(size = 0.7, alpha = 0.6) +
  scale_x_log10() + ylim(-8, 8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#F8766D"), guide = "none") +
  theme_bw() +
  labs(x = "mean of normalized counts", y = "log2 fold change",
       title = bquote(bold("MA:" ~ bolditalic("nup211-so") ~ "(48h) vs. wild-type")))

###########################################
# Summary table
###########################################
filtered_annotated <- res_48h %>% filter(gene_type == "protein coding gene")

summarize_DE <- function(data, label, cutoff) {
  keep <- !is.na(data$padj) & data$padj <= cutoff
  data.frame(
    Dataset = label,
    `p.adjusted cutoff` = cutoff,
    `2-Fold Up`   = sum(keep & data$log2FoldChange >= 1),
    `2-Fold Down` = sum(keep & data$log2FoldChange <= -1),
    `4-Fold Up`   = sum(keep & data$log2FoldChange >= 2),
    `4-Fold Down` = sum(keep & data$log2FoldChange <= -2),
    `2-Fold DE`   = sum(keep & abs(data$log2FoldChange) >= 1),
    `4-Fold DE`   = sum(keep & abs(data$log2FoldChange) >= 2),
    check.names = FALSE
  )
}

summary_table <- bind_rows(
  summarize_DE(res_48h,            "All Transcripts",   0.05),
  summarize_DE(res_48h,            "All Transcripts",   10e-6),
  summarize_DE(filtered_annotated, "Coding Genes Only", 0.05),
  summarize_DE(filtered_annotated, "Coding Genes Only", 10e-6)
)

summary_table

ft <- flextable(summary_table) %>%
  theme_box() %>%
  set_caption("Differential Expression Summary: nup211-so (48h) vs. wild-type",
              style = fp_text(font.size = 12, font.family = "Arial"),
              align_with_table = FALSE) %>%
  font(part = "all", fontname = "Arial") %>%
  bg(bg = "#f0f0f0", part = "header")

# print(read_docx() %>% body_add_flextable(ft),
#       target = "~/Rstudio/expression_analysis_48h/summary_table.docx")
