###############################################################
# RNA-seq Analysis Pipeline (48h) - Script 2b
# Comparing nup211-so (48h) vs. non-shutoff cells (no thiamine)
# Purpose: Creating Plots / Data Visualization
#################################################
# Setup and Library Loading
###########################################
# Load required libraries
libraries <- c(
  "readr", "dplyr", "tidyr", "tidyverse", "tximport",
  "GenomicFeatures", "AnnotationDbi", "AnnotationForge",
  "GenomicAlignments", "BiocParallel", "DESeq2",
  "apeglm", "ggplot2", "plotly", "gplots", "Glimma",
  "pheatmap", "RColorBrewer", "viridis",
  "EnhancedVolcano", "PoiClaClu", "ggpubr"
)

for (lib in libraries) {
  if (!require(lib, character.only = TRUE)) {
    BiocManager::install(lib)
    library(lib, character.only = TRUE)
  }
}
############################################################################
# 0.1.a Read data in (if necessary)
results_annotated_48h <- read_csv("~/Desktop/nup211_GEO_sumbission_20221014/processed data files/DEseq2Results.csv") %>%
  column_to_rownames(var = names(.)[1])

filtered_annotated_48h <- results_annotated_48h %>% filter(Type == "protein coding gene") %>%
  rename(gene_type = Type)

# 0.1.a Read data in (if necessary) -- Salmon version
# results_annotated_48h_salmon <- read_csv("~/Rstudio/expression_analysis_48h/salmon_pipeline/results/DE_spreadsheets/DESeq2_results_all_transcripts_salmon48h.csv") %>% 
#   mutate(Gene_ID = `...1`) %>% select(-c(1))
# filtered_annotated_48h_salmon <- read_csv("~/Rstudio/expression_analysis_48h/salmon_pipeline/results/DE_spreadsheets/DESeq2_results_protein_coding_salmon48h.csv")%>% 
#   mutate(Gene_ID = `...1`) %>% select(-c(1))


######### Volcano Plot function ###########

create_48h_volcano <- function(results_df, highlighted_genes = NULL, title = NULL, selectLab = "", subtitle =  NULL) {
  
  keyvals <- ifelse(
    results_df$log2FoldChange < -1 & results_df$padj < 10e-6, 'royalblue',
    ifelse(results_df$log2FoldChange >  1 & results_df$padj < 10e-6, 'red2',
           'grey30'))
  keyvals[is.na(keyvals)] <- 'grey30'
  names(keyvals)[keyvals == 'red2']      <- 'Up-regulated'
  names(keyvals)[keyvals == 'royalblue'] <- 'Down-regulated'
  names(keyvals)[keyvals == 'grey30']    <- 'NS'
  
  plot <- EnhancedVolcano(
    results_df,
    lab = results_df$gene_name,
    x = 'log2FoldChange',
    y = 'padj',
    #ylim = c(0, max(-log10(results_df[["padj"]]), na.rm = TRUE) + 0.5),
    pCutoff = 10e-6,
    FCcutoff = 1,
    pointSize = 3.0,
    labSize = 6.0,
    labFace = "bold",
    legendLabels= my_legend_labels,
    colCustom = keyvals,
    legendPosition = 'right',
    legendLabSize = 16,
    legendIconSize = 7,
    selectLab = selectLab,
    title = title,
    titleLabSize = 20,
    boxedLabels = T,
    colAlpha = 3/5,
    vline = c(-2, -1, 1, 2),
    vlineCol = 'black',
    vlineType = c('twodash', 'longdash', 'longdash', 'twodash'),
    vlineWidth = 0.4,
    drawConnectors = TRUE,
    gridlines.major = T,
    gridlines.minor = T ,
    maxoverlapsConnectors = 20,  # Added this parameter
    widthConnectors = 0.5,  # Optional: adjust connector line width
    subtitle = subtitle,
    subtitleLabSize = 16,
    captionLabSize = 16
  )
  
  return(plot +
           coord_cartesian(xlim=c(-6, 7)) +
           ggplot2::scale_x_continuous(
             breaks=seq(-6, 6, 2)) +
           ggplot2::scale_y_continuous(limits = c(0,300),
                                       breaks = seq(0, 300, 50)) +
           xlab(bquote(bold(log[2]~ 'Fold Change'))) + ylab(bquote(bold(-log[10]~ '(P)'))) +
           theme(legend.text = element_text(family= 'Arial', face = 'bold', size = 16),
                 axis.title.y = element_text(family= 'Arial', face = 'bold', size = 20),
                 axis.title.x = element_text(family= 'Arial', face = 'bold', size = 20),
                 #axis.title = element_text(family= 'Arial', face="bold"),
                 axis.text.x=element_text(family= 'Arial', size = 20, face= "plain",colour = 'black'),
                 axis.text.y=element_text(family= 'Arial', size = 20, face= "plain",colour = 'black'),
                 axis.line = element_line(size = 1.0), # Control axis thickness
                 axis.ticks = element_line(size = 1.2), # Control tick mark thickness
                 axis.ticks.length = unit(3, "mm")))
}

top_genes <- results_annotated %>%
  filter(gene_type == "protein coding gene" & pvalue < 10e-6) %>%  # Filter for protein coding genes
  mutate(
    # Create absolute fold change for easier ranking
    abs_log2FC = abs(log2FoldChange),
    # Identify direction
    direction = ifelse(log2FoldChange > 0, "up", "down")
  ) %>%
  group_by(direction) %>%
  # Rank within each direction by absolute fold change and p-value
  arrange(pvalue, desc(abs_log2FC)) %>%
  slice_head(n = 10) %>%
  pull(gene_name)  # Or whatever your gene name column is called

rtqpcr_genes <- c("chp1", "mit1", "mus7", "for3", "tea2", "pho1", "urg2")


mRNA_volcano_48h_blank <- create_48h_volcano(filtered_annotated_48h, title = "Differential expression analysis - protein coding (48h)",
                                              selectLab = "", subtitle = bquote(italic("nup211-so") ~ "(+T) vs. control"))

mRNA_volcano_48h_blank <- mRNA_volcano_48h_blank+ scale_color_manual(
                        values = c('Up-regulated' = 'red2',
                                   'Down-regulated' = 'royalblue',
                                   'NS' = 'grey30'),
                        breaks = c('Up-regulated', 'Down-regulated', 'NS')
                      )

mRNA_volcano_48h_blank
#1080x735

# ggsave(mRNA_volcano_16h_blank,file=paste0(outputFolder, "volcano_protein_coding.tiff"), width = 11, height = 7, units = "in", dpi = 300)

###########################################
# Results Summary Table 
###########################################

# Get summary data (if not done previously)
# Function to summarize DE results into a dataframe
# summarize_DE <- function(data, label, padj_cutoff) {
#   data.frame(
#     Dataset = label,
#     padj_Cutoff = padj_cutoff,
#     `2-Fold Up` = sum(data$padj <= padj_cutoff & data$log2FoldChange >= 1),
#     `2-Fold Down` = sum(data$padj <= padj_cutoff & data$log2FoldChange <= -1),
#     `4-Fold Up` = sum(data$padj <= padj_cutoff & data$log2FoldChange >= 2),
#     `4-Fold Down` = sum(data$padj <= padj_cutoff & data$log2FoldChange <= -2),
#     `2-Fold DE` = sum(data$padj <= padj_cutoff & abs(data$log2FoldChange) >= 1),
#     `4-Fold DE` = sum(data$padj <= padj_cutoff & abs(data$log2FoldChange) >= 2)
#   )
# }
# 
# summary_table_48h <- bind_rows(
#   summarize_DE(results_annotated_48h, "All Transcripts", 0.05),
#   summarize_DE(results_annotated_48h, "All Transcripts", 10e-6),
#   summarize_DE(filtered_annotated_48h, "Coding Genes Only", 0.05),
#   summarize_DE(filtered_annotated_48h, "Coding Genes Only", 10e-6)
# )
# 
# summary_table_48h <- summary_table_48h %>%
#   rename( `p.adjusted cutoff` = padj_Cutoff,
#           `2-Fold Up` = `X2.Fold.Up`,
#           `2-Fold Down` = `X2.Fold.Down`,
#           `4-Fold Up` = `X4.Fold.Up`,
#           `4-Fold Down` = `X4.Fold.Down`,
#           `2-Fold DE` = `X2.Fold.DE`,
#           `4-Fold DE` = `X4.Fold.DE`
#   )
# 
# # Display as a nicely formatted table in RMarkdown
# knitr::kable(summary_table_48h, caption = bquote("Differential Expression Summary: *nup211-so* (48h) vs. control"))
# 
# # Optional: Use flextable for enhanced formatting
# library(tidyverse)
# library(flextable)
# library(officer)
# ft2 <- flextable::flextable(summary_table_48h) %>%
#   flextable::theme_box() %>%
#   set_caption("S1 Table: Differential Expression Summary: nup211-so (48h) vs. control",
#               style = fp_text(font.size = 12, font.family = "Arial", italic = F),
#               align_with_table = FALSE) %>%
#   font(part = "header", fontname = "Arial") %>%
#   font(part = "body", fontname = "Arial")  %>%
#   bg(bg = "#f0f0f0", part = "header") #%>%  # or this for greenish color: #d4edda  or #f0f0f0 for classic gray
# #flextable::autofit()
# 
# doc <- read_docx() %>%
#   body_add_flextable(ft2) %>%
#   body_add_par("Generated using R.", style = "Normal")
# 
# # Save the document
# print(doc, target = "summary_table_48h.docx")