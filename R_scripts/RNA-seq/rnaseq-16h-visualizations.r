#################################################
# RNA-seq Analysis Pipeline (16h) - Script 2
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

###########################################
# Data Import
###########################################

# skip if you've just run DESeq2-16h-expression.r

results_annotated <- read_csv("~/Rstudio/expression_analysis_16h/results/DE_spreadsheets/DESeq2_results_all_transcripts.csv") %>%
                     column_to_rownames(var = names(.)[1])

filtered_annotated <- read_csv("~/Rstudio/expression_analysis_16h/results/DE_spreadsheets/DESeq2_results_protein_coding.csv") %>%
                      column_to_rownames(var = names(.)[1])

###########################################
# Visualization Functions
###########################################
outputFolder <- "~/Rstudio/expression_analysis_16h/results/visualizations/volcano_plots/"

######### Volcano Plot function ###########

create_16h_volcano <- function(results_df, highlighted_genes = NULL, title = NULL, selectLab = "", subtitle =  NULL) {
  
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
    ylim = c(0, max(-log10(results_df[["padj"]]), na.rm = TRUE) + 0.5),
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
           coord_cartesian(xlim=c(-4.5, 10)) +
           ggplot2::scale_x_continuous(
             breaks=seq(-6,10, 2)) +
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

rtqpcr_genes <- c("chp1", "mit1", "swi6", "for3", "tea2", "mus7", "not1")

# Make a volcano plot (no genes labeled)
mRNA_volcano_16h_blank <- create_16h_volcano(filtered_annotated, title = "16h Depletion RNA-seq",
                                             selectLab = "", subtitle = bquote(italic("nup211-so") ~ "(+T) vs. wild-type (+T)"))

# Re-level the key/legend 
mRNA_volcano_16h_blank <- mRNA_volcano_16h_blank + scale_color_manual(
  values = c('Up-regulated' = 'red2',
             'Down-regulated' = 'royalblue',
             'NS' = 'grey30'),
  breaks = c('Up-regulated', 'Down-regulated', 'NS')
)

# Label the genes analyzed by RT-qPCR
mRNA_volcano_16h_PCR <- create_16h_volcano(filtered_annotated, title = "16h Depletion RNA-seq",
                                           selectLab = rtqpcr_genes, subtitle = bquote(italic("nup211-so") ~ "(+T) vs. wild-type (+T)"))

# Then change the order of the legend labels 
mRNA_volcano_16h_PCR <- mRNA_volcano_16h_PCR + scale_color_manual(
  values = c('Up-regulated' = 'red2',
             'Down-regulated' = 'royalblue',
             'NS' = 'grey30'),
  breaks = c('Up-regulated', 'Down-regulated', 'NS')
)


mRNA_volcano_16h_blank
mRNA_volcano_16h_PCR
# 1080x735 could be good dimensions to use for saving

# ggsave(mRNA_volcano_16h_blank,file=paste0(outputFolder, "volcano_protein_coding.tiff"), width = 11, height = 7, units = "in", dpi = 300)
# ggsave(mRNA_volcano_16h_PCR,file=paste0(outputFolder, "volcano_coding_RTqPCR_labels.tiff"), width = 11, height = 7, units = "in", dpi = 300)


######### For CUNYSciCom ################## 

# Simplified plot for a more general audience 
create_16h_volcano_simple <- function(results_df, title = NULL, subtitle = NULL) {
  
  keyvals <- ifelse(
    results_df$log2FoldChange < -1 & results_df$padj < 10e-6, 'royalblue',
    ifelse(results_df$log2FoldChange >  1 & results_df$padj < 10e-6, 'red2',
           'grey70'))
  keyvals[is.na(keyvals)] <- 'grey70'
  names(keyvals)[keyvals == 'red2']      <- 'Up'
  names(keyvals)[keyvals == 'royalblue'] <- 'Down'
  names(keyvals)[keyvals == 'grey70']    <- 'NS'
  
  # Symmetric fold-change labels: ÷ on the left, 0 at center, × on the right
  fc_breaks <- -5:5
  fc_labels <- c("1/32", "1/16", "1/8", "1/4", "1/2",
                 "0",
                 "2\u00d7", "4\u00d7", "8\u00d7", "16\u00d7", "32\u00d7")
  
  plot <- EnhancedVolcano(
    results_df,
    lab = results_df$gene_name,
    x = 'log2FoldChange',
    y = 'padj',
    ylim = c(0, max(-log10(results_df[["padj"]]), na.rm = TRUE) + 0.5),
    pCutoff = 10e-6,
    FCcutoff = 1,
    pointSize = 3.0,
    selectLab = "",
    colCustom = keyvals,
    colAlpha = 3/5,
    title = title,
    titleLabSize = 20,
    subtitle = subtitle,
    subtitleLabSize = 16,
    caption = "",
    legendPosition = 'top',
    legendLabSize = 16,
    legendIconSize = 7,
    gridlines.major = TRUE,    # gridlines back on
    gridlines.minor = TRUE,
    drawConnectors = FALSE
  )
  
  return(plot +
           # Bold center reference line at "no change"
           #geom_vline(xintercept = 0, linewidth = 1.0, colour = "black") +
           coord_cartesian(xlim = c(-5, 5)) +
           ggplot2::scale_x_continuous(breaks = fc_breaks, labels = fc_labels) +
           xlab("Fold change") +
           ylab("Significance") +
           theme(
             legend.position   = "right",
             legend.title      = element_blank(),
             legend.text       = element_text(family = 'Arial', face = 'bold', size = 16),
             axis.title.x      = element_text(family = 'Arial', face = 'bold', size = 20),
             axis.title.y      = element_text(family = 'Arial', face = 'bold', size = 20),
             axis.text.x       = element_text(family = 'Arial', size = 18, colour = 'black'),
             axis.text.y       = element_blank(),
             axis.ticks.y      = element_blank(),
             axis.ticks.x      = element_line(size = 1.2),
             axis.ticks.length = unit(3, "mm"),
             # Box around the whole panel
             panel.border      = element_rect(colour = "black", fill = NA, linewidth = 2.0),
             axis.line         = element_blank()  # disable separate axis lines so the box is clean
           ))
}

mRNA_volcano_16h_simple <- create_16h_volcano_simple(filtered_annotated, title = "16h Depletion RNA-seq",
                                                     subtitle = bquote(italic("nup211-so") ~ "(+T) vs. wild-type (+T)"))
#1080x735


mRNA_volcano_16h_simple
# ggsave(plot = mRNA_volcano_16h_simple,
#        filename = "simple_volcano_16h-2.png",
#        width = 7, height = 6.8, units = "in", dpi = 600, device = cairo_pdf)
# 
# ggsave(plot = mRNA_volcano_16h_blank,
#        filename = "SciCom_volcano.png",
#        width = 9.5, height = 7.5, units = "in", dpi = 600, device = cairo_pdf)



########### Heatmap functions ############

create_heatmap <- function(vsd, res, results_annotated, top_n, filter_coding = FALSE) {
  # Get initial expression matrix and center by row means
  mat <- assay(vsd)[order(res$padj), ]
  mat <- mat - rowMeans(mat)
  
  # Standardize C/c in gene names
  rownames(mat) <- sapply(rownames(mat), function(x) {
    if (substr(x, nchar(x), nchar(x)) == "C") {
      paste0(substr(x, 1, nchar(x) - 1), "c")
    } else {
      x
    }
  })
  
  # Filter and merge with annotations
  if (filter_coding) {
    # Filter for protein coding genes only
    coding_genes <- results_annotated %>%
      filter(gene_type == "protein coding gene")
    
    mat <- merge(mat, coding_genes, by = 0)
  } else {
    # Merge with full annotations
    mat <- merge(mat, results_annotated, by = 0)
  }
  
  # Set row names and select only expression columns
  rownames(mat) <- mat$gene_name
  mat <- mat[order(res$padj), ] # Maintain ordering by p-value
  mat <- mat[!is.na(mat), ]     # Remove NA values
  mat <- mat[, c(2:7)]          # Keep only expression columns
  mat <- head(mat, top_n)       # Select top n genes
  
  # Rename columns for clarity
  names(mat) <- c("WT1", "WT2", "WT3", "16h1", "16h2", "16h3")
  
  # Create annotation dataframe
  df <- data.frame(
    condition = factor(c(rep("wild_type", 3), rep("16h_shutoff", 3)),
                       levels = c("wild_type", "16h_shutoff")),
    row.names = c("WT1", "WT2", "WT3", "16h1", "16h2", "16h3")
  )
  
  # Create and return the heatmap
  pheatmap(mat, 
           annotation_col = df,
           show_rownames = TRUE,
           cluster_rows = T,
           cluster_cols = TRUE,
           annotation_names_col = TRUE,
           annotation_colors = list(
             condition = c(wild_type = "white", `16h_shutoff` = "black")
           ),
           angle_col = 0,
           treeheight_col = 20,
           treeheight_row = 40,
           scale = "none",
           color = colorRampPalette(rev(brewer.pal(n = 256, name = "RdBu")))(100))
}

unfiltered_heatmap <- create_heatmap(vsd = vsd, res = resNorm,results_annotated = results_annotated, 
                                     top_n = 50, filter_coding = F) # change to true for filtering out ncRNAs

filtered_heatmap <- create_heatmap(vsd = vsd, res = resNorm,results_annotated = results_annotated, 
                                   top_n = 50, filter_coding = T) 

outputFolder <- "~/Rstudio/expression_analysis_16h/results/visualizations/heatmaps/"

# ggsave(unfiltered_heatmap,file=paste0(outputFolder, "expression_heatmap_all.tiff"), width = 7, height = 8, units = "in", dpi = 300)
# ggsave(filtered_heatmap,file=paste0(outputFolder, "expression_heatmap_protein_coding.tiff"), width = 7, height = 8, units = "in", dpi = 300)


###########################################
# Glimma Plot
###########################################

library(Glimma)

create_glimma_volcano <- function(data, index) {
  # Create status vector for coloring
  status <- ifelse(
    data$padj < 10e-6 & abs(data$log2FoldChange) > 1,
    "DE", "Not DE"
  )
  
  # Create PomBase links
  pombase_links <- paste0(
    '<a href="https://www.pombase.org/gene/', 
    rownames(data), 
    '" target="_blank">', 
    data$gene_name, 
    '</a>'
  )
  
  # Create the annotation data frame
  anno <- data.frame(
    Gene = pombase_links,
    GeneID = data$Gene_ID,
    Name = data$gene_name,
    LogFC = round(data$log2FoldChange, 2),
    AdjPValue = signif(data$padj, 3),
    Description = data$gene_product,
    BaseMean = signif(data$baseMean, 3),
    stringsAsFactors = FALSE
  )
  
  # Display columns in the order you want
  display_order <- c("Gene", "Description", "log2 Fold Change","neglog10(padj)","BaseMean")
  
  # Create the Glimma plot
  glXYPlot(
    x = data$log2FoldChange,
    y = -log10(data$padj),
    xlab = "log2 Fold Change",
    ylab = "neglog10(padj)",
    main = "Interactive Volcano Plot",
    anno = anno,
    groups = status,
    #status = status,
    display.columns = display_order,
    side.main = "Gene",  # Changed to match one of the display.columns
    status = status,
    cols = c("gray", "firebrick"),
    folder = "/Rstudio/expression_analysis_16h/results/visualizations/Glimma_plots/",
    html = index,
    launch = TRUE
  )
}

# create_glimma_volcano(results_annotated, index = "all_transcripts")
# create_glimma_volcano(filtered_annotated, index = "protein_coding")

###########################################
#            MA Plot function              -- use after running "DESeq2-16h-expression.r" script
###########################################
create_ma_plot <- function(res, title) {
  DESeq2::plotMA(res, ylim = c(-10, 10), main = title, colSig = "#F8766D", colNonSig = "gray",
                 cex = 0.7,      # Increase point size (default is 1)
                 ylab = bquote(bold(log[2]~ 'Fold Change')),
                 xlab = bquote(bold('mean of normalized counts')))
}
resLFC1 <- lfcShrink(dds_16h, coef = "condition_nup211so_16h_vs_wild_type", 
                     type = "normal", lfcThreshold = 1) 

create_ma_plot(resLFC1, bquote(bold("MA Plot:" ~ bolditalic("nup211-so") ~ "(16h) vs. wild-type")))



###########################################
# MA Glimma Plot
###########################################

create_glimma_MA <- function(data, index) {
  # Create status vector for coloring
  status <- ifelse(
    data$padj < 10e-6 & abs(data$log2FoldChange) > 1,
    "DE", "Not DE"
  )
  
  # Create PomBase links
  pombase_links <- paste0(
    '<a href="https://www.pombase.org/gene/', 
    rownames(data), 
    '" target="_blank">', 
    data$gene_name, 
    '</a>'
  )
  
  # Create the annotation data frame
  anno <- data.frame(
    Gene = pombase_links,
    GeneID = data$Gene_ID,
    Name = data$gene_name,
    LogFC = round(data$log2FoldChange, 2),
    AdjPValue = signif(data$padj, 3),
    Description = data$gene_product,
    BaseMean = signif(data$baseMean, 3),
    stringsAsFactors = FALSE
  )
  
  # Display columns in the order you want
  display_order <- c("Gene", "Description", "log2 Fold Change","neglog10(padj)","BaseMean")
  
  # Create the Glimma plot
  glXYPlot(
    x = log10(data$baseMean + 1),
    y = data$log2FoldChange,
    ylab = "log2 Fold Change",
    xlab = "log10(baseMean+1)",
    main = "Interactive MA Plot",
    anno = anno,
    groups = status,
    #status = status,
    display.columns = display_order,
    side.main = "Gene",  # Changed to match one of the display.columns
    status = status,
    cols = c("gray", "firebrick"),
    folder = "/Rstudio/expression_analysis_16h/results/visualizations/Glimma_plots/",
    html = index,
    launch = TRUE
  )
}

# create_glimma_MA(results_annotated, index = "all_transcripts_MA")
# create_glimma_MA(filtered_annotated, index = "protein_coding_MA")
