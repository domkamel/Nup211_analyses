#####################################################
# RNA-seq Analysis Pipeline (16h) - Script 5
# Purpose: Compare RNA-seq data for 16h and 48h
#
#     *** 48h Using STAR alignment ***
#             Published Dataset
#####################################################

# Clear workspace
# rm(list=ls())

###########################################
# Setup and Library Loading
###########################################
# Load required libraries
libraries <- c(
  "readr", "dplyr", "tidyr", "tidyverse", "tximport",
  "GenomicFeatures", "AnnotationDbi", "AnnotationForge",
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

outputFolder <- "~/Rstudio/expression_analysis_16h/results/visualizations/correlation/"

###########################################
# Compare 16h and 48h data 
###########################################

shutoff_48h_all <- read_csv("~/Desktop/nup211_GEO_sumbission_20221014/processed data files/DEseq2Results.csv")
shutoff_48h_filtered <- shutoff_48h_all %>% filter(Type == "protein coding gene") %>%
  rename(gene_type = Type)

shutoff_16h_all <- read_csv("~/Rstudio/expression_analysis_16h/results/DE_spreadsheets/DESeq2_results_all_transcripts.csv") %>% 
  mutate(Gene_ID = `...1`) %>% select(-c(1))
shutoff_16h_filtered <- read_csv("~/Rstudio/expression_analysis_16h/results/DE_spreadsheets/DESeq2_results_protein_coding.csv")%>% 
  mutate(Gene_ID = `...1`) %>% select(-c(1))


###############################
# Set cutoff values
###############################

log2fc_threshold = 1.0
padj_threshold = 10e-6

###############################
# Clean and merge datasets
###############################
  
  merged_seq_data <- shutoff_16h_filtered %>%
    select(Gene_ID, gene_name, gene_product, baseMean, log2FoldChange, padj) %>%
    rename(log2fc_16h = log2FoldChange, baseMean_16h = baseMean, padj_16h = padj) %>%
    filter(!is.infinite(padj_16h), !is.na(padj_16h)) %>%
    inner_join(
      shutoff_48h_filtered %>%
        select(Gene_ID, baseMean, log2FoldChange, padj) %>%
        rename(log2fc_48h = log2FoldChange, baseMean_48h = baseMean, padj_48h = padj) %>%
        filter(!is.infinite(padj_48h), !is.na(padj_48h)),
      by = "Gene_ID"
    ) %>%
    filter(!is.infinite(log2fc_16h), !is.infinite(log2fc_48h),
           !is.na(log2fc_16h), !is.na(log2fc_48h))

###############################
# Add significance status
###############################

  merged_seq_data <- merged_seq_data %>%
      mutate(
        is_significant = (
          abs(log2fc_16h) > log2fc_threshold &
            abs(log2fc_48h) > log2fc_threshold &
            padj_16h < padj_threshold &
            padj_48h < padj_threshold &
            sign(log2fc_16h) == sign(log2fc_48h) # Check if signs match
        ),  
        sig_16h = (abs(log2fc_16h) > log2fc_threshold &
                    padj_16h < padj_threshold), 
        sig_48h = (abs(log2fc_48h) > log2fc_threshold &
                   padj_48h < padj_threshold), 
        sig_category = case_when(
          sig_16h & sig_48h & sign(log2fc_16h) == sign(log2fc_48h) ~ "same_direction",
          sig_16h & sig_48h & sign(log2fc_16h) != sign(log2fc_48h) ~ "opposite_direction",
          sig_16h & !sig_48h ~ "16h_only",
          !sig_16h & sig_48h ~ "48h_only",
          TRUE ~ "not_significant"
        ),
        thresholds_met = (
          abs(log2fc_16h) > log2fc_threshold &
            padj_16h < padj_threshold |
          abs(log2fc_48h) > log2fc_threshold &
            padj_48h < padj_threshold
    )
  )

# Differentially expressed coding transcripts (meeting the log2 fold change and significance cutoffs and represented in the volcano plots shown in figure 1) at 
# the two time points were compared to see similarities/differences in the transcriptome profiles at the two time points. Pearson's correlation was performed on the log2 fold change values. 
# The r value (~0.20) shows ____. Each point represents a coding transcript. Points colored red are transcripts that are differentially expressed in the same direction (up or down) at both time points.

nrow(merged_seq_data)
# The number of transcripts present in both analyses is 5005


merged_seq_data %>% filter(thresholds_met == T) %>% nrow()
# The number of transcripts affected at either time-point is 1973

merged_seq_data_2 <- merged_seq_data %>% filter(thresholds_met == T)


merged_seq_data %>% filter(is_significant == TRUE) %>% nrow()
# The number of transcripts that are differentially expressed at either time-point
# in the same direction is 195

genes_same_dir <- merged_seq_data %>% filter(sig_category == "same_direction") %>% pull(gene_name)
genes_opposite <- merged_seq_data %>% filter(sig_category == "opposite_direction") %>% pull(gene_name)
genes_16h_only <- merged_seq_data %>% filter(sig_category == "16h_only") %>% pull(gene_name)
genes_48h_only <- merged_seq_data %>% filter(sig_category == "48h_only") %>% pull(gene_name)

genes_16h <- merged_seq_data %>% filter(sig_16h == TRUE) %>% pull(gene_name)
genes_48h <- merged_seq_data %>% filter(sig_48h == TRUE) %>% pull(gene_name)

genes_16h %>% length()
genes_48h %>% length()

###############################
# Optional: Get labels 
###############################
  label_data <- merged_seq_data %>% #filter(gene_name %in% custom_labels) %>%
    filter(is_significant) %>%
    slice_max(abs(log2fc_16h), n = 20)  # Adjust 'n' to control the number of labeled points

top_up <- merged_seq_data %>%
  filter(is_significant) %>%
  slice_max(log2fc_16h, n = 15)

top_down <- merged_seq_data %>%
  filter(is_significant) %>%
  slice_min(log2fc_16h, n = 10)

top_10_up_down <- bind_rows(top_up, top_down)

# write.csv(label_data,file=paste0(outputFolder,"labeled_genes.csv"))
# write.csv(merged_seq_data,file=paste0(outputFolder,"correlation_df.csv"))
###############################
# Calculate correlation
###############################

  cor_result <- cor.test(merged_seq_data$log2fc_16h, merged_seq_data$log2fc_48h)

# This is pearson correlation

##################################
# Get significant correlated genes
##################################

significant <- merged_seq_data %>%
  filter(is_significant) %>%
  mutate(across(starts_with("baseMean"), ~ formatC(., format = "f", digits = 2))) # %>% view()

significant %>% nrow()
# The number of transcripts that are differentially expressed at either timepoint
# is 195

###############################
# Filtered count
###############################

## this doesn't quite make sense so ignore for now
# filtered_count = nrow(shutoff_16h_filtered) + nrow(shutoff_48h_filtered) - nrow(merged_seq_data)

###############################
# Create correlation plot
###############################
  
  plot <- ggplot(merged_seq_data, aes(x = log2fc_16h, y = log2fc_48h)) +
    geom_abline(intercept = 0, linewidth = 1.0, slope = 1, color = "black", linetype = "dashed") +
    geom_hline(yintercept = 0, linewidth = 1.0, linetype = "solid", color = "black") +
    geom_vline(xintercept = 0, linewidth = 1.0, linetype = "solid", color = "black") +
    geom_point(aes(fill = is_significant), alpha = 0.8, shape = 21, color = "black", stroke = 0.4, size = 3) +
    # geom_label_repel(
    #   data = label_data,
    #   aes(label = gene_name),
    #   size = 5,
    #   max.overlaps = 50,  # Reduce overlaps further
    #   box.padding = 0.5,
    #   label.padding = unit(0.2, "lines")  # Adjust padding for labels
    # ) +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "lightgrey"), 
                      breaks = c(TRUE, FALSE),
                      labels = c("Yes", "No"),
                      guide = guide_legend(override.aes = list(shape = 21, size = 5))) +
    labs(
      title = paste0("RNA-seq: 48h vs 16h Results"),
      subtitle = paste0("(r = ", round(cor_result$estimate, 3), 
                        ", p=", format(cor_result$p.value, scientific = T, digits = 3), ")"),
      caption = "Coding genes; padj ≤ 10e-5",
      x = bquote(atop(bold(~log[2]~ 'Fold Change (16h)'))),
      y = bquote(atop(bold(~log[2]~ 'Fold Change (48h)'))),
      fill = "Correlated:"
    ) + 
    theme_bw()+ 
    theme( plot.title = element_text(size = 18, face = "bold"),
           plot.subtitle = element_text(size = 14, face = "bold"),
           plot.caption = element_text(size = 14, face = "plain"),
           axis.text.x = element_text(size = 16, face= "plain",colour = 'black'),    # Axis labels
           axis.text.y = element_text(size = 16, face= "plain",colour = 'black'),
           axis.title.x = element_text(size = 18, face= "bold",colour = 'black'),   # Axis titles
           axis.title.y = element_text(size = 18, face= "bold",colour = 'black'),
           legend.title = element_text(size = 14, face = "bold"),   # Adjust legend title font
           legend.text = element_text(size = 14),                  # Adjust legend labels font
           panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5)
    ) +
    scale_x_continuous(#limits = c(-10,10),
                       breaks = seq(-10, 10, by = 2)) +
    scale_y_continuous(breaks = seq(-10, 10, by = 2)) 

print(plot)

plot +
  geom_label_repel(
    data = top_10_up_down, #or rtqpcr genes or top_10_up_down or label_data
    aes(label = gene_name),
    size = 5,
    max.overlaps = 50,  # Reduce overlaps further
    box.padding = 0.5,
    label.padding = unit(0.2, "lines")  # Adjust padding for labels
  )

rtqpcr_genes <- c("chp1", "mit1", "swi6", "for3", "tea2", "mus7", "not1")
rtqpcr_labels <- merged_seq_data %>% filter(gene_name %in% rtqpcr_genes)

plot +
  geom_label_repel(
    data = rtqpcr_labels, #or rtqpcr_labelsor top_10_up_down or label_data
    aes(label = gene_name),
    size = 5,
    max.overlaps = 50,  # Reduce overlaps further
    box.padding = 0.5,
    label.padding = unit(0.2, "lines")  # Adjust padding for labels
  )
###########################################
# Venn Diagram 16h and 48h  
###########################################

# Venn Diagram 16h and 48h  using ggven
###########################################

library(ggvenn)

genes_16h_2 <- merged_seq_data_2 %>% filter(sig_16h == TRUE) %>% pull(gene_name)
genes_48h_2 <- merged_seq_data_2 %>% filter(sig_48h == TRUE) %>% pull(gene_name)

# For a simple plot limited to just three categories (e.g., 16h, 48h, overlap), simplify:

venn_list <- list(
  `16h` = genes_16h,
  `48h` = genes_48h
)

ggvenn(venn_list, fill_color = c("#F8766D", "#619CFF"), stroke_size = 0.4)


# length(genes_16h_only)    # should be ~716
# length(genes_16h)    # should be ~945
# length(genes_48h_only)    # should be ~928
# length(genes_48h)         # should be ~1157
# length(genes_same_dir)    # should be 195



#####################
# alternate approach
#####################

# Filter significant genes with direction
up_16h <- shutoff_16h_filtered %>%
  filter(log2FoldChange > 1 &
         padj < padj_threshold) %>%
  pull(gene_name)

down_16h <- shutoff_16h_filtered %>%
  filter(log2FoldChange < -1 &
           padj < padj_threshold
         ) %>%
  pull(gene_name)

up_48h <- shutoff_48h_filtered %>%
  filter(log2FoldChange > 1 &
         padj < padj_threshold
         ) %>%
  pull(gene_name)

down_48h <- shutoff_48h_filtered %>%
  filter(log2FoldChange < -1 &
         padj < padj_threshold
         ) %>%
  pull(gene_name)


library(eulerr)

# Define sets using gene vectors
venn_list <- list(
  `Up 16h` = up_16h,
  `Down 16h` = down_16h,
  `Up 48h` = up_48h,
  `Down 48h` = down_48h
)

# Create Euler diagram object
fit <- euler(venn_list)


#list(fill = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")

# Plot with options
plot(
  fit, scaling = 1,
  fills = list(fill = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a"), alpha = 0.7),
  edges = list(lwd = 1.5, col = "black"),
  quantities = list(
    cex = 1.7,               # Size of overlap numbers
    font = 1,                # Bold font
    col = "black"            # Text color
  ),
  #labels = list(font = 2, cex = 1),
  main = list(
    label = "Differentially Expression Summary",
    cex = 1.6,
    font = 2,
    col = "black"
  ),
  legend = list(cex = 1.2, font = 2, col = "black"),
  shape = "circle"  # Tries to use circles if it can
) 


summary_data <- c(
  "Up 16h total" = 661 + 163 + 28,
  "Up 16h only" = 661,
  "Up 16h&Up 48h" = 163,
  "Down 16h" = 55 + 6 +32,
  "Down 16h only" = 55,
  "Down 16h&Down 48h" = 32,
  "Up 48h" = 505 + 163 + 6,
  "Up 48h only" = 505,
  "Down 16h&Up 48h" = 6,
  "Down 48h" = 423 + 32 +28,
  "Down 48h only" = 423,
  "Up 16h&Down 48h" = 28
)

as.data.frame(summary_data)



###########################################################
#  Now, just focusing on all genes at different timepoint analyses
# (not considering if present in both datasets!)
###########################################################

genes_16h_3 <- shutoff_16h_filtered %>%
  mutate(
    is_significant = (
      abs(log2FoldChange) > log2fc_threshold &
        padj < padj_threshold)
  )  


genes_48h_3 <- shutoff_48h_filtered %>%
  mutate(
    is_significant = (
      abs(log2FoldChange) > log2fc_threshold &
        padj < padj_threshold)
  )  


genes_16h_names <-genes_16h_3 %>% filter(is_significant == TRUE) %>% pull(gene_name)

genes_48h_names <- genes_48h_3 %>% filter(is_significant == TRUE) %>% pull(gene_name)

venn_list_3 <- list(
  `16h shutoff` = genes_16h_names,
  `48h shutoff` = genes_48h_names
)

# Create Euler diagram object
fit_3 <- euler(venn_list_3)


#list(fill = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")

# Plot with options
plot(
  fit_3, scaling = 1,
  fills = list(fill = c("#1b9e77", "#d95f02"), alpha = 0.7),
  edges = list(lwd = 1.6, col = "black"),
  quantities = list(
    cex = 2,               # Size of overlap numbers
    font = 1,                # Bold font
    col = "black"            # Text color
  ),
  #labels = list(font = 2, cex = 1),
  main = list(
    label = "Number of differentially expressed transcripts",
    cex = 2,
    font = 2,
    col = "black"
  ),
  legend = list(cex = 1.6, font = 2, col = "black",
                side = "bottom", nrow = 1, ncol = 2),
  shape = "circle"  # Tries to use circles if it can
) 



geneIDs_descriptions <- read.delim("~/Downloads/gene_IDs_names_products (1).tsv") %>%
  mutate_if(is.character, list(~na_if(.,""))) %>%
  mutate(
    gene_name = coalesce(gene_name, gene_systematic_id),
    gene_systematic_id = toupper(gene_systematic_id),
    gene_systematic_id_with_prefix = gsub("PomBase:", "", gene_systematic_id_with_prefix)
  )



# Master list of all significant genes from any set
all_genes <- unique(c(up_16h, down_16h, up_48h, down_48h))

# Supplemental table: include stats + membership
supp_table <- tibble(gene_name = all_genes) %>%
  # Merge in gene_product from geneIDs_descriptions
  left_join(
    geneIDs_descriptions %>%
      select(gene_name, gene_product),
    by = "gene_name"
  ) %>%
  # Add stats from 16h
  left_join(
    shutoff_16h_filtered %>%
      select(gene_name, baseMean_16h = baseMean,log2fc_16h = log2FoldChange, padj_16h = padj),
    by = "gene_name"
  ) %>%
  # Add stats from 48h
  left_join(
    shutoff_48h_filtered %>%
      select(gene_name, baseMean_48h = baseMean,log2fc_48h = log2FoldChange, padj_48h = padj),
    by = "gene_name"
  ) %>%
  # Membership + categorization
  mutate(
    Up_16h   = gene_name %in% up_16h,
    Down_16h = gene_name %in% down_16h,
    Up_48h   = gene_name %in% up_48h,
    Down_48h = gene_name %in% down_48h,
    
    sig_16h = (abs(log2fc_16h) > log2fc_threshold & padj_16h < padj_threshold),
    sig_48h = (abs(log2fc_48h) > log2fc_threshold & padj_48h < padj_threshold),
    
    sig_category = case_when(
      sig_16h & sig_48h & sign(log2fc_16h) == sign(log2fc_48h) ~ "same_direction",
      sig_16h & sig_48h & sign(log2fc_16h) != sign(log2fc_48h) ~ "opposite_direction",
      sig_16h & !sig_48h ~ "16h_only",
      !sig_16h & sig_48h ~ "48h_only",
      TRUE ~ "not_significant"
    ),
    
    # Option 2: TRUE only if significant in both timepoints, same direction
    thresholds_met = (sig_category == "same_direction")
  ) %>%
  arrange(gene_name)

# write.csv(supp_table,file=paste0(outputFolder,"correlation_df.csv"))


###########################################
# Compare 16h and 48h data -- Correlation plot
###########################################

correlation_df <- read_csv("~/Rstudio/expression_analysis_16h/results/visualizations/correlationcorrelation_df.csv")%>% select(-c(1))

# Calculate correlation

cor_result <- cor.test(correlation_df$log2fc_16h, correlation_df$log2fc_48h)

significant <- correlation_df %>%
  filter(is_significant) %>%
  mutate(across(starts_with("baseMean"), ~ formatC(., format = "f", digits = 2))) # %>% view()

# -- Define DE status and categories -------------------------------------------
padj_threshold   <- 1e-5
log2fc_threshold <- 1

correlation_df <- correlation_df %>%
  mutate(
    de_16h = padj_16h < padj_threshold & abs(log2fc_16h) >= log2fc_threshold,
    de_48h = padj_48h < padj_threshold & abs(log2fc_48h) >= log2fc_threshold
  )

# -- Three panels --------------------------------------------------------------
panel_levels <- c("DE at 16h", "DE at 48h", "Shared")

# Grey background: full data replicated to each panel
background_df <- bind_rows(lapply(panel_levels, function(p) {
  correlation_df %>% mutate(panel = p)
})) %>%
  mutate(panel = factor(panel, levels = panel_levels))

# Colored foreground: per-panel subsets
foreground_df <- bind_rows(
  correlation_df %>% filter(de_16h)           %>% mutate(panel = "DE at 16h"),
  correlation_df %>% filter(de_48h)           %>% mutate(panel = "DE at 48h"),
  correlation_df %>% filter(de_16h & de_48h)  %>% mutate(panel = "Shared")
) %>%
  mutate(panel = factor(panel, levels = panel_levels))

panel_colors <- c(
  "DE at 16h" = "#1b9e77",
  "DE at 48h" = "#D55E00",
  "Shared"    = "black"
)

# -- Plot  -------------------------------------------------------------------
correlation_plot <- ggplot(mapping = aes(x = log2fc_16h, y = log2fc_48h)) +
  
  geom_abline(intercept = 0, slope = 1, linewidth = 1.0,
              color = "black", linetype = "dashed") +
  geom_hline(yintercept = 0, linewidth = 1.0, color = "black") +
  geom_vline(xintercept = 0, linewidth = 1.0, color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = 4, color = "black", linewidth = 0.75) +
  geom_hline(yintercept = c(-1, 1), linetype = 4, color = "black", linewidth = 0.75) +
  
  # geom_point(data = background_df,
  #            color = "darkgrey", alpha = 0.5, #696969
  #            shape = 1, stroke = 0.6, size = 2.5) +
  
  geom_point(data = foreground_df,
             aes(color = panel),
             alpha = 0.6, shape = 1,
             stroke = 1, size = 2) +
  
  facet_wrap(~ panel, nrow = 1) +
  
  scale_color_manual(values = panel_colors, guide = "none") +
  
  labs(
    title    = "RNA-seq: 48h vs 16h Results",
    subtitle = bquote("r = " * .(round(cor_result$estimate, 3)) *
                        ", p = " * 3.79 %*% 10^-48),
    caption  = bquote("Coding genes; padj ≤ " * 1 %*% 10^-5),
    x        = bquote(bold(~log[2]~ "Fold Change (16h)")),
    y        = bquote(bold(~log[2]~ "Fold Change (48h)"))
  ) +
  
  scale_x_continuous(breaks = seq(-10, 10, by = 2)) +
  scale_y_continuous(breaks = seq(-10, 10, by = 2)) +
  
  theme_bw() +
  theme(
    plot.title        = element_text(size = 20, face = "bold"),
    plot.subtitle     = element_text(size = 16, face = "bold"),
    plot.caption.position = "plot",
    plot.caption      = element_text(size = 16, face = "plain", hjust = 1, colour = "black"),
    axis.text.x       = element_text(family = "Arial", size = 18, colour = "black"),
    axis.text.y       = element_text(family = "Arial", size = 18, colour = "black"),
    axis.title.y      = element_text(family = "Arial", face = "bold", size = 20),
    axis.title.x      = element_text(family = "Arial", face = "bold", size = 20),
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 2),
    axis.line         = element_blank(),
    axis.ticks        = element_line(size = 1.2),
    axis.ticks.length = unit(3, "mm"),
    strip.background  = element_rect(fill = "black", colour = "black"),
    strip.text        = element_text(size = 16, face = "bold",colour = "white")
  )

print(correlation_plot)
