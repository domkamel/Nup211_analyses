###############################################################
# GSEA Analysis with genekitr for gene expression data
# Comparing nup211-so (16h) vs. control (558 + thiamine)
###############################################################

# Load required libraries
libraries <- c(
  "readr", "dplyr", "tidyr", "tidyverse", "ggplot2","plotly",
  "ggpubr", "gplots", "RColorBrewer", "viridis", "PoiClaClu",
  "GenomicFeatures", "AnnotationDbi", "AnnotationForge",
  "GenomicAlignments", "BiocParallel", "DESeq2", "pheatmap", 
  "org.Spombe.eg.db","clusterProfiler","DOSE", "enrichplot",
  "writexl", "openxlsx", "stringr","genekitr","patchwork", 
  "ggrepel", "ggtext","cowplot"
)
for (lib in libraries) {
  if (!require(lib, character.only = TRUE)) {
    BiocManager::install(lib)
    library(lib, character.only = TRUE)
  }
}
# Set output folder
outputFolder <- "~/Rstudio/expression_analysis_16h/results/GSEA/genekitr/coding_only/"

############################################################################
# Read data in (if necessary)
results_annotated <- read_csv("Rstudio/expression_analysis_16h/results/DE_spreadsheets/DESeq2_results_all_transcripts.csv") %>%
  column_to_rownames(var = names(.)[1])

filtered_annotated <- read_csv("Rstudio/expression_analysis_16h/results/DE_spreadsheets/DESeq2_results_protein_coding.csv") %>%
  column_to_rownames(var = names(.)[1])

###############################################################
#  Parameter Setup - Define filtering thresholds
###############################################################
min_setSize <- 20                # Minimum gene set size
max_setSize <- 200               # Maximum gene set size
min_genes_leading_edge <- 10     # Minimum genes in leading edge
max_genes_leading_edge <- 110    # Maximum genes in leading edge
min_enrichment_score <- 0.5      # Minimum enrichment score 
abs_NES <- 1.5                   # Minimum absolute normalized enrichment score
max_p_value <- 0.05              # Maximum adjusted p-value
min_tags <- 0.4                  # Minimum tags percentage
min_list <- 0.1                  # Minimum list percentage 
min_signal <- 0.3                # Minimum signal percentage
min_RichFactor <- 0.4             # Minimum RichFactor

# Define pathways and genes of interest
pathways <- c("GO:0000278", "GO:0035556", "GO:0000075", "GO:0051276", "GO:0006338")
rtqpcr_genes <- c("chp1", "mit1", "swi6", "for3", "tea2", "mus7", "not1")
volcano_pathways <- c("GO:0000278", "GO:0000075", "GO:0007059", "GO:0006338", "GO:0035556",
                      "GO:0006397", "GO:0008380", "GO:0006281", "GO:0006913",
                      "GO:0009060", "GO:0002181", "GO:0009150", "GO:0006259", "GO:0006091","GO:0045333")

###############################################################
#    Data Preparation
###############################################################

# Prepare gene list for coding genes only
codingList <- filtered_annotated$log2FoldChange
names(codingList) <- filtered_annotated$gene_name
codingList <- sort(codingList, decreasing = TRUE)

###############################################################
#    GSEA Analysis 
###############################################################

# Get gene sets for S. pombe
sp_gs <- geneset::getGO(org = "Schizosaccharomyces pombe", ont = 'bp')

# Run GSEA for coding genes only

library(genekitr)
sp_gs <- geneset::getGO(org = "Schizosaccharomyces pombe", ont = 'bp')
gsea_kitr_coding <- genGSEA(genelist = codingList, geneset = sp_gs, min_gset_size = 5)

# Extract results dataframe
gsea_kitr_coding_df <- gsea_kitr_coding$gsea_df

# Process leading edge information and calculate RichFactor
gsea_kitr_coding_df <- gsea_kitr_coding_df %>%
  mutate(tags_percentage   = as.numeric(str_extract(leading_edge, "(?<=tags=)(\\d+\\.?\\d*)")) / 100,
         list_percentage   = as.numeric(str_extract(leading_edge, "(?<=list=)(\\d+\\.?\\d*)")) / 100,
         signal_percentage = as.numeric(str_extract(leading_edge, "(?<=signal=)(\\d+\\.?\\d*)")) / 100,
         RichFactor = Count/setSize,
         neglog10padj = -1 * log10(p.adjust))

# Filter results based on thresholds
filtered_gsea_kitr_coding_df <- gsea_kitr_coding_df %>%
  filter(setSize >= min_setSize,
         setSize <= max_setSize,
         Count >= min_genes_leading_edge,
         Count <= max_genes_leading_edge,
         abs(enrichmentScore) >= min_enrichment_score,
         abs(NES) >= abs_NES, 
         p.adjust <= max_p_value,
         tags_percentage >= min_tags,
         list_percentage >= min_list,
         signal_percentage >= min_signal,
         RichFactor >=min_RichFactor)

# Examine top enriched pathways
top_enriched <- filtered_gsea_kitr_coding_df %>% 
  dplyr::arrange(desc(NES)) %>% 
  dplyr::mutate(neg_log_p = -10 * log10(p.adjust)) %>% # Added the transformation
  dplyr::select(Description, NES, p.adjust, neg_log_p) %>% 
  dplyr::slice(1:10)

bottom_enriched <- filtered_gsea_kitr_coding_df %>% 
  dplyr::arrange(NES) %>% 
  dplyr::select(Description, NES) %>%  # Added dplyr:: here
  head(7)

###############################################################
#   Save 'results' dataframes
###############################################################

# write.csv(gsea_kitr_coding_df,file = paste0(outputFolder, "results/gsea_results_coding.csv"))
# write.csv(filtered_gsea_kitr_coding_df,file = paste0(outputFolder, "results/filtered_gsea_results_coding.csv"))

###############################################################
# 4. Visualization
###############################################################

# Create classic GSEA plot for specific pathways

# Prepare data for plotting
filtered_gsea_kitr_coding <- gsea_kitr_coding                          # First, copy the genekitr output to a "new" object
filtered_gsea_kitr_coding$gsea_df <- filtered_gsea_kitr_coding_df      # Then modify overwrite the gsea_df object
                                                                       # (This is to prevent overwriting the original results!)


Pgseakitr_coding_division <- plotGSEA(gsea_kitr_coding, 
                                     plot_type = "classic", 
                                     show_pathway = "GO:0140014",
                                     label_by = c("description"),
                                     colour = "red",
                                     stats_metric = "p.adjust", 
                                     legend_text = 12,
                                     font_type = "sans",
                                     main_text_size = 18)

Pgseakitr_division <- Pgseakitr_coding_division + 
  # Using a cleaner title setup
  ggtitle(bquote(bold("GSEA Results: ") ~ bolditalic("nup211-so") ~ bold(" (16h) vs. wild-type"))) +
  #theme_bw(base_family = "sans") + # Sets a clean white background
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    #panel.grid.minor = element_blank(), # Removes distracting small grid lines
    legend.position = "right",
    plot.margin = margin(10, 10, 10, 10) # Adds padding around the whole graph
  )

print(Pgseakitr_division)

###################################################
# Custom barplot using genekitr built in plotting
###################################################

custom_list <- c("GO:0006325","GO:0008380","GO:0006468", "GO:0009150","GO:0000075",
                "GO:0035556", "GO:0002181","GO:0000278", "GO:0007059","GO:0016071",
                "GO:0006259", "GO:0006281","GO:0051276", "GO:0006091","GO:0009060","GO:0006913") #maybe remove "GO:0006913"
filtered_gsea_publication <- filtered_gsea_kitr_coding
filtered_gsea_publication_df <- filtered_gsea_publication$gsea_df %>% filter(fission_yeast_BP_ID %in% custom_list)
filtered_gsea_publication$gsea_df <- filtered_gsea_publication_df

gsea_kitr_coding_barplot_pub <- plotGSEA(filtered_gsea_publication, 
                                         plot_type = "bar", 
                                         #show_pathway = 5, 
                                         main_text_size = 18,
                                         label_by = "description", 
                                         stats_metric = "p.adjust",
                                         colour = c("blue", "red")) +
  ggtitle(bquote(bold("GSEA results:" ~ bolditalic("nup211-so") ~ "(16h) vs. control" ~ "\n"))) +
  labs(y = "Normalized Enrichment Score") + 
  theme(axis.text.y = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 18, face = "plain"),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.title.x = element_text(size = 18, face = "bold"),
        plot.title = element_text(size = 18, face = "bold"))

print(gsea_kitr_coding_barplot_pub)

# Save plot
# ggsave(paste0(outputFolder, "plots/GSEA_top20_barplot.pdf"), 
#       gsea_kitr_coding_barplot_t20, width = 13, height = 8, units = "in", dpi = 600)

################################################################
####  Create custom volcano plots! ####
################################################################

filtered_coding_gsea_16h <- read_csv("Rstudio/expression_analysis_16h/results/GSEA/genekitr/coding_only/results/filtered_gsea_results_coding.csv") %>%
  column_to_rownames(var = names(.)[1])


# Removed "GO:0006468" (Protein phosphorylation), "GO:0035556" (intracellular signal transduction)

pathways_16h <- c("GO:0006325","GO:0008380", "GO:0009150","GO:0000075",  #"GO:0006468","GO:0035556",
                  "GO:0002181","GO:0000278", "GO:0007059","GO:0016071",
                  "GO:0006259", "GO:0006281","GO:0051276", "GO:0006091",
                  "GO:0009060","GO:0006913")                             #maybe remove "GO:0006913"

gsea_df_16h <- filtered_coding_gsea_16h %>% 
  mutate( NES_group = case_when(
    NES > 1. ~ "Up",
    NES < -1.5  ~ "Down",
    TRUE ~ "NS"
  ),
  label = ifelse(rownames(filtered_coding_gsea_16h) %in% pathways_16h, Description, NA)  # or use a different column
  ) 

GSEA_volcano_16h <- ggplot(subset(gsea_df_16h, !is.na(label)), aes(x = NES, y = neglog10padj, fill = RichFactor, size = Count))  +
  geom_point(shape = 21, alpha = 0.8, color = "black", stroke = 1) +
  geom_vline(xintercept = c(-1.5, 1.5), linetype = 4, color = "black", linewidth = 0.75) +
  geom_hline(yintercept = -log10(0.05), linetype = 4, color = "black", linewidth = 0.75) +
  geom_label_repel(
    aes(label = label),
    fill = "white", color = "black", fontface="bold", size = 4, box.padding = 0.6,
    max.overlaps = 1000000, segment.color = "black",
    segment.size = 0.5,         # Thickness of the connecting line
    segment.curvature = 0.2,    # Add curvature to connectors if desired
    nudge_x = -0.5,
    nudge_y = 0.2,
    seed = 32
  ) +
  scale_fill_viridis_c(option = "plasma", guide = "colourbar", direction = 1,
                       trans = "sqrt",
                       limits = c(0.4, 0.87)) +
  scale_size(range =c(6,14), breaks = scales::pretty_breaks(n = 3)) +
  scale_y_continuous(
    limits = c(1, 6.75),
    breaks = seq(2, 6, 1),
    expand = c(0, 0)
  ) +
  scale_x_continuous(limits = c(-3.5, 2.5),
                     breaks = seq(-3, 3, 1))+
  labs(
    title = "GSEA results",
    subtitle = "<i>nup211-so </i>(16h) vs. wild-type (+T) \n",
    x = "Normalized Enrichment Score (NES)",
    y =  bquote(bold(-log[10]~"(P)")),
    fill = "RichFactor"
  ) +
  theme_pubr() +
  theme(
    plot.title = element_text(size = 18, hjust = 0, face = "bold"),
    panel.grid.minor = element_line(linewidth =  0.6, linetype = 3, color = "lightgray"),
    panel.grid.major = element_line(linewidth =  0.6, linetype = 3, color = "darkgrey"),
    axis.title.y = element_text(family= 'Arial', face = 'bold', size = 18,
                                margin = margin(r = 10)),
    axis.title.x = element_text(family= 'Arial', face = 'bold', size = 18,
                                margin = margin(t = 10)),
    axis.text.x = element_text(family= 'Arial', size = 18, face= "plain",colour = 'black'),
    axis.text.y = element_text(family= 'Arial', size = 18, face= "plain",colour = 'black'),
    axis.line = element_line(linewidth =  1.0), # Control axis thickness
    axis.ticks = element_line(linewidth =  1.2), # Control tick mark thickness
    axis.ticks.length = unit(3, "mm"),
    legend.title = element_text(size = 14, face = "bold"),
    plot.subtitle = ggtext::element_markdown(size = 16, face = "plain"),
    legend.position = "right",
    legend.text = element_text(size = 14)
  ) #+ coord_fixed(ratio = 1)


print(GSEA_volcano_16h)   #(save as 9.32 in x 7.83 in landscape) or 7.95 by 9.85 or 885 px by 739 px ratio


###################################################
# Custom barplot using ggplot2, stringr, etc. 
###################################################

keep_ids <- c(
  # Down
  "GO:0006091",  # generation of precursor metabolites and energy
  "GO:0009060",  # aerobic respiration
  "GO:0002181",  # cytoplasmic translation
  # Up
  "GO:0000278",  # mitotic cell cycle
  "GO:0006281",  # DNA repair
  "GO:0006325",  # chromatin organization
  "GO:0016071",  # mRNA metabolic process
  "GO:0007059",  # chromosome segregation
  "GO:0006913",  # nucleocytoplasmic transport
  "GO:0035556"   # intracellular signal transduction  (optional)
)

gsea_bar_df <- gsea_df_16h[keep_ids, ]

library(stringr)

# sentence case, don't mess up mRNA or DNA

gsea_df_16h <- gsea_df_16h[keep_ids, ] |>
  mutate(
    label_cap = case_when(
      label == "mRNA metabolic process"        ~ "mRNA metabolic process",
      label == "DNA repair"                    ~ "DNA repair",
      TRUE ~ str_replace(label, "^([a-z])", toupper)
    ),
    is_headline = label %in% c("chromatin organization",
                               "mRNA metabolic process",
                               "mitotic cell cycle"),
    label_styled = ifelse(is_headline,
                          paste0("**", label_cap, "**"),
                          label_cap)
  )

library(ggtext)
gsea_bar <-  ggplot(gsea_df_16h, aes(x = NES, y = reorder(label_styled, NES),  fill = NES_group)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.6) +
  geom_vline(xintercept = 0, linewidth = 1, color = "black") +
  scale_fill_manual(values = c(Up = "red2", Down = "royalblue")) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  labs(x = "Normalized Enrichment Score (NES)", y = NULL) +
  #theme_classic(base_family = "Arial", base_size = 14) +
  theme(
    legend.position   = "none",
    legend.title      = element_blank(),
    legend.text       = element_text(family = 'Arial', face = 'bold', size = 16),
    axis.title.x      = element_text(family = 'Arial', face = 'bold', size = 20),
    axis.title.y      = element_text(family = 'Arial', face = 'bold', size = 20),
    axis.text.x       = element_text(family = 'Arial', size = 18, colour = 'black'),
    axis.text.y       = element_markdown(family = 'Arial', size = 16, colour = 'black'),
    axis.ticks.y      = element_line(linewidth = 1.2), #element_blank(),
    axis.ticks.x      = element_line(linewidth = 1.2),
    axis.ticks.length = unit(3, "mm"),
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 2.0),
    axis.line         = element_blank()
  )

# ggsave(plot = gsea_bar,
#        filename = "SciCom_gsea_bar.pdf",
#        width = 11, height = 5, units = "in", dpi = 600, device = cairo_pdf)


###############################################################
#
#    GSEA analysis - Analyzing "leading edge" genes
#
###############################################################
outputFolder = "~/Rstudio/expression_analysis_16h/results/GSEA/genekitr/"

leading_edge_coding_only <- gsea_kitr_coding_df %>% select(c(1,2,5:8,11)) %>%
  rename(genes = geneID,
         go_id=fission_yeast_BP_ID,
         term=Description) %>%
  mutate(genes = gsub("/", ",", genes))

leading_edge_coding_only <- leading_edge_coding_only %>%
  mutate(gene_name = strsplit(genes, ",")) %>%
  unnest(gene_name) %>% select(-genes) %>% left_join(results_annotated, by = "gene_name") %>% 
  mutate(FoldChange = 2^(log2FoldChange)) %>%  # Calculate linear fold change
  select(c(1:3,5:7,23,10,9,14,15)) %>% 
  rename(padj.GSEA=p.adjust,
         padj.FC=padj)

# write_csv(leading_edge_coding_only, paste0(outputFolder,"coding_only/results/leading_edge_coding_only.csv"))

###########################################
# Organize information by individual genes  
###########################################

leading_edge_coding_only_by_gene <- leading_edge_coding_only %>%  # Join with mapping
  group_by(gene_name) %>%
  summarise(term = paste(term, collapse = ", "),
            go_id = paste(go_id, collapse = ", "), # Recombine names into a single string
            .groups = 'drop') %>% left_join(results_annotated, by = "gene_name") %>%
  mutate(terms = str_count(term, ",") + 1,
         FoldChange = 2^(log2FoldChange)) %>%
  rename(term_names = term,
         padj.FC=padj) %>% 
  select(c(1,20,6,5,10,11,19,3,2))  # 13:16

# write_csv(leading_edge_coding_only, paste0(outputFolder,"coding_only/results/leading_edge_coding_only_by_gene.csv"))

########################################################
# Filter out non-significantly affected genes (optional)  
########################################################

leading_edge_coding_only_by_gene_filtered <-filter(leading_edge_coding_only_by_gene, abs(log2FoldChange) > 1  & padj.FC < 10e-6)
# write_csv(leading_edge_coding_only_by_gene_filtered, paste0(outputFolder,"coding_only/results/significant_leading_edge_coding_only.csv"))

########################################################
# Analyze genes associated with specific pathways         (example)
########################################################

leading_edge_coding_only_by_gene_filtered %>% filter(str_detect(term_names, "cell wall")) %>%
  select(c(gene_name,FoldChange,padj.FC,gene_product))


########################################################
# Cherry-pick genes of interest for examination  
########################################################

manual_list <- c("swr1", "snf21", "fzr2", "nup189", "mug5", "mok11", "slr1", "snu66", "Tf2-2", "shd1", "chp1", "cut4", "set3", 
                 "sck2", "rtr1", "prp10", "ino80", "snt2", "myp2", "spo3", "mst1", "mus7", "nog2", "pcp1", "slp1", "brr2", 
                 "cut3", "klp6", "tf2-10", "dbp7", "rrp2", "rad13", "prt1", "nop2", "spt2", "spc7", "prp11", "snf22", "set1", 
                 "pcm1", "rhp42", "prp28", "for3", "mdn1", "sfc3","bub1","mad1","sad1", "sfi1","spo15","dma1","SPBC1539.02")


leading_edge_coding_only_by_gene_filtered %>% filter(gene_name %in% manual_list) %>%
  select(c(gene_name,FoldChange,padj.FC, baseMean,gene_product)) # %>% view()

genes_of_interest <- results_annotated %>% filter(gene_name %in% manual_list) %>%
  mutate(FoldChange = 2^(log2FoldChange)) %>%
  select(c(gene_name,FoldChange,padj,baseMean,gene_product)) # %>% view()

common_genes <- intersect(
  manual_list,
  leading_edge_coding_only_by_gene_filtered$gene_name
)

genes_of_interest2 <- results_annotated %>% filter(gene_name %in% common_genes) %>%
  mutate(FoldChange = 2^(log2FoldChange)) %>%
  select(c(gene_name,FoldChange,padj,baseMean,gene_product)) # %>% view()


# write.csv(genes_of_interest2, paste0(outputFolder,"coding_only/results/potential_genes_of_interest.csv"))


###########################################
# Save everything to an excel file 
###########################################

write_gsea_sheets <- function(original_df, expanded_df, by_gene_df, output_filename = paste0(outputFolder,"results/GSEA_Subsets.xlsx")) {
  # Create a new workbook
  wb <- createWorkbook()
  
  # Get unique GO Terms
  go_terms <- expanded_df$go_id[!duplicated(expanded_df$go_id)]
  
  # Add original unexpanded dataframe as first sheet
  addWorksheet(wb, "gsea_by_term")
  writeData(wb, "gsea_by_term", original_df)
  
  # Create a subset for each GO Term in the expanded dataframe
  gsea_term_subsets <- split(expanded_df, expanded_df$go_id)
  
  addWorksheet(wb, "gsea_unnested")
  writeData(wb, "gsea_unnested", expanded_df)
  
  addWorksheet(wb, "gsea_by_gene")
  writeData(wb, "gsea_by_gene", by_gene_df)
  
  # Write each GO Term subset to a separate sheet
  for (term in names(gsea_term_subsets)) {
    # Create a safe sheet name (remove invalid characters)
    safe_term <- gsub(":", "_", term)
    
    # Extract metadata from original_df
    term_info <- unique(original_df[original_df$go_id == term, c("term", "NES", "p.adjust")])
    
    # Convert to vector
    term_header <- c(
      paste0("Term: ", term_info$term),
      paste0("NES: ", term_info$NES),
      paste0("p.adj: ", term_info$p.adjust)
    )
    # Ensure the vector has exactly 3 elements (fill with NA if needed)
    if (length(term_header) < 3) {
      term_header <- c(term_header, rep("NA", 3 - length(term_header)))
    }
    
    # Add worksheet
    addWorksheet(wb, safe_term)
    
    # Write metadata in the first row
    writeData(wb, safe_term, t(data.frame(term_header)), startRow = 1, colNames = FALSE)
    
    # Write the rest of the data starting from the second row
    writeData(wb, safe_term, gsea_term_subsets[[term]], startRow = 2)
  }
  
  # Reorder sheets to match the original GO Term order
  sheet_names <- c("gsea_by_term","gsea_unnested","gsea_by_gene", gsub(":", "_", go_terms))
  worksheetOrder(wb) <- match(sheet_names, wb$sheet_names)
  
  # Save the workbook
  saveWorkbook(wb, output_filename, overwrite = TRUE)
  
  return(wb)
}

# write_gsea_sheets(leading_edge_genes_df,leading_edge_genes_expanded,leading_edge_by_gene)

