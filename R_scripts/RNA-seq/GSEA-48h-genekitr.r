###############################################################
# GSEA Analysis with genekitr for gene expression data
# Comparing nup211-so (48h) vs. non-shutoff cells (no thiamine)
###############################################################

# Load required libraries
libraries <- c(
  "readr", "dplyr", "tidyr", "tidyverse", "ggplot2","plotly",
  "ggpubr", "gplots", "RColorBrewer", "viridis", "PoiClaClu",
  "GenomicFeatures", "AnnotationDbi", "AnnotationForge",
  "GenomicAlignments", "BiocParallel", "DESeq2", "pheatmap", 
  "org.Spombe.eg.db","clusterProfiler","DOSE", "enrichplot",
  "writexl", "openxlsx", "stringr","genekitr","patchwork", 
  "ggrepel","ggtext","cowplot"
)
for (lib in libraries) {
  if (!require(lib, character.only = TRUE)) {
    BiocManager::install(lib)
    library(lib, character.only = TRUE)
  }
}
# Set output folder
outputFolder <- "~/Rstudio/expression_analysis_48h/star_pipeline/results/GSEA/genekitr/coding_only/"

############################################################################
# 0.1.a Read data in (if necessary)
results_annotated_48h <- read_csv("Desktop/nup211_GEO_sumbission_20221014/processed data files/DEseq2Results.csv") %>%
  column_to_rownames(var = names(.)[1])

filtered_annotated_48h <- results_annotated_48h %>% filter(Type == "protein coding gene") %>%
  rename(gene_type = Type)

# Function to summarize DE results into a dataframe
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
  
summary_table_48h <- bind_rows(
  summarize_DE(results_annotated_48h, "All Transcripts", 0.05),
  summarize_DE(results_annotated_48h, "All Transcripts", 10e-6),
  summarize_DE(filtered_annotated_48h, "Coding Genes Only", 0.05),
  summarize_DE(filtered_annotated_48h, "Coding Genes Only", 10e-6)
)

summary_table_48h <- summary_table_48h %>%
  rename( `p.adjusted cutoff` = padj_Cutoff,
          `2-Fold Up` = `X2.Fold.Up`,
          `2-Fold Down` = `X2.Fold.Down`,
          `4-Fold Up` = `X4.Fold.Up`,
          `4-Fold Down` = `X4.Fold.Down`,
          `2-Fold DE` = `X2.Fold.DE`,
          `4-Fold DE` = `X4.Fold.DE`
  )

# Display as a nicely formatted table in RMarkdown
knitr::kable(summary_table_48h, caption = bquote("Differential Expression Summary: *nup211-so* (48h) vs. control"))


###############################################################
# 1. Parameter Setup - Define filtering thresholds
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
pathways <- c("GO:0071852", "GO:0045893", "GO:0045944", "GO:0002181")
rtqpcr_genes <- c("chp1", "mit1", "mus7", "for3", "tea2", "pho1", "urg2",
                  "atf1","pcr1","mei2","wee1", "knh1","pom1","mbx1", "bgs1","pxl1", "agn2")
volcano_pathways <- c("GO:0002181","GO:0045893", "GO:0071852","GO:0071554","GO:0009250",
                      "GO:0006066", "GO:0046395", "GO:0051253","GO:0043087"
)

# Decoding pathway order:cytoplasmic translation, positive regulation of DNA-templated transcription
# fungal-type cell wall organization or biogenesis, glucan biosynthetic process,
# alcohol metabolic process, carboxylic acid catabolic process, 
# negative regulation of RNA metabolic process, regulation of GTPase activity

###############################################################
#    Data Preparation
###############################################################

#### Coding Only ####
# Prepare gene list for coding genes only
codingList <- filtered_annotated_48h$log2FoldChange
names(codingList) <- filtered_annotated_48h$gene_name
codingList <- sort(codingList, decreasing = TRUE)


###############################################################
#  GSEA Analysis 
###############################################################

# Get gene sets for S. pombe
sp_gs <- geneset::getGO(org = "Schizosaccharomyces pombe", ont = 'bp')


# Run GSEA for coding genes only
gsea_kitr_coding <- genGSEA(genelist = codingList, geneset = sp_gs, min_gset_size = 5, max_gset_size = 250)

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
  filter(setSize >= 10,
         setSize <= max_setSize,
         Count >= min_genes_leading_edge,
         Count <= max_genes_leading_edge,
         abs(enrichmentScore) >= min_enrichment_score,
         abs(NES) >= abs_NES, 
         p.adjust <= max_p_value,
         tags_percentage >= min_tags,
         list_percentage >= min_list,
         signal_percentage >= min_signal,
         RichFactor >=min_RichFactor) %>%
  bind_rows(gsea_kitr_coding_df %>% slice(c(17,20,56))) %>% arrange(p.adjust)


# Examine top enriched pathways
top_enriched <- filtered_gsea_kitr_coding_df %>% 
  dplyr::arrange(desc(NES)) %>% 
  select(Description, NES) %>% 
  slice(1:10)

bottom_enriched <- filtered_gsea_kitr_coding_df %>% 
  dplyr::arrange((NES)) %>% 
  select(Description, NES) %>% 
  head(7)

###############################################################
#  Save 'results' dataframes
###############################################################
# write.csv(gsea_kitr_df,file = "~/Rstudio/expression_analysis_48h/star_pipeline/results/GSEA/genekitr/all_transcripts/results/gsea_results_all_transcripts.csv")
# write.csv(filtered_gsea_kitr,file = "~/Rstudio/expression_analysis_48h/star_pipeline/results/GSEA/genekitr/all_transcripts/results/filtered_gsea_results_all_transcripts.csv")
# write.csv(gsea_kitr_coding_df,file = paste0(outputFolder, "results/gsea_results_coding.csv"))
# write.csv(filtered_gsea_kitr_coding_df,file = paste0(outputFolder, "results/filtered_gsea_results_coding.csv"))


###############################################################
# Visualization
###############################################################


####  Top 20 Bar Plot ####
# Select top and bottom pathways based on NES
filtered_gsea_top20 <- filtered_gsea_kitr_coding
filtered_gsea_top20$gsea_df <- filtered_gsea_kitr_coding_df %>%
  slice(-13)

# Create barplot
gsea_kitr_coding_barplot_t20 <- plotGSEA(filtered_gsea_top20, 
                                         plot_type = "bar", 
                                         show_pathway = 5, 
                                         main_text_size = 18,
                                         label_by = "description", 
                                         stats_metric = "p.adjust",
                                         colour = c("blue", "red")) +
  ggtitle(bquote(bold("GSEA results:" ~ bolditalic("nup211-so") ~ "(48h) vs. control" ~ "\n"))) +
  labs(y = "Normalized Enrichment Score") + 
  theme(axis.text.y = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 18, face = "plain"),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.title.x = element_text(size = 18, face = "bold"),
        plot.title = element_text(size = 18, face = "bold"))

print(gsea_kitr_coding_barplot_t20)
# Save plot
# ggsave(paste0(outputFolder, "plots/GSEA_top20_barplot.pdf"), 
#       gsea_kitr_coding_barplot_t20, width = 13, height = 8, units = "in", dpi = 600)


###
# corresponding barplot
filtered_gsea_publication <- filtered_gsea_kitr_coding
filtered_gsea_publication$gsea_df <- filtered_gsea_kitr_coding_df
gsea_kitr_coding_barplot_t20 <- plotGSEA(filtered_gsea_publication, 
                                         plot_type = "bar", 
                                         #show_pathway = 5, 
                                         main_text_size = 18,
                                         label_by = "description", 
                                         stats_metric = "p.adjust",
                                         colour = c("blue", "red")) +
  ggtitle(bquote(bold("GSEA results:" ~ bolditalic("nup211-so") ~ "(48h) vs. control" ~ "\n"))) +
  labs(y = "Normalized Enrichment Score") + 
  theme(axis.text.y = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 18, face = "plain"),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.title.x = element_text(size = 18, face = "bold"),
        plot.title = element_text(size = 18, face = "bold"))

print(gsea_kitr_coding_barplot_t20)


###############################################################
# Custom Volcano Plot
###############################################################

# 48h GSEA results (Supplemental Figure) ----------------------------------

pathways_48h <- c("GO:0044282", "GO:0006066", "GO:0005996", "GO:0002181", "GO:0006468", "GO:0035556", 
                  "GO:0045893", "GO:0071852", "GO:0009250", "GO:0050801", "GO:0043087","GO:0006979",
                  "GO:2001141","GO:1901575")

filtered_coding_gsea_48h <- read_csv("~/Rstudio/expression_analysis_48h/star_pipeline/results/GSEA/genekitr/coding_only/results/filtered_gsea_results_coding.csv")

gsea_df_48h <- filtered_coding_gsea_48h %>% 
  mutate( RichFactor = Count / setSize, 
          NES_group = case_when(
            NES > 1. ~ "Up",
            NES < -1.5  ~ "Down",
            TRUE ~ "NS"
          ),
          label = ifelse(NA_BP_ID %in% pathways_48h, Description, NA)  # or use a different column
  ) 


GSEA_volcano_48h <- ggplot(subset(gsea_df_48h, !is.na(label)), aes(x = NES, y = neglog10padj, fill = RichFactor, size = Count))  +
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
  )+
  # Added order = 2 to push Count below RichFactor
  scale_size(
    range = c(6, 14), 
    breaks = scales::pretty_breaks(n = 3),
    guide = guide_legend(order = 2) 
  ) + 
  # Added order = 1 to pull RichFactor to the top
  scale_fill_viridis_c(
    option = "magma", 
    direction = 1, 
    trans = "sqrt", 
    limits = c(0.4, 0.87),
    guide = guide_colorbar(order = 1) 
  ) +
  scale_y_continuous(
    limits = c(1, 3.6),
    breaks = seq(1.5, 4, 0.5),
    expand = c(0, 0)
  ) +
  scale_x_continuous(limits = c(-3.5, 2.5),
                     breaks = seq(-3, 3, 1))+
  labs(
    title = "GSEA results",
    subtitle = "<i>nup211-so </i>(48h) vs. <i>nup211-so </i>(-T) \n",
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
    axis.line = element_line(size = 1.0), # Control axis thickness
    axis.ticks = element_line(size = 1.2), # Control tick mark thickness
    axis.ticks.length = unit(3, "mm"),
    legend.title = element_text(size = 14, face = "bold"),
    plot.subtitle = ggtext::element_markdown(size = 16, face = "plain"),
    legend.position = "right",
    legend.text = element_text(size = 14)
  )  #+ coord_fixed(ratio = 1)


print(GSEA_volcano_48h)   #(save as 8 in x 11 in landscape) or 7.95 by 9.85 or 885 px by 739 px ratio


###############################################################
#
#
#    GSEA analysis - Analyzing "leading edge" genes
#
#
###############################################################
outputFolder = "~/Rstudio/expression_analysis_48h/star_pipeline/results/GSEA/genekitr/"


leading_edge_all_transcripts <- gsea_kitr_coding_df %>% select(c(1,2,5:8,11)) %>%
                                rename(genes = geneID,
                                       go_id=fission_yeast_BP_ID,
                                       term=Description) %>%
                                mutate(genes = gsub("/", ",", genes))

leading_edge_all_transcripts <- leading_edge_all_transcripts %>%
  mutate(gene_name = strsplit(genes, ",")) %>%
  unnest(gene_name) %>% select(-genes) %>% left_join(results_annotated_48h, by = "gene_name") %>% 
  mutate(FoldChange = 2^(log2FoldChange)) %>%  # Calculate linear fold change
  select(c(1,2,3,5:7,24,10,9,15,17)) %>% 
  rename(padj.GSEA=p.adjust,
         padj.FC=padj)

# write_csv(leading_edge_all_transcripts, paste0(outputFolder,"all_transcripts/results/leading_edge_all_transcripts.csv"))

###########################################
# Organize information by individual genes  
###########################################

leading_edge_all_transcripts_by_gene <- leading_edge_all_transcripts %>%  # Join with mapping
  group_by(gene_name) %>%
  summarise(term = paste(term, collapse = ", "),
            go_id = paste(go_id, collapse = ", "), # Recombine names into a single string
            .groups = 'drop') %>% left_join(results_annotated_48h, by = "gene_name") %>%
  mutate(terms = str_count(term, ",") + 1,
         FoldChange = 2^(log2FoldChange)) %>%
  rename(term_names = term,
         padj.FC=padj) %>% 
  select(c(1,13,21,6,11,5,20,3,2,))  # 15:18

# write_csv(leading_edge_all_transcripts, paste0(outputFolder,"all_transcripts/results/leading_edge_all_transcripts_by_gene.csv"))

########################################################
# Filter out non-significantly affected genes (optional)  
########################################################

leading_edge_all_transcripts_by_gene_filtered <-filter(leading_edge_all_transcripts_by_gene, abs(log2FoldChange) > 1  & padj.FC < 10e-6)
# write_csv(leading_edge_all_transcripts_by_gene_filtered, paste0(outputFolder,"all_transcripts/results/significant_leading_edge_all_transcripts.csv"))

########################################################
# Analyze genes associated with specific pathways  
########################################################

leading_edge_all_transcripts_by_gene %>% filter(str_detect(term_names, "cell wall")) %>%
  select(c(gene_name,FoldChange,padj.FC,Description)) 

leading_edge_all_transcripts_by_gene_filtered %>% filter(str_detect(term_names, "cell wall")) %>%
  select(c(gene_name,FoldChange,padj.FC,Description))


########################################################
# Cherry-pick genes of interest for examination  
########################################################

manual_list <- c("swr1", "snf21", "fzr2", "nup189", "mug5", "mok11", "slr1", "snu66", "Tf2-2", "shd1", "chp1", "cut4", "set3", 
                 "sck2", "rtr1", "prp10", "ino80", "snt2", "myp2", "spo3", "mst1", "mus7", "nog2", "pcp1", "slp1", "brr2", 
                 "cut3", "klp6", "tf2-10", "dbp7", "rrp2", "rad13", "prt1", "nop2", "spt2", "spc7", "prp11", "snf22", "set1", 
                 "pcm1", "rhp42", "prp28", "for3", "mdn1", "sfc3","bub1","mad1","sad1", "sfi1","spo15","dma1","SPBC1539.02")


leading_edge_all_transcripts_by_gene_filtered %>% filter(gene_name %in% manual_list) %>%
  select(c(gene_name,FoldChange,padj.FC, baseMean,Description)) # %>% view()

genes_of_interest <- results_annotated_48h %>% filter(gene_name %in% manual_list) %>%
  mutate(FoldChange = 2^(log2FoldChange)) %>%
  select(c(gene_name,FoldChange,padj,pvalue,Description)) # %>% view()

common_genes <- intersect(
  manual_list,
  leading_edge_all_transcripts_by_gene_filtered$gene_name
)

genes_of_interest2 <- results_annotated_48h %>% filter(gene_name %in% common_genes) %>%
  mutate(FoldChange = 2^(log2FoldChange)) %>%
  select(c(gene_name,FoldChange,padj,baseMean,Description)) # %>% view()


# write.csv(genes_of_interest2, paste0(outputFolder,"all_transcripts/results/potential_genes_of_interest.csv"))


###############################################################
#
#    Section 2 - "coding_only" GSEA Output
###############################################################
leading_edge_coding_only <- filtered_gsea_kitr_coding_df %>% select(c(1,2,5:8,11)) %>%      # or gsea_kitr_coding_df %>% select(c(1,2,5:8,11)) 
  rename(genes = geneID,
         go_id=fission_yeast_BP_ID,
         term=Description) %>%
  mutate(genes = gsub("/", ",", genes))

leading_edge_coding_only <- leading_edge_coding_only %>%
  mutate(gene_name = strsplit(genes, ",")) %>%
  unnest(gene_name) %>% select(-genes) %>% left_join(results_annotated_48h, by = "gene_name") %>% 
  mutate(FoldChange = 2^(log2FoldChange)) %>%  # Calculate linear fold change
  select(c(1,2,3,5:7,24,10,9,15,17)) %>% 
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
            .groups = 'drop') %>% left_join(results_annotated_48h, by = "gene_name") %>%
  mutate(terms = str_count(term, ",") + 1,
         FoldChange = 2^(log2FoldChange)) %>%
  rename(term_names = term,
         padj.FC=padj) %>% 
  select(c(1,13,21,6,11,5,20,3,2))  # 15:18

# write_csv(leading_edge_coding_only_by_gene, paste0(outputFolder,"results/leading_edge_coding_only_by_gene.csv"))

########################################################
# Filter out non-significantly affected genes (optional)  
########################################################

leading_edge_coding_only_by_gene_filtered <-filter(leading_edge_coding_only_by_gene, abs(log2FoldChange) > 1  & padj.FC < 10e-6)
# write_csv(leading_edge_coding_only_by_gene_filtered, paste0(outputFolder,"results/significant_leading_edge_coding_only.csv"))
########################################################
# Analyze genes associated with specific pathways  
########################################################

leading_edge_coding_only_by_gene %>% filter(str_detect(term_names, "cell wall")) %>%
  select(c(gene_name,FoldChange,padj.FC,Description)) 

leading_edge_coding_only_by_gene_filtered %>% filter(str_detect(term_names, "cell wall")) %>%
  select(c(gene_name,FoldChange,padj.FC,Description))


########################################################
# Cherry-pick genes of interest for examination  
########################################################

manual_list <- c("swr1", "snf21", "fzr2", "nup189", "mug5", "mok11", "slr1", "snu66", "Tf2-2", "shd1", "chp1", "cut4", "set3", 
                 "sck2", "rtr1", "prp10", "ino80", "snt2", "myp2", "spo3", "mst1", "mus7", "nog2", "pcp1", "slp1", "brr2", 
                 "cut3", "klp6", "tf2-10", "dbp7", "rrp2", "rad13", "prt1", "nop2", "spt2", "spc7", "prp11", "snf22", "set1", 
                 "pcm1", "rhp42", "prp28", "for3", "mdn1", "sfc3","bub1","mad1","sad1", "sfi1","spo15","dma1","SPBC1539.02")


leading_edge_coding_only_by_gene_filtered %>% filter(gene_name %in% manual_list) %>%
  select(c(gene_name,FoldChange,padj.FC, baseMean,Description)) # %>% view()

genes_of_interest <- results_annotated_48h %>% filter(gene_name %in% manual_list) %>%
  mutate(FoldChange = 2^(log2FoldChange)) %>%
  select(c(gene_name,FoldChange,padj,baseMean,Description)) # %>% view()

common_genes <- intersect(
  manual_list,
  leading_edge_coding_only_by_gene_filtered$gene_name
)

genes_of_interest2 <- results_annotated_48h %>% filter(gene_name %in% common_genes) %>%
  mutate(FoldChange = 2^(log2FoldChange)) %>%
  select(c(gene_name,FoldChange,padj,baseMean,Description)) # %>% view()


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


