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

############################################################################
# Read data in (if necessary)
results_annotated <- read_csv("Rstudio/expression_analysis_16h/results/DE_spreadsheets/DESeq2_results_all_transcripts.csv") %>%
  column_to_rownames(var = names(.)[1])

filtered_annotated <- read_csv("Rstudio/expression_analysis_16h/results/DE_spreadsheets/DESeq2_results_protein_coding.csv") %>%
  column_to_rownames(var = names(.)[1])

############################################################################
# Set output folder
outputFolder <- "~/Rstudio/expression_analysis_16h/results/GO/genekitr/coding_only/"

###############################################################
# Data Preparation
###############################################################

#### Coding Only ####
# Prepare gene list for coding genes only
DE_coding <- filtered_annotated %>%                                               #think about changing to results_annotated object
             filter(filtered_annotated$padj < 10e-6) %>%                          #think about changing to 0.05
             mutate(regulation = case_when( log2FoldChange > 2 ~ "up4Fold",
                                           log2FoldChange > 1 ~ "up2Fold",
                                           log2FoldChange < -2 ~ "down4Fold",
                                           log2FoldChange < -1 ~ "down2Fold",
                                           TRUE ~ NA_character_)) %>%
            filter(!is.na(regulation))

codingList <- DE_coding$log2FoldChange
names(codingList) <- DE_coding$gene_name
codingList <- sort(codingList, decreasing = TRUE)

LFC1_coding <- names(codingList)[abs(codingList) > 1]
LFC2_coding <- names(codingList)[abs(codingList) > 2]

up2_coding <- names(codingList)[codingList > 1] #head(names(codingList), 100)
up4_coding <- names(codingList)[codingList > 2] #head(names(codingList), 100)
down2_coding <- names(codingList)[codingList < -1]#names(codingList)[codingList < -1]
down4_coding <- names(codingList)[codingList < -2]#names(codingList)[codingList < -1]


###############################################################
# GO Analysis 
###############################################################
library(genekitr)
sp_gs <- geneset::getGO(org = "Schizosaccharomyces pombe", ont = 'bp')
str(sp_gs)

up2_go_coding <- genORA(up2_coding, geneset = sp_gs,min_gset_size = 5, max_gset_size = 200)
up4_go_coding <- genORA(up4_coding, geneset = sp_gs,min_gset_size = 5, max_gset_size = 200)
down2_go_coding <- genORA(down2_coding, geneset = sp_gs, min_gset_size = 5 )
down4_go_coding <- genORA(down4_coding, geneset = sp_gs, min_gset_size = 5 )


LFC1_go_coding <- genORA(LFC1_coding, geneset = sp_gs,min_gset_size = 5, max_gset_size = 200)
LFC2_go_coding <- genORA(LFC2_coding, geneset = sp_gs,min_gset_size = 5, max_gset_size = 200)

up_go_test <- bind_rows(up4_go_coding,up2_go_coding) %>%
  distinct(fission_yeast_BP_ID, .keep_all = T)

down_go_coding <- bind_rows(down4_go_coding,down2_go_coding) %>%
  distinct(fission_yeast_BP_ID, .keep_all = T)


up2_go_coding$geneID_symbol <- up2_go_coding$geneID
up4_go_coding$geneID_symbol <- up4_go_coding$geneID
down2_go_coding$geneID_symbol <- down2_go_coding$geneID
down4_go_coding$geneID_symbol <- down4_go_coding$geneID
down_go_coding$geneID_symbol <- down_go_coding$geneID

up4_go_coding$group  <- "Activated"
down_go_coding$group <- "Repressed"



#############################################################################
# Select relevant columns from DESeq2 expression results spreadsheets                                 
#############################################################################

expression_data_16h <- results_annotated %>%
  dplyr::select(c(2,9,3:8,11:14))

expression_data_16h_filtered <-  filtered_annotated %>%
  dplyr::select(c(2,9,3:8,11:14))


#############################################
# Unnest the genes / gene names                                
#############################################
expand_gene_data <- function(df, results_data) {
  df %>%
    mutate(
      gene_name = strsplit(geneID, "/"),
      BP_ID = trimws(fission_yeast_BP_ID)
    ) %>%
    unnest(gene_name) %>%
    select(
      BP_ID,
      Description,
      gene_name
    ) %>%
    # Optional: trim whitespace
    mutate(
      gene_name = trimws(gene_name)
    ) %>%
    left_join(results_data, by = "gene_name") %>%  #First join with results data
    mutate(FoldChange = 2^(log2FoldChange),
           Gene_ID = sub("C$", "c", Gene_ID),
           Gene_ID = trimws(Gene_ID)) %>%
    select(c(2,1,4,3,11,19,6,5,10))
}

up2_go_coding_genes <- up2_go_coding %>% expand_gene_data(results_annotated)
up4_go_coding_genes <- up4_go_coding %>% expand_gene_data(results_annotated) 
down2_go_coding_genes <- down2_go_coding %>% expand_gene_data(results_annotated)
down4_go_coding_genes <- down4_go_coding %>% expand_gene_data(results_annotated)


##################################################
# Group / Annotate / add normal FC info (not log2)                               
##################################################
# summarize and annotate
summarize_and_annotate_GO <- function(data, results_data) {
  data %>%
    group_by(gene_name) %>% rename(terms = Description) %>%
    summarise(BP_ID = paste(BP_ID, collapse = ", "),
              terms = paste(terms, collapse = ", "), # Recombine names into a single string
              .groups = "drop") %>%
    left_join(results_data, by = "gene_name") %>%
    mutate(count = str_count(terms, ",") + 1,
           FoldChange = 2^(log2FoldChange),
           Gene_ID = sub("C$", "c", Gene_ID)) %>% 
    select(c(4,1,11,20,6,5,10,19,2,3)) %>% 
    arrange(padj)  # Order by increasing padj
}

up2_go_coding_summanno <- up2_go_coding_genes %>% summarize_and_annotate_GO(results_annotated)
up4_go_coding_summanno <- up4_go_coding_genes %>% summarize_and_annotate_GO(results_annotated)
down2_go_coding_summanno <- down2_go_coding_genes %>% summarize_and_annotate_GO(results_annotated)
down4_go_coding_summanno <- down4_go_coding_genes %>% summarize_and_annotate_GO(results_annotated)


##################################################
# Save data!                                
##################################################

outputFolder <- "~/Rstudio/expression_analysis_16h/results/GO_analysis/genekitr/coding_only/results/"

write_go_term_sheets <- function(original_df, summary_df, expanded_df, output_filename = file) {
  # Create a new workbook
  wb <- createWorkbook()
  
  # Get unique GO Terms
  go_terms <- original_df$NA_BP_ID[!duplicated(original_df$NA_BP_ID)]
  
  # Add original unexpanded dataframe as first sheet
  addWorksheet(wb, "GO Results")
  writeData(wb, "GO Results", original_df)
  
  # Add summary dataframe as second sheet
  addWorksheet(wb, "highlighted genes")
  writeData(wb, "highlighted genes", summary_df)
  
  # Create a subset for each GO Term in the expanded dataframe
  go_term_subsets <- split(expanded_df, expanded_df$BP_ID) %>%
    lapply(function(df) df %>% arrange(padj))  # Sort by increasing padj
  
  # Add expanded dataframe as third sheet
  addWorksheet(wb, "GO Results by gene")
  writeData(wb, "GO Results by gene", expanded_df)
  
  # Write each GO Term subset to a separate sheet
  for (term in names(go_term_subsets)) {
    # Create a safe sheet name (remove invalid characters)
    safe_term <- gsub(":", "_", term)
    
    # Add worksheet
    addWorksheet(wb, safe_term)
    
    # Write data to worksheet
    writeData(wb, safe_term, go_term_subsets[[term]])
  }
  
  # Reorder sheets to match the original GO Term order
  sheet_names <- c("GO Results","highlighted genes","GO Results by gene", gsub(":", "_", go_terms))
  worksheetOrder(wb) <- match(sheet_names, wb$sheet_names)
  
  # Save the workbook
  saveWorkbook(wb, output_filename, overwrite = TRUE)
  
  return(wb)
}

# # 2-Fold Upregulated Excel sheet:
#  write_go_term_sheets(up2_go_coding,up2_go_coding_summanno, up2_go_coding_genes, paste0(outputFolder,"GO_subsets_2fold_upregulated.xlsx"))
# # 4-Fold Upregulated Excel sheet:
#  write_go_term_sheets(up4_go_coding,up4_go_coding_summanno, up4_go_coding_genes, paste0(outputFolder,"GO_subsets_4fold_upregulated.xlsx"))
#  # 2-Fold Downregulated Excel sheet:
#  write_go_term_sheets(down2_go_coding,down2_go_coding_summanno, down2_go_coding_genes, paste0(outputFolder,"GO_subsets_2fold_downregulated.xlsx"))
#  # 4-Fold Downregulated Excel sheet:
#  write_go_term_sheets(down4_go_coding,down4_go_coding_summanno, down4_go_coding_genes, paste0(outputFolder,"GO_subsets_4fold_downregulated.xlsx"))


##################################################
# Making Plots                                
##################################################

# Gene Ontology ORA plots -------------------------------------------------


library(readxl)
library(dplyr)
library(ggplot2)
library(stringr)
library(DOSE)

# ---- Read sheet 1 of each file ----
up   <- read_excel("~/Rstudio/expression_analysis_16h/results/GO_analysis/genekitr/coding_only/results/GO_subsets_4fold_upregulated.xlsx",   sheet = 1)
down <- read_excel("~/Rstudio/expression_analysis_16h/results/GO_analysis/genekitr/coding_only/results/GO_subsets_2fold_downregulated.xlsx", sheet = 1)

# ---- Redundancy flagging helper ----
# Flags a term as redundant if >= overlap_threshold of its genes are already
# covered by a higher-ranked (kept) term. Walk order = the row order passed in,
# so arrange BEFORE calling this.

flag_redundant <- function(df, gene_col = "geneID", overlap_threshold = 0.5) {
  gene_sets <- str_split(df[[gene_col]], "/")
  redundant      <- logical(nrow(df))
  redundant_with <- rep(NA_character_, nrow(df))
  kept <- integer(0)
  for (i in seq_len(nrow(df))) {
    gi <- gene_sets[[i]]
    flagged <- FALSE
    for (k in kept) {
      overlap <- length(intersect(gi, gene_sets[[k]])) / length(gi)
      if (overlap >= overlap_threshold) {
        redundant[i]      <- TRUE
        redundant_with[i] <- df$Description[k]
        flagged <- TRUE
        break
      }
    }
    if (!flagged) kept <- c(kept, i)
  }
  df$redundant      <- redundant
  df$redundant_with <- redundant_with
  df
}

# ---- Filter and order the UPregulated set ----
up_filtered <- up %>%
  mutate(
    bg_num      = as.numeric(sub("/.*", "", BgRatio)),
    bg_den      = as.numeric(sub(".*/", "", BgRatio)),
    bg_fraction = bg_num / bg_den
  ) %>%
  filter(
    FoldEnrich  > 1.5,
    bg_fraction < 0.20,
    p.adjust    < 0.05
  ) %>%
  arrange(bg_fraction) %>%        # most specific (smallest background) first
  flag_redundant(overlap_threshold = 0.5) %>%
  arrange(p.adjust, desc(FoldEnrich))   # re-sort for display/plotting after flagging

# ---- Downregulated: no filtering, same ordering ----
down_kept <- down %>%
  mutate(
    bg_num      = as.numeric(sub("/.*", "", BgRatio)),
    bg_den      = as.numeric(sub(".*/", "", BgRatio)),
    bg_fraction = bg_num / bg_den
  ) %>%
  filter(p.adjust < 0.05) %>%
  arrange(p.adjust, desc(FoldEnrich)) #%>%
#flag_redundant(overlap_threshold = 0.5)

# ---- Inspect, including the redundancy flag ----
print(up_filtered %>% select(Description, FoldEnrich, bg_fraction, RichFactor, Count, p.adjust, redundant, redundant_with), n = Inf)

# ---- The non-redundant sets to be used for plotting ----
up_nonredundant <- up_filtered %>% 
  filter(
    # Override to keep this specific description
    (!redundant | Description == "regulatory ncRNA-mediated heterochromatin formation") &
      # Override to discard this uninformative process 
      Description != "negative regulation of cellular process"
  )
down_nonredundant <- down_kept   %>% filter(!redundant)



library(dplyr)
library(ggplot2)
library(stringr)

go_combined <- bind_rows(
  up_nonredundant   %>% mutate(direction = "Upregulated (\u22654-fold)"),
  down_kept %>% mutate(direction = "Downregulated (\u22652-fold)")
) %>%
  mutate(
    direction = factor(direction, levels = c("Upregulated (\u22654-fold)",
                                             "Downregulated (\u22652-fold)")),
    
    # order so up terms group together above down terms, each block by RichFactor
    Description = factor(Description,
                         levels = unique(Description[order(direction, RichFactor)]))
  ) 

levels(go_combined$Description) <- levels(go_combined$Description) %>%
  gsub("positive regulation of transcription by RNA polymerase II", 
       "positive regulation of transcription\nby RNA polymerase II", 
       .) %>%
  gsub("organonitrogen compound biosynthetic process", 
       "organonitrogen compound\nbiosynthetic process", 
       .)

ggplot(go_combined, aes(x = RichFactor, y = Description, fill = p.adjust, size = Count)) +
  geom_point(shape = 21, color = "black", stroke = 1) +
  scale_size(range = c(4, 12), name = "Count", breaks = scales::pretty_breaks(n = 3)) + 
  scale_fill_viridis_c(option = "plasma", guide = "colourbar", name = "p.adjust", direction = -1) + 
  #scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) +
  theme_dose()+
  facet_wrap(~ direction, ncol = 1, scales = "free_y") + 
  labs(x = "RichFactor", y = NULL, title = "Gene Ontology ORA Results",
       #subtitle = bquote(italic("points ordered by fold enrichment"))) +
       subtitle = "<i>nup211-so </i>(16h) vs. wild-type (+T) \n") +
  
  xlim(0.01, 0.81) +
  theme( text = element_text(family = "Arial"),
         panel.border = element_rect(colour = "black", fill = NA, linewidth = 2), 
         strip.background = element_rect(fill = "black", colour = "black"), 
         strip.text = element_text(size = 16, face = "bold", colour = "white"),
         axis.text.y = element_text(size = 14, face = "bold"),
         axis.text.x = element_text(size = 16, face = "bold"),
         axis.title.x = element_text(size = 16, face = "bold",
                                     margin = margin(t = 10)),
         legend.text = element_text(size = 14),
         legend.title = element_text(size = 14, face = "bold", family = "Arial"),
         plot.title = element_text(size = 16, face = "bold"),
         plot.subtitle = ggtext::element_markdown(size = 14, face = "plain",
                                                  margin = margin(b = 10)),
         legend.position = "right" 
  )

