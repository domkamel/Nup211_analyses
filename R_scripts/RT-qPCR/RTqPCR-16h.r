#################################################
# RNA-seq Analysis Pipeline (16h) - Script 2
# Purpose: Creating Plots / Data Visualization
#################################################
# Clear workspace
# rm(list=ls())
#################################################
# Setup and Library Loading
###########################################
# Load required libraries
libraries <- c(
  "readr", "readxl", "dplyr", "tidyr", "tidyverse","pcr",
  "ggplot2","ggpubr", "ggExtra" ,"plotly", "gplots", "Glimma",
  "pheatmap", "RColorBrewer", "viridis", "PoiClaClu", 
  "kableExtra", "pcr","broom","rstatix", "ggsignif"
)

for (lib in libraries) {
  if (!require(lib, character.only = TRUE)) {
    BiocManager::install(lib)
    library(lib, character.only = TRUE)
  }
}

# devtools::install_github('MahShaaban/pcr')

# Gene display order used throughout this script
gene_levels <- c("chp1", "mit1", "mus7", "for3", "tea2", "pho1", "urg2")
###########################################
# Data Import
###########################################

top_transcripts_16h <- read_excel("~/Desktop/qPCR data/16h RT-qPCR/CT_values_16h_TopTranscripts.xlsx")

top_transcripts_16h <- top_transcripts_16h %>% 
  mutate(condition = ifelse(grepl("558_T", sample), "control", "shutoff_16h"),
         bio_rep = as.numeric(factor(gsub(".*_(\\d)$", "\\1", sample))))

####################################################
# Compute average Ct for each biological replicate
####################################################
tech_rep_avg_16h <- top_transcripts_16h %>%
  group_by(gene, condition, bio_rep) %>%
  summarise(mean_Ct = mean(Ct), .groups = "drop") %>%
  arrange(gene, condition, bio_rep)

###############################################################
# Put data in format for Prism if licence available (optional) 
###############################################################
#  topTranscripts_16h_prism <-tech_rep_avg_16h %>%
#   mutate(row_id = paste(condition, bio_rep, sep = "_")) %>%  # Create new column IDs
#   select(gene, row_id, mean_Ct) %>%  # Keep relevant columns
#   pivot_wider(names_from = row_id, values_from = mean_Ct)

#   write_csv(topTranscripts_16h_prism, file = "~/Desktop/qPCR data/16h RT-qPCR/topTranscripts_16h_Ct_prism.csv")
####################################################
# Format data for getting individual replicate data 
####################################################
topTranscripts_16h_Ct <- tech_rep_avg_16h %>%
  mutate(row_id = paste(condition, bio_rep, sep = "_")) %>%  # Create new row ID
  select(row_id, gene, mean_Ct) %>%  # Keep relevant columns
  pivot_wider(names_from = gene, values_from = mean_Ct)

####################################################
# Format data for use with the PCR R program 
####################################################
topTranscripts_16h_wide <- tech_rep_avg_16h %>%
  mutate(row_id = paste(condition, bio_rep, sep = "_")) %>%  # Create new row ID
  select(row_id, gene, mean_Ct) %>%  # Keep relevant columns
  pivot_wider(names_from = gene, values_from = mean_Ct) %>%
  select(-row_id) # blank this out in case you want to see if the Ct values are in the right order

# **** First Approach: Using the PCR program   **** #
#####################################################
# Run the PCR R program 
#####################################################
# Provide a grouping variable:
  group <- rep(c('control', 'shutoff_16h'), each = 3)

library('pcr')
# Get expression information
  topTranscripts_16h_FC <-  pcr_ddct(topTranscripts_16h_wide,
                                      group_var = group,
                                      reference_gene = 'act1',
                                      reference_group = 'control',
                                      mode = "separate_tube",
                                      plot = F)
  
# Get significance information 
  topTranscripts_16h_stats <- pcr_ttest(topTranscripts_16h_wide,
                                         group_var = group,
                                         reference_gene = 'act1',
                                         reference_group = 'control',
                                         tidy = T)
# Correct for multiple hypothesis testing
# (BH adjustment):
  topTranscripts_16h_stats <- topTranscripts_16h_stats %>% 
                                     mutate(padj = p.adjust(p_value, method = "BH"),
                                            significant = case_when(
                                                          padj < 0.001 ~ "***",
                                                          padj < 0.01  ~ "**",
                                                          padj < 0.05  ~ "*",
                                                          TRUE          ~ "ns"
                                                        )
                                            )
  
# Merge significance and expression information 
  topTranscripts_16h_data <- topTranscripts_16h_FC %>% 
                              left_join(topTranscripts_16h_stats %>% 
                                        select(gene, estimate, p_value, padj, significant), 
                                        by = "gene",
                                        keep = F) %>%
                              rename(FoldChange = relative_expression,
                                     pvalue = p_value,
                                     condition = group) %>%
                              mutate(log2FoldChange = log2(FoldChange),
                                     log2Upper = log2(upper),  # Log2 transform upper value directly
                                     log2Lower = log2(lower))  # Log2 transform lower value directly

#################################################################
# **** Second Approach: Using manual/prism calculations   ***** #
#################################################################

  # Calculate delta Cts for each transcript per condition 
  ###############################################################
  
  deltaCt_df <- topTranscripts_16h_Ct %>%
    pivot_longer(cols = -row_id, names_to = "gene", values_to = "Ct") %>%
    group_by(row_id) %>%
    mutate(act1_Ct = Ct[gene == "act1"]) %>%  # Store Act1 Ct for each sample
#    filter(!(gene %in% c("act1","adh1","urg1","gcd1"))) %>%  # Remove Act1 from final output
    mutate(deltaCt = Ct - act1_Ct) %>% 
    ungroup()

  # Calculate the mean of the deltaCts for only the control condition
  ###################################################################    
  mean_deltaCt_controls <- deltaCt_df %>%
    filter(grepl("control", row_id)) %>%
    group_by(gene) %>%
    summarise(mean_deltaCt_control = mean(deltaCt, na.rm = TRUE))
  
  # Merge the dataframes and calculate DeltaDeltaCt for each gene/condition
  #########################################################################   
  df_ddCt <- deltaCt_df %>%
    left_join(mean_deltaCt_controls, by = "gene") %>%
    mutate(DeltaDeltaCt = deltaCt - mean_deltaCt_control)
  
  # Then calculate FoldChange (2^-DeltaDeltaCt) for each gene/condition
  # and calculate log2FoldChange
  ######################################################################### 
  
  df_FC <- df_ddCt %>% mutate(FoldChange = 2^(-DeltaDeltaCt),
                              log2FoldChange = log2(FoldChange),
                              condition = case_when(
                                grepl("control", row_id, ignore.case = TRUE) ~ "control",
                                grepl("shutoff_16", row_id, ignore.case = TRUE) ~ "shutoff_16h",
                                TRUE ~ "other" # In case there are unexpected values
                              ))
  ###############################################################
  # Put data in format for Prism if licence available (optional) 
  ###############################################################
    df_FC_prism <-df_FC %>%
          select(gene, row_id, FoldChange) %>%  # Keep relevant columns
          pivot_wider(names_from = row_id, values_from = FoldChange)
  
#  write_csv(df_FC_prism, file = "~/Desktop/qPCR data/16h RT-qPCR/topTranscripts_16h_FC_prism.csv")
  
  ################################################################################ 
  # Now, that we have the expression values for each replicate, we can calculate
  # the  average/mean values for each transcript in each condition and check them
  # against the pcr R program output to verify everything was done correctly 
  #          Both a sanity check and useful for plotting results
  ################################################################################ 
  mean_values <- df_FC %>%
                  group_by(gene, condition) %>%
                  summarise(mean_log2FC = mean(log2FoldChange, na.rm = TRUE), 
                            mean_FC = 2^(mean_log2FC), 
                            .groups = "drop") %>%
                  left_join(
                    topTranscripts_16h_data %>%
                      select(gene, estimate, pvalue, padj, significant, upper, lower, log2Upper, log2Lower, condition) %>%
                      distinct(),  # Ensure uniqueness before joining
                    by = c("gene", "condition")
                  )

  # After verifying this approach, I could also calculate pvalues and padj's 
  # independently of the pcr R program just to show either data can be used for 
  # plotting significance! 
  #     Again, a sanity check and potentially useful for plotting results
  ################################################################################ 
  t_test_results <- df_FC %>%
                      group_by(gene) %>%
                      summarise(
                        t_test = list(t.test(log2FoldChange ~ condition, 
                                             data = cur_data()))) %>%
                      mutate(tidy_results = map(t_test, broom::tidy)) %>%
                      unnest(tidy_results) %>%
                      select(gene, estimate, conf.low, conf.high, p.value) %>%
                      mutate(FoldChange = 2^((-1)*(estimate)),
                             log2FoldChange = log2(FoldChange),
                             group1 = "control", 
                             group2 = "shutoff_16h",
                             condition = "shutoff_16h",
                             padj = p.adjust(p.value, method = "BH"),
                             significant = case_when(
                               padj < 0.001 ~ "***",
                               padj < 0.01  ~ "**",
                               padj < 0.05  ~ "*",
                               TRUE          ~ "ns")
                             ) %>% group_by(gene) %>%
                      mutate(y.position = max(df_FC$log2FoldChange[df_FC$gene == gene]) + 0.5) %>%
                      ungroup()
  ###############################################################
  # Put data in format for Prism if licence available (optional) 
  ###############################################################
  #  topTranscripts_16h_prism <-tech_rep_avg_16h %>%
  #   mutate(row_id = paste(condition, bio_rep, sep = "_")) %>%  # Create new column IDs
  #   select(gene, row_id, mean_Ct) %>%  # Keep relevant columns
  #   pivot_wider(names_from = row_id, values_from = mean_Ct)
  
  #   write_csv(topTranscripts_16h_prism, file = "~/Desktop/qPCR data/16h RT-qPCR/topTranscripts_16h_Ct_prism.csv")

 shutoff_FC <- df_FC %>%
   filter(condition == "shutoff_16h",
          !gene %in% c("act1","adh1","urg1","gcd1","clp1","cdc12"))  
 
 shutoff_datapoints <- topTranscripts_16h_data %>%
   filter(condition == "shutoff_16h",
          !gene %in% c("act1","adh1","urg1","gcd1","clp1","cdc12")) 
 
 shutoff_mean_values_16h <- mean_values %>%
   mutate(log2FoldChange = mean_log2FC,
     fill_color = case_when(
     mean_log2FC > 1 & padj < 0.5  ~ "salmon",
     mean_log2FC < -1 & padj < 0.5 ~ "skyblue",
     TRUE ~ "gray80"  # optional: neutral color for values in between
   )) %>%
   filter(condition == "shutoff_16h",
          !gene %in% c("act1","adh1","urg1","gcd1","clp1","cdc12")) 
 
 t_test_signif <- t_test_results %>%
   filter(!gene %in% c("act1","adh1","urg1","gcd1","clp1","cdc12")) 
 
 shutoff_FC_16h <- shutoff_FC %>%
   mutate(gene = factor(gene, levels = gene_levels)) %>%
   mutate(fill_color = case_when(
     log2FoldChange > 1  ~ "salmon",
     log2FoldChange < -1 ~ "skyblue",
     TRUE ~ "gray80"  # optional: neutral color for values in between
   ))
 
 shutoff_datapoints <- shutoff_datapoints %>%
   mutate(gene = factor(gene, levels = gene_levels))  %>%
   mutate(fill_color = case_when(
     log2FoldChange > 1  ~ "salmon",
     log2FoldChange < -1 ~ "skyblue",
     TRUE ~ "gray80"  # optional: neutral color for values in between
   ))
 
 t_test_signif <- t_test_signif %>%
   mutate(gene = factor(gene, levels = gene_levels)) %>%
   mutate(fill_color = case_when(
     log2FoldChange > 1  ~ "salmon",
     log2FoldChange < -1 ~ "skyblue",
     TRUE ~ "gray80"  # optional: neutral color for values in between
   ))

 # FINAL 16h figure -- version 1: error bars & stars sourced from shutoff_mean_values_16h 
 ggplot(shutoff_FC_16h, aes(x = gene, y = log2FoldChange, fill = fill_color)) +
   geom_errorbar(data = shutoff_mean_values_16h, aes(x = gene, ymin = log2Lower, ymax = log2Upper), width = 0.2, size = 1) +
   geom_boxplot(alpha = 1, position = "dodge2", coef = 0, outlier.shape = NA, size = 0.8, color = "black") +
   #geom_point(size = 2, alpha = 0.7) +
   geom_jitter(width = 0.28, size = 3, alpha = 0.9, shape = 21, fill = "white", color = "black", stroke = 1.3) +
   # geom_text(data = shutoff_mean_values_16h, aes(x = gene, y = log2Upper + 0.7, 
   #                                           label = round(mean_log2FC, 2)), position = position_dodge(0.9),
   #           size = 5, fontface = "bold") +
   geom_text(data = shutoff_mean_values_16h, aes(x = gene, y = log2Upper + 0.2, #max(df_FC$log2FoldChange) - 0.1, 
                                            label = paste(significant)), position = position_dodge(0.9),
             size = 12, fontface = "bold") +
   geom_hline(yintercept = 0, linewidth = 1.0, linetype = "solid", color = "black") +
   theme_gray() +
   labs(title = "RT-qPCR Results (16h)",
        subtitle = bquote(bolditalic("nup211-so") ~ bold("vs. wild-type (+T)") ~ "\n"),
        legend = "",
        y = bquote(bold(log[2]~ 'Fold Change')),
        x = "") +
   scale_y_continuous(limits = c(-5,4),
                      breaks=seq(-5,4, 1))+
   # scale_fill_manual(name = "",
   #                   values = c("shutoff_16h" = "salmon"), 
   #                   labels = c(bquote(bolditalic("nup211-so") ~ bold("(16h)")))) +
   theme(plot.title = element_text(size = 18, face = "bold"),
         plot.subtitle = element_text(size = 16, face = "bold"),
         axis.line = element_line(size = 1.0), # Control axis thickness
         axis.ticks = element_line(size = 1.2), # Control tick mark thickness
         axis.ticks.length = unit(3, "mm"), # Control tick mark length
         axis.text.x = element_text(size = 18, face= "bold",colour = 'black',
                                    angle = 30, hjust = 1),    # Axis labels
         axis.text.y = element_text(size = 16, face= "bold",colour = 'black'),
         axis.title.y = element_text(size = 18, face= "bold",colour = 'black'),   # Axis titles
         legend.title = element_text(size = 14, face = "bold"),   # Adjust legend title font
         legend.text = element_text(size = 14, face = "bold" ), # Adjust legend labels font
         legend.position = "none") # Adjust legend labels font

 # Save shutoff mean values and shutoff FC 16 hour data

# write_csv(x = shutoff_mean_values_16h, file = "~/Desktop/qPCR data/16h RT-qPCR/shutoff_mean_values_16h.csv")
# write_csv(x = shutoff_FC_16h, file = "~/Desktop/qPCR data/16h RT-qPCR/shutoff_FC_16h.csv")


 # Then read in the corresponding data for 48 hours  (and eventually the rescue as well!)
 
#  If Available: Read in 48h qPCR analysis of same transcripts -------------------------------------------------
#
 shutoff_FC_48h <- read_csv("~/Desktop/qPCR data/48h RT-qPCR/shutoff_FC_48h.csv")
 shutoff_mean_values_48h <- read_csv("~/Desktop/qPCR data/16h RT-qPCR/shutoff_mean_values_48h.csv")%>%
   mutate(FoldChange = mean_FC)

#  # Merge the replicate data
 merged_FC_data <- shutoff_FC_16h %>%
   select(row_id, gene, FoldChange, log2FoldChange, condition, fill_color) %>%
   bind_rows(
     shutoff_FC_48h %>%
       select(row_id, gene, FoldChange, log2FoldChange, condition, fill_color)
   ) %>% as.data.frame()

 merged_FC_data$fill_color <- factor(merged_FC_data$fill_color,
                                     levels = c("salmon","gray80",  "skyblue"))

#  # Merge the mean values and the statistical info for condition vs. wild-type
 merged_mean_values <- shutoff_mean_values_16h %>%
   select(gene, condition, mean_FC, log2FoldChange, pvalue, padj, significant,
          upper, lower, log2Upper, log2Lower, fill_color) %>%
   bind_rows(
     shutoff_mean_values_48h %>%
       select(gene, condition, mean_FC, log2FoldChange, pvalue, padj, significant,
          upper, lower, log2Upper, log2Lower, fill_color)
     ) %>%
   mutate(significant = replace_na(significant, ""))   %>% 
       rename(FoldChange = mean_FC) %>% 
       as.data.frame() 

#  # Calculate significant differences between the two groups(16 vs 48h)
#  
 t_test_results_16h_48h <- merged_FC_data %>%
   filter(condition %in% c("shutoff_16h", "shutoff_48h")) %>%
   group_by(gene) %>%
   summarise(
     t_test = list(t.test(log2FoldChange ~ condition, data = cur_data()))
   ) %>%
   mutate(tidy_results = map(t_test, broom::tidy)) %>%
   unnest(tidy_results) %>%
   select(gene, estimate1, estimate2, conf.low, conf.high, p.value) %>%
   rename(log2FC_16h = estimate1,
          log2FC_48h = estimate2) %>%
   mutate(group1 = "shutoff_16h",
          group2 = "shutoff_48h",
          comparison = "16h vs 48h",
          padj = p.adjust(p.value, method = "BH"),
          significant = case_when(
            padj < 0.001 ~ "***",
            padj < 0.01  ~ "**",
            padj < 0.05  ~ "*",
            TRUE         ~ "ns")
   ) %>%
   group_by(gene) %>%
   mutate(y.position = max(merged_mean_values$log2FoldChange[merged_mean_values$gene == gene]) + 1.1) %>%
   ungroup()

 t_test_results_16h_48h$y.position <-  c(3.5, 3.75, 3,3.27,3.2,-1.3,-2.92)

# # Then Time to PLOT! 

 p <- ggplot(merged_FC_data, aes(x = gene, y = log2FoldChange, fill = fill_color)) +
   geom_hline(yintercept = 0, linewidth = 1.0, linetype = "solid", color = "black") +

   geom_errorbar(data = merged_mean_values,
                 aes(x = gene, ymin = log2Lower, ymax = log2Upper, group = condition),
                 position = position_dodge(width = 0.6),
                 width = 0.3, size = 1, inherit.aes = FALSE) +

   geom_boxplot(aes(group = interaction(gene, condition)),
                position = position_dodge(width = 0.6),
                coef = 0, outlier.shape = NA,
                size = 0.8, color = "black", alpha = 1) +

   geom_jitter(aes(shape = condition),
               position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5),
               size = 2.5, stroke = 1.2, fill = "white", color = "black") +

   # geom_text(data = merged_mean_values,
   #           aes(x = gene, y = log2Upper + 0.6, label = round(log2FoldChange, 2), group = condition),
   #           position = position_dodge(0.6),
   #           size = 5, fontface = "bold", inherit.aes = FALSE) +

   geom_text(data = merged_mean_values,
             aes(x = gene, y = log2Upper + 0.1, label = paste(significant), group = condition),
             position = position_dodge(0.6),
             size = 8, fontface = "bold", inherit.aes = FALSE) +

   scale_shape_manual(values = c("shutoff_16h" = 21, "shutoff_48h" = 24),
                      labels = c("16 hours",
                                 "48 hours")) +
   scale_fill_manual(values = c("salmon" = "salmon",
                                "gray80" = "gray80",
                                "skyblue" = "skyblue"),
                     labels = c("Upregulated",
                                "NS",
                                "Downregulated")) +

   scale_y_continuous(limits = c(-5.2, 4.03), breaks = seq(-5, 4, 1)) +

   labs(title = "RT-qPCR Results (16h & 48h)",
        subtitle = bquote(bolditalic("nup211-so") ~ bold("vs. wild-type (+T)")),
        x = NULL,
        y = bquote(bold(log[2] ~ "Fold Change")),
        shape = "Timepoint",
        fill = "Expression level") +

   theme(
     plot.title = element_text(size = 18, face = "bold"),
     plot.subtitle = element_text(size = 16, face = "bold"),
     axis.line = element_line(size = 1.0),
     axis.ticks = element_line(size = 1.2),
     axis.ticks.length = unit(3, "mm"),
     axis.text.x = element_text(size = 18, face = "bold", colour = 'black', angle = 30, hjust = 1),
     axis.text.y = element_text(size = 16, face = "bold", colour = 'black'),
     axis.title.y = element_text(size = 18, face = "bold", colour = 'black'),
     legend.title = element_text(size = 14, face = "bold"),
     legend.text = element_text(size = 14, face = "bold"),
     legend.position = "right"
   )

print (p)

p +
  geom_signif(
    data = t_test_results_16h_48h,
    aes(
      x = gene,
      xmin = as.numeric(factor(gene)) - 0.15,  # left position (16h)
      xmax = as.numeric(factor(gene)) + 0.15,  # right position (48h)
      annotations = significant,
      fontface = "bold",
      y_position = y.position
    ),
    tip_length = 0.01,
    size = 1,
    textsize = 8,
    vjust = 0.2,
    inherit.aes = F,
    manual = T
  )
 
# merged_FC_data
# merged_mean_values
# t_test_results_16h_48h
# 
# replicate_data_table <-  merged_FC_data %>% select(c(1:4)) %>% pivot_wider(names_from = row_id, values_from = c(FoldChange,log2FoldChange)) %>%
#                   `colnames<-`(c("gene", "16h_rep1_FC","16h_rep2_FC","16h_rep3_FC","48h_rep1_FC","48h_rep2_FC","48h_rep3_FC",
#                                  "16h_rep1_log2FC","16h_rep2_log2FC","16h_rep3_log2FC", "48h_rep1_log2FC","48h_rep2_log2FC","48h_rep3_log2FC"))
# 
# mean_data_table <-  merged_mean_values  %>% select(-12) %>% pivot_wider(names_from = condition, values_from = c(FoldChange,log2FoldChange, pvalue, padj, 
#                                                                                                significant, upper,lower, log2Upper, log2Lower))  %>%
#   `colnames<-`(c("gene", "mean_FC_16h","mean_FC_48h","mean_log2FC_16h","mean_log2FC_48h","pvalue_16h","pvalue_48h","padj_16h","padj_48h",
#                  "annotation_16h","annotation_48h","upperCI_16h","upperCI_48h","lowerCI_16h","lowerCI_48h","log2upperCI_16h","log2upperCI_48h",
#                  "log2lowerCI_16h","log2lowerCI_48h")) %>%
#   mutate(annotation_48h = ifelse(annotation_48h == "", "ns", annotation_48h)) %>% merge(t_test_results_16h_48h %>% 
#                                                                                   select(c(gene, p.value,padj,significant)), by = 'gene') %>%
#   rename(`16h_vs_48h_pvalue`= p.value,
#          `16h_vs_48h_padj`= padj,
#          `16h_vs_48h_annotation`= significant)
# 
# supplementary_data <- replicate_data_table %>% left_join(mean_data_table, by = 'gene')
# 
# supplemental_table <-  supplementary_data %>% select(c(1,8:10,11:13,16,28,30,17,29,31,20,22,21,23,33,34)) %>%
#                       rename(padj_16h_v_control = 14,
#                              annotation_16h_v_control = 15,
#                              padj_48h_v_control = 16,
#                              annotation_48h_v_control = 17,
#                              padj_16h_vs_48h = 18,
#                              annotation_16h_vs_48h = 19)
# 
# #write_csv(x = supplemental_table, file = "~/Desktop/qPCR data/16h RT-qPCR/supplemental_table.csv")
# 
# normal_scale_supp <- supplementary_data%>% select(c(1,2:7,24,26,14,25,27,15,20,22,21,23,33,34))
