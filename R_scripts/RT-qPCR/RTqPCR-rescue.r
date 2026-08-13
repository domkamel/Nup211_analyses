#################################################
# RNA-seq Analysis Pipeline (rescue) - Script 4
# Purpose: Creating Plots / Data Visualization
#################################################
# Clear workspace
# rm(list=ls())
#################################################
# Setup and Library Loading
###########################################
# Load required libraries
libraries <- c(
  "readr", "readxl", "dplyr", "tidyr", "tidyverse",
  "ggplot2","ggpubr", "ggExtra" ,"plotly", "gplots", "Glimma",
  "pheatmap", "RColorBrewer", "viridis", "PoiClaClu", 
  "kableExtra", "pcr","broom","rstatix", "ggsignif", "ggtext"
)

for (lib in libraries) {
  if (!require(lib, character.only = TRUE)) {
    BiocManager::install(lib)
    library(lib, character.only = TRUE)
  }
}

# Gene display order used throughout this script
gene_levels <- c("chp1", "mit1", "mus7", "for3", "tea2", "pho1", "urg2")
###########################################
# Data Import
###########################################
top_transcripts_rescue_1 <- read_excel("~/Desktop/qPCR data/Rescue RT-qPCR/CT_values_TopTranscripts_fullLength_rescue_1.xlsx")
top_transcripts_rescue_2 <- read_excel("~/Desktop/qPCR data/Rescue RT-qPCR/CT_values_TopTranscripts_fullLength_rescue_2.xlsx")

top_transcripts_rescue_1 <- top_transcripts_rescue_1 %>% 
  mutate(condition = ifelse(grepl("control", sample), "control", "FL_rescue"),
         bio_rep = as.numeric(factor(gsub(".*_(\\d)$", "\\1", sample))))
top_transcripts_rescue_2 <- top_transcripts_rescue_2 %>% 
  mutate(condition = ifelse(grepl("control", sample), "control", "FL_rescue"),
         bio_rep = as.numeric(factor(gsub(".*_(\\d)$", "\\1", sample))))
####################################################
# Compute average Ct for each biological replicate
####################################################
tech_rep_avg_rescue_1 <- top_transcripts_rescue_1 %>%
  group_by(gene, condition, bio_rep) %>%
  summarise(mean_Ct = mean(Ct), .groups = "drop") %>%
  arrange(gene, condition, bio_rep)

tech_rep_avg_rescue_2 <- top_transcripts_rescue_2 %>%
  group_by(gene, condition, bio_rep) %>%
  summarise(mean_Ct = mean(Ct), .groups = "drop") %>%
  arrange(gene, condition, bio_rep)
####################################################
# Format data for getting individual replicate data 
####################################################
topTranscripts_rescue_Ct_1 <- tech_rep_avg_rescue_1 %>%
  mutate(row_id = paste(condition, bio_rep, sep = "_")) %>%  # Create new row ID
  select(row_id, gene, mean_Ct) %>%  # Keep relevant columns
  pivot_wider(names_from = gene, values_from = mean_Ct) 

topTranscripts_rescue_Ct_2 <- tech_rep_avg_rescue_2 %>%
  mutate(row_id = paste(condition, bio_rep, sep = "_")) %>%  # Create new row ID
  select(row_id, gene, mean_Ct) %>%  # Keep relevant columns
  pivot_wider(names_from = gene, values_from = mean_Ct)

####################################################
# Format data for use with the PCR R program 
####################################################
topTranscripts_rescue_wide_1 <- tech_rep_avg_rescue_1 %>%
  mutate(row_id = paste(condition, bio_rep, sep = "_")) %>%  # Create new row ID
  select(row_id, gene, mean_Ct) %>%  # Keep relevant columns
  pivot_wider(names_from = gene, values_from = mean_Ct) %>%
  select(-row_id) # blank this out in case you want to see if the Ct values are in the right order

topTranscripts_rescue_wide_2 <- tech_rep_avg_rescue_2 %>%
  mutate(row_id = paste(condition, bio_rep, sep = "_")) %>%  # Create new row ID
  select(row_id, gene, mean_Ct) %>%  # Keep relevant columns
  pivot_wider(names_from = gene, values_from = mean_Ct) %>%
  select(-row_id) # blank this out in case you want to see if the Ct values are in the right order

                      #***#
#####################################################
# **** First Approach: Using the PCR program   **** #
#####################################################
# Run the PCR R program 
#####################################################
# Provide a grouping variable:
group <- rep(c('FL_rescue','control'), each = 3)

# Get expression information
topTranscripts_rescue_FC_1 <-  pcr_ddct(topTranscripts_rescue_wide_1,
                                     group_var = group,
                                     reference_gene = 'act1',
                                     reference_group = 'control',
                                     mode = "separate_tube",
                                     plot = F)
topTranscripts_rescue_FC_2 <-  pcr_ddct(topTranscripts_rescue_wide_2,
                                     group_var = group,
                                     reference_gene = 'act1',
                                     reference_group = 'control',
                                     mode = "separate_tube",
                                     plot = F)

# Get significance information 
topTranscripts_rescue_stats_1 <- pcr_ttest(topTranscripts_rescue_wide_1,
                                        group_var = group,
                                        reference_gene = 'act1',
                                        reference_group = 'control',
                                        tidy = T)
topTranscripts_rescue_stats_2 <- pcr_ttest(topTranscripts_rescue_wide_2,
                                        group_var = group,
                                        reference_gene = 'act1',
                                        reference_group = 'control',
                                        tidy = T)
# Correct for multiple hypothesis testing
# (BH adjustment):
topTranscripts_rescue_stats_1 <- topTranscripts_rescue_stats_1 %>% 
  mutate(padj = p.adjust(p_value, method = "BH"),
         significant = case_when(
           padj < 0.001 ~ "***",
           padj < 0.01  ~ "**",
           padj < 0.05  ~ "*",
           TRUE          ~ "ns"
         )
  )
topTranscripts_rescue_stats_2 <- topTranscripts_rescue_stats_2 %>% 
  mutate(padj = p.adjust(p_value, method = "BH"),
         significant = case_when(
           padj < 0.001 ~ "***",
           padj < 0.01  ~ "**",
           padj < 0.05  ~ "*",
           TRUE          ~ "ns"
         )
  )
# Merge significance and expression information 
topTranscripts_FLrescue_data_1 <- topTranscripts_rescue_FC_1 %>% 
  left_join(topTranscripts_rescue_stats_1 %>% 
              select(gene, estimate, p_value, padj, significant), 
            by = "gene",
            keep = F) %>%
  rename(FoldChange = relative_expression,
         pvalue = p_value,
         condition = group) %>%
  mutate(log2FoldChange = log2(FoldChange),
         log2Upper = log2(upper),  # Log2 transform upper value directly
         log2Lower = log2(lower))  # Log2 transform lower value directly

topTranscripts_FLrescue_data_2 <- topTranscripts_rescue_FC_2 %>% 
  left_join(topTranscripts_rescue_stats_2 %>% 
              select(gene, estimate, p_value, padj, significant), 
            by = "gene",
            keep = F) %>%
  rename(FoldChange = relative_expression,
         pvalue = p_value,
         condition = group) %>%
  mutate(log2FoldChange = log2(FoldChange),
         log2Upper = log2(upper),  # Log2 transform upper value directly
         log2Lower = log2(lower))  # Log2 transform lower value directly

topTranscripts_FLrescue_data <- bind_rows(topTranscripts_FLrescue_data_1, topTranscripts_FLrescue_data_2)  

# Get confidence interval information

FLrescue_CI_bars <- topTranscripts_FLrescue_data %>% select(gene,condition,upper,lower,log2Upper,log2Lower) %>% 
                    filter(condition == "FL_rescue")

                            #***#
#################################################################
# **** Second Approach: Using manual/prism calculations   ***** #
#################################################################

# Calculate delta Cts for each transcript per condition 
###############################################################
deltaCt_df_1 <- topTranscripts_rescue_Ct_1 %>%
  pivot_longer(cols = -row_id, names_to = "gene", values_to = "Ct") %>%
  group_by(row_id) %>%
  mutate(act1_Ct = Ct[gene == "act1"]) %>%  # Store Act1 Ct for each sample
  #    filter(!(gene %in% c("act1","adh1","urg1","gcd1"))) %>%  # Remove Act1 from final output
  mutate(deltaCt = Ct - act1_Ct) %>% 
  ungroup()

deltaCt_df_2 <- topTranscripts_rescue_Ct_2 %>%
  pivot_longer(cols = -row_id, names_to = "gene", values_to = "Ct") %>%
  group_by(row_id) %>%
  mutate(act1_Ct = Ct[gene == "act1"]) %>%  # Store Act1 Ct for each sample
  #    filter(!(gene %in% c("act1","adh1","urg1","gcd1"))) %>%  # Remove Act1 from final output
  mutate(deltaCt = Ct - act1_Ct) %>% 
  ungroup()
# Calculate the mean of the deltaCts for only the control condition
###################################################################    
mean_deltaCt_controls_1 <- deltaCt_df_1 %>%
  filter(grepl("control", row_id)) %>%
  group_by(gene) %>%
  summarise(mean_deltaCt_control = mean(deltaCt, na.rm = TRUE))

mean_deltaCt_controls_2 <- deltaCt_df_2 %>%
  filter(grepl("control", row_id)) %>%
  group_by(gene) %>%
  summarise(mean_deltaCt_control = mean(deltaCt, na.rm = TRUE))

# Merge the dataframes and calculate DeltaDeltaCt for each gene/condition
#########################################################################   
df_ddCt_1 <- deltaCt_df_1 %>%
  left_join(mean_deltaCt_controls_1, by = "gene") %>%
  mutate(DeltaDeltaCt = deltaCt - mean_deltaCt_control)

df_ddCt_2 <- deltaCt_df_2 %>%
  left_join(mean_deltaCt_controls_2, by = "gene") %>%
  mutate(DeltaDeltaCt = deltaCt - mean_deltaCt_control)

# Then calculate FoldChange (2^-DeltaDeltaCt) for each gene/condition
# and calculate log2FoldChange
######################################################################### 

df_FC_1 <- df_ddCt_1 %>% mutate(FoldChange = 2^(-DeltaDeltaCt),
                                log2FoldChange = log2(FoldChange),
                                condition = case_when(
                                  grepl("control", row_id, ignore.case = TRUE) ~ "control",
                                  grepl("FL_rescue", row_id, ignore.case = TRUE) ~ "FL_rescue",
                                  TRUE ~ "other" # In case there are unexpected values
                                )) %>% 
  filter(gene != "act1")
df_FC_2 <- df_ddCt_2 %>% mutate(FoldChange = 2^(-DeltaDeltaCt),
                                log2FoldChange = log2(FoldChange),
                                condition = case_when(
                                  grepl("control", row_id, ignore.case = TRUE) ~ "control",
                                  grepl("FL_rescue", row_id, ignore.case = TRUE) ~ "FL_rescue",
                                  TRUE ~ "other" # In case there are unexpected values
                                )) %>% 
  filter(gene != "act1")

df_FC <- bind_rows(df_FC_1, df_FC_2)  

# ###############################################################
# # Put data in format for Prism if licence available (optional) 
# ###############################################################
# df_FC_prism <-df_FC %>%
#   select(gene, row_id, FoldChange) %>%  # Keep relevant columns
#   pivot_wider(names_from = row_id, values_from = FoldChange)
# 
# # write_csv(df_FC_prism, file = "~/Desktop/qPCR data/Rescue RT-qPCR/topTranscripts_rescue_FC_prism.csv")

################################################################################ 
# Now, that we have the expression values for each replicate, we can calculate
# the  average/mean values for each transcript in both conditions 
#
#          Both a sanity check and useful for plotting results
################################################################################ 
mean_values <- df_FC %>%
  group_by(gene, condition) %>%
  summarise(log2FoldChange = mean(log2FoldChange, na.rm = TRUE), 
            FoldChange = 2^(log2FoldChange), 
            .groups = "drop") 

# After, I can calculate pvalues and padj's independently of the pcr R program 
# just to show either data can be used for plotting significance! 
#
#  --> Again, a sanity check and potentially useful for plotting results
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
         group2 = "FL_rescue",
         condition = "FL_rescue",
         padj = p.adjust(p.value, method = "BH"),
         significant = case_when(
           padj < 0.001 ~ "***",
           padj < 0.01  ~ "**",
           padj < 0.05  ~ "*",
           TRUE          ~ "ns")
  ) %>% group_by(gene) %>%
  mutate(y.position = max(df_FC$log2FoldChange[df_FC$gene == gene]) + 0.5) %>%
  ungroup()

################################################################################ 
# 
################################################################################ 

rescue_FC <- df_FC %>%
  filter(condition == "FL_rescue") %>%
  mutate(gene = factor(gene, levels = gene_levels)) %>%
  mutate(fill_color = case_when(
    log2FoldChange < -1 ~ "skyblue",
    TRUE ~ "gray80"  # optional: neutral color for values in between
  ))

t_test_signif <- t_test_results %>%
  mutate(gene = factor(gene, levels = gene_levels)) %>%
  mutate(significant = ifelse(significant == "ns", "", significant),
         fill_color = case_when(
           log2FoldChange < -1 ~ "skyblue",
           TRUE ~ "gray80"  # optional: neutral color for values in between
         )) %>% as.data.frame()

rescue_mean_values <- mean_values %>%
  mutate(gene = factor(gene, levels = gene_levels))  %>%
  filter(condition == "FL_rescue") %>%
  left_join(FLrescue_CI_bars,
            by = c("gene", "condition")
  ) %>% select(c(1:2,4:6,3,7:8))  %>%
  left_join(
    t_test_signif %>%
      select(gene, padj, significant, y.position, fill_color) %>%
      distinct(),  # Ensure uniqueness before joining
    by ="gene"
  ) %>% as.data.frame()

###########################################
# Make the plot!
###########################################

ggplot(rescue_FC, aes(x = gene, y = log2FoldChange, fill = fill_color)) +
  geom_hline(yintercept = 0, linewidth = 1.0, linetype = "solid", color = "black", alpha = 0.7) +
  geom_errorbar(data = rescue_mean_values, aes(x = gene, ymin = log2Lower, ymax = log2Upper),
                width = 0.2, size = 1) +
  geom_boxplot(alpha = 1,  coef = 0, outlier.shape = NA, size = 0.8, color = "black") +
  geom_jitter(width = 0.28, size = 3, alpha = 0.9, shape = 21, fill = "white", color = "black", stroke = 1.3) +
  # geom_text(data = rescue_mean_values, aes(x = gene, log2Upper + 0.7,
  #                                           label = round(log2FoldChange, 2)), position = position_dodge(0.9),
  #           size = 6, fontface = "bold") +
  geom_text(data = rescue_mean_values, aes(x = gene, y = log2Upper + 0.2, #max(df_FC$log2FoldChange) - 0.1, 
                                            label = paste(significant)), position = position_dodge(1),
            size = 12, fontface = "bold") +
  theme_gray() +
  labs(title = "RT-qPCR: Full-length Nup211 rescue results",
       subtitle = bquote(italic("nup211-so") ~ "+ pARS1 vs. wild-type (+T)" ~ "\n"),
       legend = "",
       y = bquote(bold(log[2]~ 'Fold Change')),
       x = "",
       fill = "Expression") +
  scale_y_continuous(limits = c(-4.6,1.2),
                     breaks=seq(-4,1, 1)) +
  scale_fill_manual(values = c("gray80" = "lightgray","skyblue" = "skyblue"),
                      labels = c("NS",
                                  "Downregulated")) + 
  theme(plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 14, face = "plain"),
        axis.line = element_line(size = 1.0), # Control axis thickness
        axis.ticks = element_line(size = 1.2), # Control tick mark thickness
        axis.ticks.length = unit(3, "mm"), # Control tick mark length
        axis.text.x = element_text(size = 18, face= "bold",colour = 'black',
                                   angle = 30, hjust = 1),    # Axis labels
        axis.text.y = element_text(size = 16, face= "bold",colour = 'black'),
        axis.title.y = element_text(size = 18, face= "bold",colour = 'black'),   # Axis titles
        legend.title = element_text(size = 14, face = "bold"),   # Adjust legend title font
        legend.text = element_text(size = 14, face = "bold" ), # Adjust legend labels font
        legend.position = "right") # Adjust legend labels font

# write_csv(rescue_FC, file = "Desktop/qPCR data/Rescue RT-qPCR/rescue_FC_value_points.csv")
# write_csv(rescue_mean_values, file = "Desktop/qPCR data/Rescue RT-qPCR/rescue_mean_values_significance.csv" )

###########################################
# Now read in the 16h data!
###########################################

shutoff_FC_16h <- read_csv("~/Desktop/qPCR data/16h RT-qPCR/shutoff_FC_16h.csv") 
shutoff_mean_values_16h <- read_csv("~/Desktop/qPCR data/16h RT-qPCR/shutoff_mean_values_16h.csv") %>%
  mutate(FoldChange = mean_FC)

# Merge the replicate data
merged_FC_data <- shutoff_FC_16h %>%
  select(row_id, gene, FoldChange, log2FoldChange, condition, fill_color) %>%
  bind_rows(
    rescue_FC %>%
      select(row_id, gene, FoldChange, log2FoldChange, condition, fill_color)
  ) %>% mutate(row_id = factor(row_id, levels = c("shutoff_16h_1", "shutoff_16h_2", "shutoff_16h_3", 
                                                  "FL_rescue_1", "FL_rescue_2", "FL_rescue_3")),
               condition = factor(condition, levels = c("shutoff_16h","FL_rescue")))%>% 
  as.data.frame()

merged_FC_data$fill_color <- factor(merged_FC_data$fill_color, 
                                    levels = c("salmon","gray80",  "skyblue"))

# Merge the mean values and the statistical info for condition vs. wild-type
merged_mean_values <- shutoff_mean_values_16h %>%
  select(gene, condition, FoldChange, log2FoldChange, padj, significant,
         upper, lower, log2Upper, log2Lower, fill_color) %>%
  bind_rows(
    rescue_mean_values %>%
      select(gene, condition, FoldChange, log2FoldChange, padj, significant,
             upper, lower, log2Upper, log2Lower, fill_color)
  ) %>% mutate(condition = factor(condition, levels = c("shutoff_16h","FL_rescue"))) %>% as.data.frame()

###########################################
# Calculate significant differences between the two groups (16h vs rescue)
###########################################

t_test_results_16h_rescue <- merged_FC_data %>%
  filter(condition %in% c("FL_rescue", "shutoff_16h")) %>%
  group_by(gene) %>%
  summarise(
    t_test = list(t.test(log2FoldChange ~ condition, data = cur_data()))
  ) %>%
  mutate(tidy_results = map(t_test, broom::tidy)) %>%
  unnest(tidy_results) %>%
  select(gene, estimate1, estimate2, conf.low, conf.high, p.value) %>%
  rename(log2FC_16h = estimate1,
         log2FC_rescue = estimate2) %>%
  mutate(group1 = "FL_rescue", 
         group2 = "shutoff_16h",
         comparison = "rescue vs 16h",
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

###########################################
# manually set the y.position if necessary
# t_test_results_16h_rescue$y.position <-  c(3.5, 3.75, 3,3.27,3.2,-1.3,-2.92)
# or 
t_test_results_16h_rescue <- t_test_results_16h_rescue %>%
  mutate(y.position = y.position + 0.4)
###########################################

###########################################
# Then Time to PLOT! 

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

  scale_shape_manual(values = c("shutoff_16h" = 21, "FL_rescue" = 24),
                     labels = c("16 hours", 
                                "Full-length Nup211 rescue")) +
  scale_fill_manual(values = c("salmon" = "salmon",
                               "gray80" = "gray80",
                               "skyblue" = "skyblue"),
                    labels = c("Upregulated",
                               "NS",
                               "Downregulated")) +
  
  scale_y_continuous(limits = c(-4.6, 4.0), breaks = seq(-4, 4, 1)) +
  
  labs(title = "RT-qPCR Results (Rescue vs. 16h shutoff)",
       subtitle = bquote(bolditalic("nup211-so") ~ bold("+ pARS1 vs. wild-type (+T)")),
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
    data = t_test_results_16h_rescue,
    aes(
      x = gene,
      xmin = as.numeric(factor(gene)) - 0.15,  # left position (16h)
      xmax = as.numeric(factor(gene)) + 0.15,  # right position (rescue)
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

###############################################################
# Plotting data for now
###############################################################

# ###########################################
# 16h vs rescue
# ###########################################

rescue_16h_FC <- merged_FC_data

merged_mean_values_3 <- merged_mean_values %>%
  mutate(gene = factor(gene, levels = gene_levels)) %>%
  as.data.frame()

merged_mean_values_3 <- merged_mean_values_3 %>%
  mutate(gene = factor(gene, levels = gene_levels)) %>%
  as.data.frame()

rescue_16h_FC <- rescue_16h_FC %>%
  mutate(gene = factor(gene, levels = gene_levels)) %>%
  as.data.frame()

plot_rescue_16h <- ggplot(rescue_16h_FC, aes(x = gene, y = log2FoldChange, fill = condition)) + #fill_color omitted to make plot more legible
  geom_hline(yintercept = 0, linewidth = 1.0, linetype = "solid", color = "black") +
  
  geom_errorbar(data = merged_mean_values_3,
                aes(x = gene, ymin = log2Lower, ymax = log2Upper, group = condition),
                position = position_dodge(width = 0.85),
                width = 0.3, size = 1, inherit.aes = FALSE) +
  
  geom_boxplot(aes(group = interaction(gene, condition)),
               position = position_dodge(width = 0.85),
               coef = 0, outlier.shape = NA,
               size = 0.8, color = "black", alpha = 1) +
  
  geom_jitter(aes(shape = condition),
              position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.85, jitter.height = 0.01),
              size = 2.7, stroke = 1.2, fill = "white", color = "black") +
  
  # geom_text(data = merged_mean_values,
  #           aes(x = gene, y = log2Upper + 0.6, label = round(log2FoldChange, 2), group = condition),
  #           position = position_dodge(0.6),
  #           size = 5, fontface = "bold", inherit.aes = FALSE) +
  
  geom_text(data = merged_mean_values_3,
            aes(x = gene, y = log2Upper + 0.2, label = paste(significant), group = condition),
            position = position_dodge(0.85),
            size = 8, fontface = "bold", inherit.aes = FALSE) +

  scale_shape_manual(values = c("shutoff_16h" = 21, "FL_rescue" = 24),
                     labels = c("<b><i>nup211-so </i></b><b>(16h)</b>",
                                "<b>Nup211<sub>full-length </sub></b><b>rescue</b>")) +
  scale_fill_manual(values = c("shutoff_16h" = "salmon",
                               "FL_rescue" = "turquoise"),
                    labels = c("<b><i>nup211-so </i></b><b>(16h)</b>",
                               "<b>Nup211<sub>full-length </sub></b><b>rescue</b>")) +
  
  scale_y_continuous(limits = c(-4.7, 4.2), breaks = seq(-4, 4, 1)) +
  
  labs(title = "RT-qPCR Results",
       subtitle = bquote(bolditalic("nup211-so") ~ bold("vs. wild-type (+T)")),
       x = NULL,
       y = bquote(bold(log[2] ~ "Fold Change")),
       shape = "Condition: ",
       fill = "Condition: ") +
  
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
    legend.text = ggtext::element_markdown(size = 14, face = "bold"),
    legend.position = "right",
    legend.key.width = unit(1, "cm"),
    legend.spacing.x = unit(0.5, "cm")
    )

plot_rescue_16h

t_test_results_16h_rescue <- t_test_results_16h_rescue %>%
  mutate(gene = factor(gene, levels = gene_levels)) %>%
  as.data.frame()

t_test_results_16h_rescue <-  t_test_results_16h_rescue %>%
  mutate(y.position_2 = y.position - 0.4)

# CANDIDATE final "16h vs rescue" figure (A, "poster"): plot_rescue_16h + brackets from t_test_results_16h_rescue; legend bottom.
plot_rescue_poster <- plot_rescue_16h + 
  geom_signif(
    data = t_test_results_16h_rescue,
    aes(
      x = gene,
      xmin = as.numeric(factor(gene)) - 0.215,  # left position (16h)
      xmax = as.numeric(factor(gene)) + 0.215,  # right position (rescue)
      annotations = significant,
      fontface = "bold",
      y_position = y.position_2 + 0.1
    ),
    tip_length = 0.01,
    size = 1,
    textsize = 8,
    vjust = 0.2,
    inherit.aes = F,
    manual = T
  ) 

##############################################
# Supplemental table building

replicate_data_table <- merged_FC_data %>%
  select(c(1:4)) %>%
  pivot_wider(names_from = row_id, values_from = c(FoldChange, log2FoldChange)) %>%
  `colnames<-`(c("gene",
                 "16h_rep1_FC", "16h_rep2_FC", "16h_rep3_FC",
                 "rescue_rep1_FC", "rescue_rep2_FC", "rescue_rep3_FC",
                 "16h_rep1_log2FC", "16h_rep2_log2FC", "16h_rep3_log2FC",
                 "rescue_rep1_log2FC", "rescue_rep2_log2FC", "rescue_rep3_log2FC"))

mean_data_table <- merged_mean_values %>%
  select(-fill_color) %>%
  pivot_wider(names_from = condition,
              values_from = c(FoldChange, log2FoldChange, padj,
                              significant, upper, lower, log2Upper, log2Lower)) %>%
  `colnames<-`(c("gene",
                 "mean_FC_16h", "mean_FC_rescue",
                 "mean_log2FC_16h", "mean_log2FC_rescue",
                 "padj_16h", "padj_rescue",
                 "annotation_16h", "annotation_rescue",
                 "upperCI_16h", "upperCI_rescue",
                 "lowerCI_16h", "lowerCI_rescue",
                 "log2upperCI_16h", "log2upperCI_rescue",
                 "log2lowerCI_16h", "log2lowerCI_rescue")) %>%
  mutate(annotation_16h    = ifelse(is.na(annotation_16h)    | annotation_16h    == "", "ns", annotation_16h),
         annotation_rescue = ifelse(is.na(annotation_rescue) | annotation_rescue == "", "ns", annotation_rescue))

supplementary_data <- replicate_data_table %>% left_join(mean_data_table, by = "gene")

supplemental_table <- supplementary_data %>%
  rename(padj_16h_v_control       = padj_16h,
         padj_rescue_v_control    = padj_rescue,
         annotation_16h_v_control = annotation_16h,
         annotation_rescue_v_control = annotation_rescue)

supplemental_table <- supplemental_table %>%
  merge(t_test_results_16h_rescue %>%
          select(c(gene, padj, significant)), by = "gene") %>%
  rename(`rescue_vs_16h_padj`       = padj,
         `rescue_vs_16h_annotation` = significant)

rescue_vs_16h_labels <- supplemental_table %>%
  select(gene, rescue_vs_16h_padj, rescue_vs_16h_annotation) %>%
  mutate(gene = factor(gene, levels = gene_levels)) %>%
  left_join(
    merged_mean_values_3  %>%         # was rescue_16h_FC
      group_by(gene) %>%
      summarise(max_y = max(log2Upper, na.rm = TRUE)),
    by = "gene"
  ) %>% 
  mutate(
    y_bracket = max_y + 1.6,
    clean_comp_label = case_when(
      rescue_vs_16h_annotation == "ns" ~ "",
      is.na(rescue_vs_16h_annotation)  ~ "",
      TRUE                             ~ rescue_vs_16h_annotation
    )
  )

#write_csv(x = supplemental_table, file = "~/Desktop/qPCR data/supplemental_files/supplemental_16h_rescue.csv")

## For supplemental thesis stuff: ----------------------------------------------------------------------------------

###########################################
# Data Import
###########################################
nup211_16h <- read_excel("~/Desktop/qPCR data/additional candidate testing/CT_values_20260628-1.xlsx")
nup211_rescue <- read_excel("~/Desktop/qPCR data/additional candidate testing/CT_values_20260628-2.xlsx")

nup211_16h <- nup211_16h %>% 
  mutate(condition = ifelse(grepl("control", sample), "control", "shutoff_16h"),
         bio_rep = as.numeric(factor(gsub(".*_(\\d)$", "\\1", sample))))
nup211_rescue <- nup211_rescue %>% 
  mutate(condition = ifelse(grepl("control", sample), "control", "FL_rescue"),
         bio_rep = as.numeric(factor(gsub(".*_(\\d)$", "\\1", sample))))
####################################################
# Compute average Ct for each biological replicate
####################################################
nup211_16h_rep_avg <- nup211_16h %>%
  group_by(gene, condition, bio_rep) %>%
  summarise(mean_Ct = mean(Ct), .groups = "drop") %>%
  arrange(gene, condition, bio_rep)

nup211_rescue_rep_avg <- nup211_rescue %>%
  group_by(gene, condition, bio_rep) %>%
  summarise(mean_Ct = mean(Ct), .groups = "drop") %>%
  arrange(gene, condition, bio_rep)
####################################################
# Format data for getting individual replicate data 
####################################################
nup211_16h_Ct <- nup211_16h_rep_avg %>%
  mutate(row_id = paste(condition, bio_rep, sep = "_")) %>%  # Create new row ID
  select(row_id, gene, mean_Ct) %>%  # Keep relevant columns
  pivot_wider(names_from = gene, values_from = mean_Ct) 

nup211_rescue_Ct <- nup211_rescue_rep_avg %>%
  mutate(row_id = paste(condition, bio_rep, sep = "_")) %>%  # Create new row ID
  select(row_id, gene, mean_Ct) %>%  # Keep relevant columns
  pivot_wider(names_from = gene, values_from = mean_Ct)

#################################################################
# **** Second Approach: Using manual/prism calculations   ***** #
#################################################################

deltaCt_16h <- nup211_16h_Ct %>%
  pivot_longer(cols = -row_id, names_to = "gene", values_to = "Ct") %>%
  group_by(row_id) %>%
  mutate(act1_Ct = Ct[gene == "act1"]) %>%  # Store Act1 Ct for each sample
  #    filter(!(gene %in% c("act1","adh1","urg1","gcd1"))) %>%  # Remove Act1 from final output
  mutate(deltaCt = Ct - act1_Ct) %>% 
  ungroup()

deltaCt_rescue <- nup211_rescue_Ct %>%
  pivot_longer(cols = -row_id, names_to = "gene", values_to = "Ct") %>%
  group_by(row_id) %>%
  mutate(act1_Ct = Ct[gene == "act1"]) %>%  # Store Act1 Ct for each sample
  #    filter(!(gene %in% c("act1","adh1","urg1","gcd1"))) %>%  # Remove Act1 from final output
  mutate(deltaCt = Ct - act1_Ct) %>% 
  ungroup()

# Calculate the mean of the deltaCts for only the control condition
###################################################################    
mean_deltaCt_16h_controls <- deltaCt_16h %>%
  filter(grepl("control", row_id)) %>%
  group_by(gene) %>%
  summarise(mean_deltaCt_control = mean(deltaCt, na.rm = TRUE))

mean_deltaCt_rescue_controls <- deltaCt_rescue %>%
  filter(grepl("control", row_id)) %>%
  group_by(gene) %>%
  summarise(mean_deltaCt_control = mean(deltaCt, na.rm = TRUE))

# Merge the dataframes and calculate DeltaDeltaCt for each gene/condition
#########################################################################   
df_ddCt_16h <- deltaCt_16h %>%
  left_join(mean_deltaCt_16h_controls, by = "gene") %>%
  mutate(DeltaDeltaCt = deltaCt - mean_deltaCt_control)

df_ddCt_resue <- deltaCt_rescue %>%
  left_join(mean_deltaCt_rescue_controls, by = "gene") %>%
  mutate(DeltaDeltaCt = deltaCt - mean_deltaCt_control)

# Then calculate FoldChange (2^-DeltaDeltaCt) for each gene/condition
# and calculate log2FoldChange
######################################################################### 

df_FC_16h <- df_ddCt_16h %>% mutate(FoldChange = 2^(-DeltaDeltaCt),
                                log2FoldChange = log2(FoldChange),
                                condition = case_when(
                                  grepl("control", row_id, ignore.case = TRUE) ~ "control",
                                  grepl("shutoff_16h", row_id, ignore.case = TRUE) ~ "shutoff_16h",
                                  TRUE ~ "other" # In case there are unexpected values
                                )) %>% 
  filter(gene != "act1")
df_FC_rescue <- df_ddCt_resue %>% mutate(FoldChange = 2^(-DeltaDeltaCt),
                                log2FoldChange = log2(FoldChange),
                                condition = case_when(
                                  grepl("control", row_id, ignore.case = TRUE) ~ "control",
                                  grepl("FL_rescue", row_id, ignore.case = TRUE) ~ "FL_rescue",
                                  TRUE ~ "other" # In case there are unexpected values
                                )) %>% 
  filter(gene != "act1",
         condition != "control")

df_FC_nup211 <- bind_rows(df_FC_16h, df_FC_rescue)  

df_FC_rescue_2 <- df_ddCt_resue %>% mutate(FoldChange = 2^(-DeltaDeltaCt),
                                         log2FoldChange = log2(FoldChange),
                                         condition = case_when(
                                           grepl("control", row_id, ignore.case = TRUE) ~ "control",
                                           grepl("FL_rescue", row_id, ignore.case = TRUE) ~ "FL_rescue",
                                           TRUE ~ "other" # In case there are unexpected values
                                         )) %>% 
  filter(gene != "act1")

################################################################################ 
# Now, that we have the expression values for each replicate, we can calculate
# the  average/mean values for each transcript in both conditions 
#
#          Both a sanity check and useful for plotting results
################################################################################ 

mean_values_nup211 <- df_FC_nup211 %>%
  group_by(gene, condition) %>%
  summarise(log2FoldChange = mean(log2FoldChange, na.rm = TRUE), 
            FoldChange = 2^(log2FoldChange), 
            .groups = "drop") 

# After, I can calculate pvalues and padj's independently of the pcr R program 
# just to show either data can be used for plotting significance! 
#
#  --> Again, a sanity check and potentially useful for plotting results

t_test_results_16h <- df_FC_16h %>%
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
         group2 = "FL_rescue",
         condition = "FL_rescue",
         padj = p.adjust(p.value, method = "BH"),
         significant = case_when(
           padj < 0.001 ~ "***",
           padj < 0.01  ~ "**",
           padj < 0.05  ~ "*",
           TRUE          ~ "ns")
  ) %>% group_by(gene) %>%
  mutate(y.position = max(df_FC_16h$log2FoldChange[df_FC_16h$gene == gene]) + 0.5) %>%
  ungroup()

################################################################################ 
t_test_results_rescue <- df_FC_rescue_2 %>%
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
         group2 = "FL_rescue",
         condition = "FL_rescue",
         padj = p.adjust(p.value, method = "BH"),
         significant = case_when(
           padj < 0.001 ~ "***",
           padj < 0.01  ~ "**",
           padj < 0.05  ~ "*",
           TRUE          ~ "ns")
  ) %>% group_by(gene) %>%
  mutate(y.position = max(df_FC_rescue_2$log2FoldChange[df_FC_rescue_2$gene == gene]) + 0.5) %>%
  ungroup()

###########################################
# Calculate significant differences between the two groups (16h vs rescue)
###########################################

t_test_results_16h_rescue <- df_FC_nup211 %>%
  filter(condition %in% c("FL_rescue", "shutoff_16h")) %>%
  group_by(gene) %>%
  summarise(
    t_test = list(t.test(log2FoldChange ~ condition, data = cur_data()))
  ) %>%
  mutate(tidy_results = map(t_test, broom::tidy)) %>%
  unnest(tidy_results) %>%
  select(gene, estimate1, estimate2, conf.low, conf.high, p.value) %>%
  rename(log2FC_16h = estimate1,
         log2FC_rescue = estimate2) %>%
  mutate(group1 = "FL_rescue", 
         group2 = "shutoff_16h",
         comparison = "rescue vs 16h",
         padj = p.adjust(p.value, method = "BH"),
         significant = case_when(
           padj < 0.001 ~ "***",
           padj < 0.01  ~ "**",
           padj < 0.05  ~ "*",
           TRUE         ~ "ns")
  ) %>%
  group_by(gene) %>%
  mutate(y.position = max(df_FC_nup211$log2FoldChange[df_FC_nup211$gene == gene]) + 1.1) %>%
  ungroup() 

################################################################################ 
# 
################################################################################ 
condition_levels <- c("shutoff_16h", "FL_rescue")

FC_nup211 <- df_FC_nup211 %>% filter(condition %in% c("FL_rescue", "shutoff_16h")) %>% mutate(condition = factor(condition, levels = condition_levels))
mean_FC_nup211 <- mean_values_nup211 %>% filter(condition %in% c("FL_rescue", "shutoff_16h")) %>% mutate(condition = factor(condition, levels = condition_levels))

FC_nup211 
mean_FC_nup211
t_test_results_16h
t_test_results_rescue
t_test_results_16h_rescue

# Per-condition summary for the supplemental panel: 95% CI error bars, star labels
# and label positions. (Replaces the main-analysis rescue_16h_FC, which holds a
# different gene panel and uses log2Upper/log2Lower.)
nup211_summary <- FC_nup211 %>%
  group_by(gene, condition) %>%
  summarise(log2mean = mean(log2FoldChange, na.rm = TRUE),
            n        = sum(!is.na(log2FoldChange)),
            sem      = sd(log2FoldChange, na.rm = TRUE) / sqrt(n),
            .groups  = "drop") %>%
  mutate(half      = ifelse(n > 1, qt(0.975, df = n - 1) * sem, NA_real_),
         log2upper = log2mean + half,
         log2lower = log2mean - half)

nup211_summary <- nup211_summary %>%
  left_join(
    bind_rows(
      t_test_results_16h     %>% mutate(condition = "shutoff_16h"),
      t_test_results_rescue  %>% mutate(condition = "FL_rescue")
    ) %>% select(gene, condition, significant),
    by = c("gene", "condition")
  ) %>%
  mutate(condition   = factor(condition, levels = condition_levels),
         clean_label = ifelse(is.na(significant) | significant == "ns", "", significant),
         mean_label  = sprintf("%.2f", log2mean),
         mean_y_pos  = ifelse(log2mean >= 0, log2upper + 0.15, log2lower - 0.45),
         star_y_pos  = ifelse(log2mean >= 0, log2upper + 0.95, log2lower - 1.05))

plot_nup211_rescue_16h <- ggplot(FC_nup211, aes(x = gene, y = log2FoldChange, fill = condition)) +
  geom_hline(yintercept = 0, linewidth = 1.0, color = "black") +
  
  # 2. Boxplot & Jitter (Using the 30-row raw data)
  geom_boxplot(aes(group = interaction(gene, condition)),
               position = position_dodge(width = 0.85),
               coef = 0, outlier.shape = NA, size = 0.8) +
  geom_jitter(aes(shape = condition),
              position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.85),
              size = 2.7, stroke = 1.2, fill = "white", color = "black") +
  
  # 3. Error Bars (Switching to the 10-row summary data)
  geom_errorbar(data = nup211_summary,
                aes(x = gene, ymin = log2lower, ymax = log2upper, group = condition),
                position = position_dodge(width = 0.85),
                width = 0.3, linewidth = 1, inherit.aes = FALSE) +
  
  # 4. Condition-specific Stars (Using the 10-row summary data)
  # LAYER A: The Numeric Mean (Black, bold, 2 sig figs)
  # geom_text(data = nup211_summary,
  #           aes(x = gene, y = mean_y_pos, label = mean_label, group = condition),
  #           color = "black",
  #           position = position_dodge(0.85),
  #           size = 4.5, fontface = "bold.italic", vjust = 0, inherit.aes = FALSE) +
  
  # LAYER B: The Significance Stars (Black, bold)
  geom_text(data = nup211_summary,
            aes(x = gene, y = star_y_pos-0.7, label = clean_label, group = condition),
            color = "black",
            position = position_dodge(0.85),
            size = 8, fontface = "bold", vjust = 0, inherit.aes = FALSE) +
  
  # 5. Scales & Styling (From your template)
  scale_shape_manual(values = c("16h" = 21, "rescue" = 24),
                     labels = c("<b><i>nup211-so </i></b><b>(16h)</b>",
                                "<b>Ectopic Nup211 rescue</b>")) +
  scale_fill_manual(values = c("16h" = "salmon", "rescue" = "turquoise"),
                    labels = c("<b><i>nup211-so </i></b><b>(16h)</b>",
                               "<b>Ectopic Nup211 rescue</b>")) +
  scale_y_continuous(limits = c(-1.8, 5), breaks = seq(-1, 5, 1)) +
  labs(title = "RT-qPCR Results",
       subtitle = bquote(bolditalic("nup211-so") ~ bold("vs. wild-type (+T)")),
       x = NULL,
       y = bquote(bold(log[2] ~ "Fold Change")),
       shape = "Condition: ",
       fill = "Condition: ") +
  #theme_minimal() +
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
    legend.text = ggtext::element_markdown(size = 14, face = "bold"),
    legend.position = "right",
    legend.key.width = unit(1, "cm"),
    legend.spacing.x = unit(0.5, "cm")
  )

plot_rescue_16h

# Create a separate dataframe for the Rescue vs 16h comparison
# Bracket labels for the supplemental (7-gene) panel.
# Sourced from this section's own t_test_results_16h_rescue (df_FC_nup211-based),
# NOT supplemental_table, which holds the main-analysis gene panel.
rescue_vs_16h_labels <- t_test_results_16h_rescue %>%
  select(gene, padj, significant) %>%
  rename(padj_rescue_vs_16h       = padj,
         annotation_rescue_vs_16h = significant) %>%
  mutate(gene = factor(gene, levels = c("chp1", "mit1", "swi6", "for3", "tea2", "mus7", "not1"))) %>%
  # Bracket height from the highest datapoint per gene (FC_nup211 has no CI columns)
  left_join(
    FC_nup211 %>%
      group_by(gene) %>%
      summarise(max_y = max(log2FoldChange, na.rm = TRUE)),
    by = "gene"
  ) %>%
  mutate(
    y_bracket = max_y + 1.6,
    clean_comp_label = case_when(
      annotation_rescue_vs_16h == "ns" ~ "",
      is.na(annotation_rescue_vs_16h)  ~ "",
      TRUE                             ~ annotation_rescue_vs_16h
    )
  )

rescue_vs_16h_labels_clean <- rescue_vs_16h_labels %>%
  filter(clean_comp_label != "")  # remove non-significant brackets

# 6. Add the Significance Brackets (Using the 5-row label data)
# CANDIDATE final "16h vs rescue" figure (B, "WIP"): same base + brackets from rescue_vs_16h_labels_clean.
plot_rescue_WIP <- plot_rescue_16h + 
  geom_signif(
    data = rescue_vs_16h_labels_clean,
    aes(
      x = gene,
      xmin = as.numeric(gene) - 0.22,  # Aligns with 16h bar center
      xmax = as.numeric(gene) + 0.22,  # Aligns with rescue bar center
      annotations = clean_comp_label, 
      y_position = y_bracket-0.7,           # Using your calculated column
      fontface = "bold" 
    ),
    manual = TRUE, 
    tip_length = 0.02, 
    size = 1, 
    textsize = 8, 
    inherit.aes = FALSE
  )

plot_rescue_WIP

# ---- 0. Settings --------------------------------------------------------
err_type <- "ci"      # "ci" = 95% CI, "sem", or "sd"
fnt      <- "Arial"
dodge_w  <- 0.60
box_w    <- 0.50      # each box ends up box_w / 2 wide after dodging

cond_levels <- c("shutoff_16h", "FL_rescue")
cond_cols   <- c(shutoff_16h = "salmon", FL_rescue = "turquoise")
cond_shapes <- c(shutoff_16h = 21, FL_rescue = 24)
cond_labs   <- c(shutoff_16h = "<b><i>nup211-so</i></b> <b>(16h)</b>",
                 FL_rescue   = "<b>Ectopic Nup211</b>")

FC_nup211 <- FC_nup211 %>%
  mutate(condition = factor(as.character(condition), levels = cond_levels))

# ---- 1. Summary, error bars, per-condition stars ------------------------
summ <- FC_nup211 %>%
  group_by(gene, condition) %>%
  summarise(log2mean = mean(log2FoldChange), sd = sd(log2FoldChange),
            n = n(), .groups = "drop") %>%
  mutate(se   = sd / sqrt(n),
         half = case_when(err_type == "ci"  ~ qt(0.975, n - 1) * se,
                          err_type == "sem" ~ se,
                          TRUE              ~ sd),
         log2lower = log2mean - half,
         log2upper = log2mean + half)

cond_stats <- bind_rows(
  t_test_results_16h    %>% transmute(gene, condition = "shutoff_16h", significant),
  t_test_results_rescue %>% transmute(gene, condition = "FL_rescue",   significant)
) %>% mutate(condition = factor(condition, levels = cond_levels))

summ <- summ %>%
  left_join(cond_stats, by = c("gene", "condition")) %>%
  mutate(clean_label = coalesce(significant, ""),
         star_y      = if_else(log2mean >= 0, log2upper + 0.25, log2lower - 1.05))

# ---- 2. Shared theme ----------------------------------------------------
theme_qpcr <- theme_grey(base_family = fnt) +
  theme(
    plot.title        = element_text(size = 18, face = "bold"),
    plot.subtitle     = element_text(size = 15, face = "bold"),
    panel.border      = element_blank(),
    axis.line         = element_line(linewidth = 1.0, colour = "black"),
    axis.ticks        = element_line(linewidth = 1.0, colour = "black"),
    axis.ticks.length = unit(2, "mm"),
    axis.text.x       = element_text(size = 18, face = "bold", colour = "black",
                                     angle = 0, hjust = 0.5), # or angle = 30
    axis.text.y       = element_text(size = 15, face = "bold", colour = "black"),
    axis.title.y      = element_text(size = 18, face = "bold", colour = "black"),
    legend.title      = element_text(size = 14, face = "bold"),
    legend.text       = element_markdown(size = 13),
    legend.position   = "right",
    legend.key.width  = unit(1, "cm")
  )

sub_lab <- bquote(bolditalic("nup211-so") ~ bold("(16h) vs. wild-type (+T)"))

# ---- 3. PLOT A: 16h only ------------------------------------------------
d16 <- filter(FC_nup211, condition == "shutoff_16h")
s16 <- filter(summ,      condition == "shutoff_16h")

plot_16h <- ggplot(d16, aes(x = gene, y = log2FoldChange, fill = condition)) +
  geom_hline(yintercept = 0, linewidth = 1.0, colour = "black") +
  geom_boxplot(width = box_w / 2, coef = 0, outlier.shape = NA, linewidth = 0.8) +
  geom_errorbar(data = s16, aes(x = gene, ymin = log2lower, ymax = log2upper),
                width = 0.12, linewidth = 1, inherit.aes = FALSE) +
  geom_jitter(aes(shape = condition), position = position_jitter(width = 0.06, height = 0),
              size = 2.7, stroke = 1.2, fill = "white", colour = "black") +
  geom_text(data = s16, aes(x = gene, y = star_y, label = clean_label),
            size = 8, fontface = "bold", vjust = 0, family = fnt, inherit.aes = FALSE) +
  scale_fill_manual(values = cond_cols, labels = cond_labs, drop = TRUE) +
  scale_shape_manual(values = cond_shapes, labels = cond_labs, drop = TRUE) +
  scale_x_discrete(expand = expansion(add = 0.55)) +
  scale_y_continuous(breaks = seq(-6, 2, 1)) +
  coord_cartesian(ylim = c(-6, 1.5)) +
  labs(title = "RT-qPCR", subtitle = sub_lab, x = NULL,
       y = bquote(bold(log[2] ~ "Fold Change")),
       fill = "Condition: ", shape = "Condition: ") +
  theme_qpcr

# ---- 4. PLOT B: 16h + FL rescue ----------------------------------------
# NOTE: re-derives gene order from the data for this block's bracket x-positions (overrides the display constant above)
gene_levels <- sort(unique(FC_nup211$gene))
off         <- dodge_w / 4          # dodged box centres, +/- 0.15
y_brk       <- max(summ$log2upper, na.rm = TRUE) + 1.3

brackets <- t_test_results_16h_rescue %>%
  transmute(gene, annot = significant,
            xmin = match(gene, gene_levels) - off,
            xmax = match(gene, gene_levels) + off,
            y_bracket = y_brk) %>%
  filter(!is.na(annot), annot != "")

# CANDIDATE final "16h vs rescue" figure (C, theme_qpcr rebuild): legend right; the commented ggsave() at the bottom targets this one.
plot_16h_rescue <- ggplot(FC_nup211, aes(x = gene, y = log2FoldChange, fill = condition)) +
  geom_hline(yintercept = 0, linewidth = 1.0, colour = "black") +
  geom_boxplot(aes(group = interaction(gene, condition)), width = box_w,
               position = position_dodge(width = dodge_w),
               coef = 0, outlier.shape = NA, linewidth = 0.8) +
  geom_errorbar(data = summ, aes(x = gene, ymin = log2lower, ymax = log2upper,
                                 group = condition),
                position = position_dodge(width = dodge_w),
                width = 0.12, linewidth = 1, inherit.aes = FALSE) +
  geom_jitter(aes(shape = condition),
              position = position_jitterdodge(jitter.width = 0.12, dodge.width = dodge_w),
              size = 2.7, stroke = 1.2, fill = "white", colour = "black") +
  geom_text(data = summ, aes(x = gene, y = star_y, label = clean_label, group = condition),
            position = position_dodge(dodge_w),
            size = 8, fontface = "bold", vjust = 0, family = fnt, inherit.aes = FALSE) +
  geom_signif(data = brackets,
              aes(xmin = xmin, xmax = xmax, annotations = annot, y_position = y_bracket),
              manual = TRUE, inherit.aes = FALSE,
              tip_length = 0.02, size = 1, textsize = 8, fontface = "bold", family = fnt) +
  scale_fill_manual(values = cond_cols, labels = cond_labs) +
  scale_shape_manual(values = cond_shapes, labels = cond_labs) +
  scale_x_discrete(expand = expansion(add = 0.55)) +
  scale_y_continuous(breaks = seq(-6, 10, 2)) +
  coord_cartesian(ylim = c(-5, y_brk + 0.6)) +
  labs(title = "RT-qPCR", subtitle = sub_lab, x = NULL,
       y = bquote(bold(log[2] ~ "Fold Change")),
       fill = "Condition: ", shape = "Condition: ") +
  theme_qpcr

plot_16h
plot_16h_rescue

#ggsave("nup211_16h.pdf",        plot_16h,        width = 5.5, height = 5.5, device = cairo_pdf)
#ggsave("nup211_16h_rescue.pdf", plot_16h_rescue, width = 6.5, height = 6.0, device = cairo_pdf)

