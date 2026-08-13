#################################################
# RT-qPCR Analysis  (655 CRISPR) -
# Purpose: Creating Plots / Data Visualization
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
top_transcripts_655_1 <- read_excel("~/Desktop/qPCR data/655 RT-qPCR/CT_values_TopTranscripts_655_CRISPR_1-2.xlsx")
top_transcripts_655_2 <- read_excel("~/Desktop/qPCR data/655 RT-qPCR/CT_values_TopTranscripts_655_CRISPR_2-2.xlsx")

top_transcripts_655_1 <- top_transcripts_655_1 %>% 
  mutate(condition = ifelse(grepl("control", sample), "control", "Nup211_655"),
         bio_rep = as.numeric(factor(gsub(".*_(\\d)$", "\\1", sample))))
top_transcripts_655_2 <- top_transcripts_655_2 %>% 
  mutate(condition = ifelse(grepl("control", sample), "control", "Nup211_655"),
         bio_rep = as.numeric(factor(gsub(".*_(\\d)$", "\\1", sample))))

####################################################
# Compute average Ct for each biological replicate
####################################################
tech_rep_avg_655_1 <- top_transcripts_655_1 %>%
  group_by(gene, condition, bio_rep) %>%
  summarise(mean_Ct = mean(Ct), .groups = "drop") %>%
  arrange(gene, condition, bio_rep)

tech_rep_avg_655_2 <- top_transcripts_655_2 %>%
  group_by(gene, condition, bio_rep) %>%
  summarise(mean_Ct = mean(Ct), .groups = "drop") %>%
  arrange(gene, condition, bio_rep)
####################################################
# Format data for getting individual replicate data 
####################################################
topTranscripts_655_Ct_1 <- tech_rep_avg_655_1 %>%
  mutate(row_id = paste(condition, bio_rep, sep = "_")) %>%  # Create new row ID
  select(row_id, gene, mean_Ct) %>%  # Keep relevant columns
  pivot_wider(names_from = gene, values_from = mean_Ct) 

topTranscripts_655_Ct_2 <- tech_rep_avg_655_2 %>%
  mutate(row_id = paste(condition, bio_rep, sep = "_")) %>%  # Create new row ID
  select(row_id, gene, mean_Ct) %>%  # Keep relevant columns
  pivot_wider(names_from = gene, values_from = mean_Ct)

# write_csv(topTranscripts_655_Ct_1, file = "~/Desktop/qPCR data/655 RT-qPCR/topTranscripts_Nup211655_formatted_Ct_1.csv")
# write_csv(topTranscripts_655_Ct_2, file = "~/Desktop/qPCR data/655 RT-qPCR/topTranscripts_Nup211655_formatted_Ct_2.csv")

####################################################
# Format data for use with the PCR R program 
####################################################
topTranscripts_655_wide_1 <- tech_rep_avg_655_1 %>%
  mutate(row_id = paste(condition, bio_rep, sep = "_")) %>%  # Create new row ID
  select(row_id, gene, mean_Ct) %>%  # Keep relevant columns
  pivot_wider(names_from = gene, values_from = mean_Ct) %>%
  select(-row_id) # blank this out in case you want to see if the Ct values are in the right order

topTranscripts_655_wide_2 <- tech_rep_avg_655_2 %>%
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
group <- rep(c('Nup211_655','control'), each = 3)

# Get expression information
topTranscripts_655_FC_1 <-  pcr_ddct(topTranscripts_655_wide_1,
                                     group_var = group,
                                     reference_gene = 'act1', #or actin
                                     reference_group = 'control',
                                     mode = "separate_tube",
                                     plot = F)
topTranscripts_655_FC_2 <-  pcr_ddct(topTranscripts_655_wide_2,
                                     group_var = group,
                                     reference_gene = 'act1', #or actin
                                     reference_group = 'control',
                                     mode = "separate_tube",
                                     plot = F)

# Get significance information 
topTranscripts_655_stats_1 <- pcr_ttest(topTranscripts_655_wide_1,
                                        group_var = group,
                                        reference_gene = 'act1',
                                        reference_group = 'control',
                                        tidy = T)
topTranscripts_655_stats_2 <- pcr_ttest(topTranscripts_655_wide_2,
                                        group_var = group,
                                        reference_gene = 'act1',
                                        reference_group = 'control',
                                        tidy = T)
# Correct for multiple hypothesis testing
# (BH adjustment):
topTranscripts_655_stats_1 <- topTranscripts_655_stats_1 %>% 
  mutate(padj = p.adjust(p_value, method = "BH"),
         significant = case_when(
           padj < 0.001 ~ "***",
           padj < 0.01  ~ "**",
           padj < 0.05  ~ "*",
           TRUE          ~ "ns"
         )
  )
topTranscripts_655_stats_2 <- topTranscripts_655_stats_2 %>% 
  mutate(padj = p.adjust(p_value, method = "BH"),
         significant = case_when(
           padj < 0.001 ~ "***",
           padj < 0.01  ~ "**",
           padj < 0.05  ~ "*",
           TRUE          ~ "ns"
         )
  )
# Merge significance and expression information 
topTranscripts_655_data_1 <- topTranscripts_655_FC_1 %>% 
  left_join(topTranscripts_655_stats_1 %>% 
              select(gene, estimate, p_value, padj, significant), 
            by = "gene",
            keep = F) %>%
  rename(FoldChange = relative_expression,
         pvalue = p_value,
         condition = group) %>%
  mutate(log2FoldChange = log2(FoldChange),
         log2Upper = log2(upper),  # Log2 transform upper value directly
         log2Lower = log2(lower))  # Log2 transform lower value directly

topTranscripts_655_data_2 <- topTranscripts_655_FC_2 %>% 
  left_join(topTranscripts_655_stats_2 %>% 
              select(gene, estimate, p_value, padj, significant), 
            by = "gene",
            keep = F) %>%
  rename(FoldChange = relative_expression,
         pvalue = p_value,
         condition = group) %>%
  mutate(log2FoldChange = log2(FoldChange),
         log2Upper = log2(upper),  # Log2 transform upper value directly
         log2Lower = log2(lower))  # Log2 transform lower value directly

topTranscripts_655_data <- bind_rows(topTranscripts_655_data_1, topTranscripts_655_data_2)  

# Get confidence interval information

CI_bars_655 <- topTranscripts_655_data %>% select(gene,condition,upper,lower,log2Upper,log2Lower) %>% 
  filter(condition == "Nup211_655")

# write_csv(CI_bars_655, file = "~/Desktop/qPCR data/655 RT-qPCR/topTranscripts_Nup211655_CI_bars.csv")

#***#
#################################################################
# **** Second Approach: Using manual/prism calculations   ***** #
#################################################################

# Calculate delta Cts for each transcript per condition 
###############################################################
deltaCt_df_1 <- topTranscripts_655_Ct_1 %>%
  pivot_longer(cols = -row_id, names_to = "gene", values_to = "Ct") %>%
  group_by(row_id) %>%
  mutate(act1_Ct = Ct[gene == "act1"]) %>%  # Store Act1 Ct for each sample
  #    filter(!(gene %in% c("act1","adh1","urg1","gcd1"))) %>%  # Remove Act1 from final output
  mutate(deltaCt = Ct - act1_Ct) %>% 
  ungroup()

deltaCt_df_2 <- topTranscripts_655_Ct_2 %>%
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
                                  grepl("Nup211_655", row_id, ignore.case = TRUE) ~ "Nup211_655",
                                  TRUE ~ "other" # In case there are unexpected values
                                )) %>% 
  filter(gene != "act1")
df_FC_2 <- df_ddCt_2 %>% mutate(FoldChange = 2^(-DeltaDeltaCt),
                                log2FoldChange = log2(FoldChange),
                                condition = case_when(
                                  grepl("control", row_id, ignore.case = TRUE) ~ "control",
                                  grepl("Nup211_655", row_id, ignore.case = TRUE) ~ "Nup211_655",
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
# # write_csv(df_FC_prism, file = "~/Desktop/qPCR data/655 RT-qPCR/topTranscripts_Nup211655_FC_prism.csv")

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
         group2 = "Nup211_655",
         condition = "Nup211_655",
         padj = p.adjust(p.value, method = "BH"),
         significant = case_when(
           padj < 0.001 ~ "***",
           padj < 0.01  ~ "**",
           abs(log2FoldChange) >= 1 & padj < 0.05  ~ "*",
           TRUE          ~ "ns")
  ) %>% group_by(gene) %>%
  mutate(y.position = max(df_FC$log2FoldChange[df_FC$gene == gene]) + 0.5) %>%
  ungroup()

################################################################################ 
# 
################################################################################ 

FC_655 <- df_FC %>%
  filter(condition == "Nup211_655") %>%
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

mean_values_655 <- mean_values %>%
  filter(condition == "Nup211_655") %>%
  left_join(CI_bars_655, by = c("gene", "condition")) %>%
  select(c(1:2,4:6,3,7:8)) %>%
  left_join(
    t_test_signif %>%
      select(gene, padj, significant, y.position, fill_color) %>%
      distinct(),
    by ="gene"
  ) %>%
  mutate(gene = factor(gene, levels = gene_levels)) %>%
  as.data.frame()

###########################################
# Make the plot!
###########################################

CRISPR655_plot <- ggplot(FC_655, aes(x = gene, y = log2FoldChange, fill = condition)) +
  geom_hline(yintercept = 0, linewidth = 1.0, linetype = "solid", color = "black", alpha = 0.9) +
  geom_errorbar(data = mean_values_655, aes(x = gene, ymin = log2Lower, ymax = log2Upper),
                width = 0.2, size = 1.2) +
  geom_boxplot(alpha = 1, width=0.5, position = "dodge2", coef = 0, outlier.shape = NA, size = 0.8, color = "black") +
  geom_jitter(width = 0.23, size = 3, alpha = 0.9, shape = 21, fill = "white", color = "black", stroke = 1.2) +
  # geom_text(data = rescue_mean_values, aes(x = gene, log2Upper + 0.7,
  #                                           label = round(log2FoldChange, 2)), position = position_dodge(0.9),
  #           size = 6, fontface = "bold") +
  geom_text(data = mean_values_655, aes(x = gene, y = log2Upper + 0.2, #max(df_FC$log2FoldChange) - 0.1,
                                           label = paste(significant)), position = position_dodge(1),
            size = 12, fontface = "bold") +
  theme_gray() +
  labs(title = "RT-qPCR Results",
       subtitle = "<i><b>nup211<sub>1-655</sub></b></i> <b>vs. wild-type</b>",
       legend = "",
       y = bquote(bold(log[2]~ 'Fold Change')),
       x = "",
       fill = "") +
  scale_y_continuous(limits = c(-4,2),
                     breaks=seq(-4,2, 1)) +
  scale_fill_manual(values = c("Nup211_655" = "sienna2"),
                    labels = c("<i>nup211<sub>1-655</sub></i>")) + 
  theme(plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_markdown(size = 16),
        axis.line = element_line(size = 1.0), # Control axis thickness
        axis.ticks = element_line(size = 1.2), # Control tick mark thickness
        axis.ticks.length = unit(3, "mm"), # Control tick mark length
        axis.text.x = element_text(size = 18, face= "bold",colour = 'black',
                                   angle = 30, hjust = 1),    # Axis labels
        axis.text.y = element_text(size = 16, face= "bold",colour = 'black'),
        axis.title.y = element_text(size = 18, face= "bold",colour = 'black'),   # Axis titles
        legend.title = element_text(size = 14, face = "bold"),   # Adjust legend title font
        legend.text = element_markdown(size = 14, face = "bold"), # Adjust legend labels font
        legend.position = "bottom") # Adjust legend labels font

print(CRISPR655_plot)

# write_csv(rescue_FC, file = "Desktop/qPCR data/Rescue RT-qPCR/rescue_FC_value_points.csv")
# write_csv(rescue_mean_values, file = "Desktop/qPCR data/Rescue RT-qPCR/rescue_mean_values_significance.csv" )

