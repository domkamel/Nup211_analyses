#################################################
# RT-qPCR Analysis  (863 CRISPR) -
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
top_transcripts_863_1 <- read_excel("~/Desktop/qPCR data/863 RT-qPCR/CT_values_TopTranscripts_863_CRISPR_1-3.xlsx")
top_transcripts_863_2 <- read_excel("~/Desktop/qPCR data/863 RT-qPCR/CT_values_TopTranscripts_863_CRISPR_2-3.xlsx")

top_transcripts_863_1 <- top_transcripts_863_1 %>% 
  mutate(condition = ifelse(grepl("control", sample), "control", "Nup211_863"),
         bio_rep = as.numeric(factor(gsub(".*_(\\d)$", "\\1", sample))))
top_transcripts_863_2 <- top_transcripts_863_2 %>% 
  mutate(condition = ifelse(grepl("control", sample), "control", "Nup211_863"),
         bio_rep = as.numeric(factor(gsub(".*_(\\d)$", "\\1", sample))))

####################################################
# Compute average Ct for each biological replicate
####################################################
tech_rep_avg_863_1 <- top_transcripts_863_1 %>%
  group_by(gene, condition, bio_rep) %>%
  summarise(mean_Ct = mean(Ct), .groups = "drop") %>%
  arrange(gene, condition, bio_rep)

tech_rep_avg_863_2 <- top_transcripts_863_2 %>%
  group_by(gene, condition, bio_rep) %>%
  summarise(mean_Ct = mean(Ct), .groups = "drop") %>%
  arrange(gene, condition, bio_rep)
####################################################
# Format data for getting individual replicate data 
####################################################
topTranscripts_863_Ct_1 <- tech_rep_avg_863_1 %>%
  mutate(row_id = paste(condition, bio_rep, sep = "_")) %>%  # Create new row ID
  select(row_id, gene, mean_Ct) %>%  # Keep relevant columns
  pivot_wider(names_from = gene, values_from = mean_Ct) 

topTranscripts_863_Ct_2 <- tech_rep_avg_863_2 %>%
  mutate(row_id = paste(condition, bio_rep, sep = "_")) %>%  # Create new row ID
  select(row_id, gene, mean_Ct) %>%  # Keep relevant columns
  pivot_wider(names_from = gene, values_from = mean_Ct)

# write_csv(topTranscripts_863_Ct_1, file = "~/Desktop/qPCR data/655 RT-qPCR/topTranscripts_Nup211863_formatted_Ct_1.csv")
# write_csv(topTranscripts_863_Ct_2, file = "~/Desktop/qPCR data/655 RT-qPCR/topTranscripts_Nup211863_formatted_Ct_2.csv")

####################################################
# Format data for use with the PCR R program 
####################################################
topTranscripts_863_wide_1 <- tech_rep_avg_863_1 %>%
  mutate(row_id = paste(condition, bio_rep, sep = "_")) %>%  # Create new row ID
  select(row_id, gene, mean_Ct) %>%  # Keep relevant columns
  pivot_wider(names_from = gene, values_from = mean_Ct) %>%
  select(-row_id) # blank this out in case you want to see if the Ct values are in the right order

topTranscripts_863_wide_2 <- tech_rep_avg_863_2 %>%
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
group <- rep(c('Nup211_863','control'), each = 3)

# Get expression information
topTranscripts_863_FC_1 <-  pcr_ddct(topTranscripts_863_wide_1,
                                        group_var = group,
                                        reference_gene = 'act1',
                                        reference_group = 'control',
                                        mode = "separate_tube",
                                        plot = F)
topTranscripts_863_FC_2 <-  pcr_ddct(topTranscripts_863_wide_2,
                                        group_var = group,
                                        reference_gene = 'act1',
                                        reference_group = 'control',
                                        mode = "separate_tube",
                                        plot = F)

# Get significance information 
topTranscripts_863_stats_1 <- pcr_ttest(topTranscripts_863_wide_1,
                                           group_var = group,
                                           reference_gene = 'act1',
                                           reference_group = 'control',
                                           tidy = T)
topTranscripts_863_stats_2 <- pcr_ttest(topTranscripts_863_wide_2,
                                           group_var = group,
                                           reference_gene = 'act1',
                                           reference_group = 'control',
                                           tidy = T)
# Correct for multiple hypothesis testing
# (BH adjustment):
topTranscripts_863_stats_1 <- topTranscripts_863_stats_1 %>% 
  mutate(padj = p.adjust(p_value, method = "BH"),
         significant = case_when(
           padj < 0.001 ~ "***",
           padj < 0.01  ~ "**",
           padj < 0.05  ~ "*",
           TRUE          ~ "ns"
         )
  )
topTranscripts_863_stats_2 <- topTranscripts_863_stats_2 %>% 
  mutate(padj = p.adjust(p_value, method = "BH"),
         significant = case_when(
           padj < 0.001 ~ "***",
           padj < 0.01  ~ "**",
           padj < 0.05  ~ "*",
           TRUE          ~ "ns"
         )
  )
# Merge significance and expression information 
topTranscripts_863_data_1 <- topTranscripts_863_FC_1 %>% 
  left_join(topTranscripts_863_stats_1 %>% 
              select(gene, estimate, p_value, padj, significant), 
            by = "gene",
            keep = F) %>%
  rename(FoldChange = relative_expression,
         pvalue = p_value,
         condition = group) %>%
  mutate(log2FoldChange = log2(FoldChange),
         log2Upper = log2(upper),  # Log2 transform upper value directly
         log2Lower = log2(lower))  # Log2 transform lower value directly

topTranscripts_863_data_2 <- topTranscripts_863_FC_2 %>% 
  left_join(topTranscripts_863_stats_2 %>% 
              select(gene, estimate, p_value, padj, significant), 
            by = "gene",
            keep = F) %>%
  rename(FoldChange = relative_expression,
         pvalue = p_value,
         condition = group) %>%
  mutate(log2FoldChange = log2(FoldChange),
         log2Upper = log2(upper),  # Log2 transform upper value directly
         log2Lower = log2(lower))  # Log2 transform lower value directly

topTranscripts_863_data <- bind_rows(topTranscripts_863_data_1, topTranscripts_863_data_2)  

# Get confidence interval information

CI_bars_863 <- topTranscripts_863_data %>% select(gene,condition,upper,lower,log2Upper,log2Lower) %>% 
  filter(condition == "Nup211_863")

#***#
#################################################################
# **** Second Approach: Using manual/prism calculations   ***** #
#################################################################

# Calculate delta Cts for each transcript per condition 
###############################################################
deltaCt_df_1 <- topTranscripts_863_Ct_1 %>%
  pivot_longer(cols = -row_id, names_to = "gene", values_to = "Ct") %>%
  group_by(row_id) %>%
  mutate(act1_Ct = Ct[gene == "act1"]) %>%  # Store Act1 Ct for each sample
  #    filter(!(gene %in% c("act1","adh1","urg1","gcd1"))) %>%  # Remove Act1 from final output
  mutate(deltaCt = Ct - act1_Ct) %>% 
  ungroup()

deltaCt_df_2 <- topTranscripts_863_Ct_2 %>%
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
                                  grepl("Nup211_863", row_id, ignore.case = TRUE) ~ "Nup211_863",
                                  TRUE ~ "other" # In case there are unexpected values
                                )) %>% 
  filter(gene != "act1")
df_FC_2 <- df_ddCt_2 %>% mutate(FoldChange = 2^(-DeltaDeltaCt),
                                log2FoldChange = log2(FoldChange),
                                condition = case_when(
                                  grepl("control", row_id, ignore.case = TRUE) ~ "control",
                                  grepl("Nup211_863", row_id, ignore.case = TRUE) ~ "Nup211_863",
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
         group2 = "Nup211_863",
         condition = "Nup211_863",
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

FC_863 <- df_FC %>%
  filter(condition == "Nup211_863") %>%
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

mean_values_863 <- mean_values %>%
  filter(condition == "Nup211_863") %>%
  left_join(CI_bars_863,
            by = c("gene", "condition")
  ) %>% select(c(1:2,4:6,3,7:8))  %>%
  left_join(
    t_test_signif %>%
      select(gene, padj, significant, y.position, fill_color) %>%
      distinct(),  # Ensure uniqueness before joining
    by ="gene"
  ) %>%
  mutate(gene = factor(gene, levels = gene_levels)) %>%
  as.data.frame()

###########################################
# Make the plot!
###########################################

CRISPR863_plot <- ggplot(FC_863, aes(x = gene, y = log2FoldChange, fill = condition)) +
  geom_hline(yintercept = 0, linewidth = 1.0, linetype = "solid", color = "black", alpha = 0.9) +
  geom_errorbar(data = mean_values_863, aes(x = gene, ymin = log2Lower-0.1, ymax = log2Upper+0.1),
                width = 0.2, size = 1.2) +
  geom_boxplot(alpha = 1, width=0.5, position = "dodge2", coef = 0, outlier.shape = NA, size = 0.8, color = "black") +
  geom_jitter(width = 0.23, size = 3, alpha = 0.9, shape = 21, fill = "white", color = "black", stroke = 1.2) +
  # geom_text(data = rescue_mean_values, aes(x = gene, log2Upper + 0.7,
  #                                           label = round(log2FoldChange, 2)), position = position_dodge(0.9),
  #           size = 6, fontface = "bold") +
  geom_text(data = mean_values_863, aes(x = gene, y = log2Upper + 0.2, #max(df_FC$log2FoldChange) - 0.1,
                                           label = paste(significant)), position = position_dodge(1),
            size = 12, fontface = "bold") +
  theme_gray() +
  labs(title = "RT-qPCR Results",
       subtitle = "<i><b>nup211<sub>1-863</sub></b></i> <b>vs. wild-type</b>",
       legend = "",
       y = bquote(bold(log[2]~ 'Fold Change')),
       x = "",
       fill = "Strain") +
  scale_y_continuous(limits = c(-4,3),
                     breaks=seq(-4,3, 1)) +
  scale_fill_manual(values = c("Nup211_863" = "thistle"), #or violet
                    labels = c("<i>nup211<sub>1-863</sub></i>")) + 
  theme(plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_markdown(size = 16),
        axis.line = element_line(linewidth = 1.0), # Control axis thickness
        axis.ticks = element_line(linewidth = 1.2), # Control tick mark thickness
        axis.ticks.length = unit(3, "mm"), # Control tick mark length
        axis.text.x = element_text(size = 18, face= "bold",colour = 'black',
                                   angle = 30, hjust = 1),    # Axis labels
        axis.text.y = element_text(size = 16, face= "bold",colour = 'black'),
        axis.title.y = element_text(size = 18, face= "bold",colour = 'black'),   # Axis titles
        legend.title = element_text(size = 14, face = "bold"),   # Adjust legend title font
        legend.text = element_markdown(size = 14, face = "bold"), # Adjust legend labels font
        legend.position = "bottom") # Adjust legend labels font

print(CRISPR863_plot)

# write_csv(rescue_FC, file = "Desktop/qPCR data/Rescue RT-qPCR/rescue_FC_value_points.csv")
# write_csv(rescue_mean_values, file = "Desktop/qPCR data/Rescue RT-qPCR/rescue_mean_values_significance.csv" )

###########################################
# Now merge the 655 data!
###########################################

shutoff_FC_16h <- read_csv("~/Desktop/qPCR data/16h RT-qPCR/shutoff_FC_16h.csv") 
shutoff_mean_values_16h <- read_csv("~/Desktop/qPCR data/16h RT-qPCR/shutoff_mean_values_16h.csv") %>%
  mutate(FoldChange = mean_FC)

# DEPENDENCY: FC_655 / mean_values_655 / CI_bars_655 are created by RTqPCR-655.R.
# Run that script first in the same session.
if (!exists("FC_655") || !exists("mean_values_655")) {
  stop("Run RTqPCR-655.R first -- FC_655 and mean_values_655 are needed for the 863 vs 655 comparison.")
}

# Merge the replicate data
merged_FC_data <- FC_863 %>%
  select(row_id, gene, FoldChange, log2FoldChange, condition) %>%
  bind_rows(
    FC_655 %>%
      select(row_id, gene, FoldChange, log2FoldChange, condition)
  ) %>% mutate(row_id = factor(row_id, levels = c("Nup211_863_1", "Nup211_863_2", "Nup211_863_3", 
                                                  "Nup211_655_1", "Nup211_655_2", "Nup211_655_3")),
               condition = factor(condition, levels = c("Nup211_863","Nup211_655")))%>% 
  as.data.frame()

# Merge the mean values and the statistical info for condition vs. wild-type
merged_mean_values <- mean_values_863 %>%
  select(gene, condition, FoldChange, log2FoldChange, padj, significant,
         upper, lower, log2Upper, log2Lower) %>%
  bind_rows(
    mean_values_655 %>%
      select(gene, condition, FoldChange, log2FoldChange, padj, significant,
             upper, lower, log2Upper, log2Lower)
  ) %>% mutate(condition = factor(condition, levels = c("Nup211_863","Nup211_655"))) %>% as.data.frame()

###########################################
# Calculate significant differences between the two groups (863 vs 655)
###########################################

t_test_863vs655 <- merged_FC_data %>%
  filter(condition %in% c("Nup211_863","Nup211_655")) %>%
  group_by(gene) %>%
  summarise(
    t_test = list(t.test(log2FoldChange ~ condition, data = cur_data()))
  ) %>%
  mutate(tidy_results = map(t_test, broom::tidy)) %>%
  unnest(tidy_results) %>%
  select(gene, estimate1, estimate2, conf.low, conf.high, p.value) %>%
  rename(log2FC_16h = estimate1,
         log2FC_rescue = estimate2) %>%
  mutate(group1 = "Nup211_863", 
         group2 = "Nup211_655",
         comparison = "Nup211_863 vs Nup211_655",
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
t_test_863vs655 <- t_test_863vs655 %>%
  mutate(y.position = y.position + 0.4,
         significant = ifelse(significant == "ns", "", significant))
###########################################

###########################################
# Then Time to PLOT! 

p <- ggplot(merged_FC_data, aes(x = gene, y = log2FoldChange, fill = condition )) +
  geom_hline(yintercept = 0, linewidth = 1.0, linetype = "solid", color = "black") +
  
  geom_errorbar(data = merged_mean_values,
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
  
  geom_text(data = merged_mean_values,
            aes(x = gene, y = log2Upper + 0.2, label = paste(significant), group = condition),
            position = position_dodge(0.85),
            size = 8, fontface = "bold", inherit.aes = FALSE) +

  scale_shape_manual(values = c("Nup211_863" = 22, "Nup211_655" = 25),
                     labels = c("<i>nup211<sub>1-863</sub></i>",
                                "<i>nup211<sub>1-655</sub></i>")) +
  scale_fill_manual(values = c("Nup211_863" = "thistle",
                               "Nup211_655" = "sienna2"),
                    labels = c("<i>nup211<sub>1-863</sub></i>",
                               "<i>nup211<sub>1-655</sub></i>")) +
  
  scale_y_continuous(limits = c(-4.7, 4.2), breaks = seq(-4, 4, 1)) +
  
  labs(title = "RT-qPCR Results",
       subtitle =  "<b><i>nup211</i><sub><i>1-863</i></sub></b> <b>and</b> <b><i>nup211</i><sub><i>1-655</i></sub></b> <b>vs. wild-type</b>",
       x = NULL,
       y = bquote(bold(log[2] ~ "Fold Change")),
       shape = "Strain: ",
       fill = "Strain: ") +
  
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_markdown(size = 16),
    axis.line = element_line(size = 1.0),
    axis.ticks = element_line(size = 1.2),
    axis.ticks.length = unit(3, "mm"),
    axis.text.x = element_text(size = 18, face = "bold", colour = 'black', angle = 30, hjust = 1),
    axis.text.y = element_text(size = 16, face = "bold", colour = 'black'),
    axis.title.y = element_text(size = 18, face = "bold", colour = 'black'),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_markdown(size = 14, face = "bold"),
    legend.position = "bottom",
    legend.key.width = unit(1, "cm"),
    legend.spacing.x = unit(0.5, "cm")
  )

print (p)

t_test_863vs655_2 <- t_test_863vs655 %>% 
  filter(!(gene %in% c("tea2","urg2")))

plot_863_655 <- p + 
  geom_signif(
    data = t_test_863vs655,
    aes(
      x = gene,
      xmin = as.numeric(factor(gene)) - 0.215,  # left position (16h)
      xmax = as.numeric(factor(gene)) + 0.215,  # right position (48h)
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

replicate_data_table_CRISPR <-  merged_FC_data %>% select(c(1:4)) %>% pivot_wider(names_from = row_id, values_from = c(FoldChange,log2FoldChange)) %>%
  `colnames<-`(c("gene", "Nup211_863_rep1_FC","Nup211_863_rep2_FC","Nup211_863_rep3_FC","Nup211_655_rep1_FC","Nup211_655_rep2_FC","Nup211_655_rep3_FC",
                 "Nup211_863rep1_log2FC","Nup211_863_rep2_log2FC","Nup211_863_rep3_log2FC", "Nup211_655_rep1_log2FC","Nup211_655_rep2_log2FC","Nup211_655_rep3_log2FC"))

mean_data_table_CRISPR <-  merged_mean_values %>% pivot_wider(names_from = condition, values_from = c(FoldChange,log2FoldChange, padj, 
                                                                                                                significant, upper,lower, log2Upper, log2Lower))  %>%
  `colnames<-`(c("gene", "mean_FC_Nup211_863","mean_FC_Nup211_655","mean_log2FC_Nup211_863","mean_log2FC_Nup211_655","padj_Nup211_863","padj_Nup211_655",
                 "annotation_Nup211_863","annotation_Nup211_655","upperCI_Nup211_863","upperCI_Nup211_655","lowerCI_Nup211_863","lowerCI_Nup211_655","log2upperCI_Nup211_863",
                 "log2upperCI_Nup211_655","log2lowerCI_Nup211_863","log2lowerCI_Nup211_655")) %>%
  mutate(annotation_Nup211_863 = ifelse(annotation_Nup211_863 == "", "ns", annotation_Nup211_863),
         annotation_Nup211_655 = ifelse(annotation_Nup211_655 == "", "ns", annotation_Nup211_655)) %>% merge(t_test_863vs655 %>% 
                                                                                          select(c(gene,padj,significant)), by = 'gene') %>%
  rename(`Nup211_863_vs_Nup211_655_padj`= padj,
         `1Nup211_863_vs_Nup211_655_annotation`= significant)

supplementary_data_CRISPR <- replicate_data_table_CRISPR %>% left_join(mean_data_table_CRISPR, by = 'gene')

supplementary_table_CRISPR <-  supplementary_data_CRISPR %>% select(c(1,8:10,11:13,16,26,28,17,27,29,18,20,19,21,30,31)) %>%
  rename(padj_Nup211_863_v_control = 14,
         annotation_Nup211_863_v_control = 15,
         padj_Nup211_655_v_control = 16,
         annotation_Nup211_655_v_control = 17,
         padj_Nup211_863_vs_Nup211_655 = 18,
         annotation_Nup211_863_vs_Nup211_655 = 19) 

# write_csv(x = supplementary_table_CRISPR, file = "~/Desktop/qPCR data/supplemental_files/supplemental_CRISPR_RTqPCR.csv")

