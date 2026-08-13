#Intro##################################################
# RNA-seq Analysis Pipeline (16h) 
# Testing additional leading edge genes by RT-qPCR
#
# Updated: 3 conditions (control, shutoff_16h, FL_rescue)
#
# SCALABLE FUNCTIONS INCLUDED to hopefully make things easier
###################################################

#Setup##########################################
# Setup and Library Loading
###########################################
libraries <- c(
  "readr", "readxl", "dplyr", "tidyr", "tidyverse",
  "ggplot2","ggpubr", "ggExtra", "plotly", "gplots", "Glimma",
  "pheatmap", "RColorBrewer", "viridis", "PoiClaClu", 
  "kableExtra", "broom", "rstatix", "ggsignif","ggtext"
)

for (lib in libraries) {
  if (!require(lib, character.only = TRUE)) {
    BiocManager::install(lib)
    library(lib, character.only = TRUE)
  }
}

# Gene display orders used in this script
# Main dissertation panel (leading-edge RT-qPCR):
gene_levels <- c("chp1", "mit1", "swi6", "for3", "tea2", "mus7", "not1")
# Extended panel used for the CRISPR strains (adds bub1, mad1, alm1):
gene_levels_CRISPR <- c("chp1", "mit1", "swi6", "bub1", "mad1", "for3", "tea2", "mus7", "not1", "alm1")

#Import#################################################
# Data Import - read all excel sheets into a named list
#################################################
file_paths <- list(
  leading_transcripts_1 = "~/Desktop/qPCR data/additional candidate testing/CT_values_20260307_1-2.xlsx",
  leading_transcripts_2 = "~/Desktop/qPCR data/additional candidate testing/CT_values_20260307_2-2.xlsx",
  leading_transcripts_3 = "~/Desktop/qPCR data/additional candidate testing/CT_values_20260307_3-2.xlsx",
  leading_transcripts_4 = "~/Desktop/qPCR data/additional candidate testing/CT_values_20260307_4.xlsx"
)

# Read all files in one step - adding a new file is just one new line above
raw_data <- map(file_paths, read_excel)

#Cleaning##########################################
# Condition mapping - applied to all dataframes at once
###########################################
raw_data <- map(raw_data, ~ .x %>%
                  mutate(condition = case_when(
                    grepl("control",      sample) ~ "control",
                    grepl("shutoff_16h",  sample) ~ "shutoff_16h",
                    grepl("pARS1_rescue", sample) ~ "FL_rescue",
                    TRUE                          ~ NA_character_
                  ),
                  bio_rep = as.numeric(factor(gsub(".*_(\\d)$", "\\1", sample))))
)

### Manual approach turned into a pipeline / function ➡️ calculate ΔCt,meanΔCt, ΔΔCt, replicate FCs, and mean FCs ###

#Pipeline##########################################
# Pipeline function - do this ONCE, use for all dataframes
###########################################
run_qpcr_pipeline <- function(df, comparisons, condition_levels) {
  
  # Technical replicate averaging
  tech_rep_avg <- df %>%
    group_by(sample, condition, bio_rep, gene) %>%
    summarise(mean_Ct = mean(Ct, na.rm = TRUE), .groups = "drop")
  
  # Format for downstream analysis
  ct_wide <- tech_rep_avg %>%
    mutate(row_id = paste(condition, bio_rep, sep = "_")) %>%
    select(row_id, gene, mean_Ct) %>%
    pivot_wider(names_from = gene, values_from = mean_Ct)
  
  # Delta Ct
  deltaCt_df <- ct_wide %>%
    pivot_longer(cols = -row_id, names_to = "gene", values_to = "Ct") %>%
    group_by(row_id) %>%
    mutate(act1_Ct = Ct[gene == "act1"],
           deltaCt = Ct - act1_Ct) %>%
    ungroup()
  
  # Mean delta Ct for controls
  mean_deltaCt_controls <- deltaCt_df %>%
    filter(grepl("control", row_id)) %>%
    group_by(gene) %>%
    summarise(mean_deltaCt_control = mean(deltaCt, na.rm = TRUE))
  
  # DeltaDeltaCt + fold changes
  df_FC <- deltaCt_df %>%
    left_join(mean_deltaCt_controls, by = "gene") %>%
    mutate(
      DeltaDeltaCt   = deltaCt - mean_deltaCt_control,
      FoldChange     = 2^(-DeltaDeltaCt),
      log2FoldChange = log2(FoldChange),
      condition = case_when(
        grepl("control",     row_id, ignore.case = TRUE) ~ "control",
        grepl("shutoff_16h", row_id, ignore.case = TRUE) ~ "shutoff_16h",
        grepl("FL_rescue",   row_id, ignore.case = TRUE) ~ "FL_rescue",
        TRUE ~ NA_character_
      ),
      condition = factor(condition, levels = condition_levels)
    )
  
  # Pairwise Welch's t-tests
  t_test_results <- df_FC %>%
    filter(!is.na(log2FoldChange)) %>%
    group_by(gene) %>%
    reframe(
      map_dfr(comparisons, function(comp) {
        grp1 <- log2FoldChange[condition == comp[1]]
        grp2 <- log2FoldChange[condition == comp[2]]
        
        if (length(grp1) < 2 | length(grp2) < 2) {
          return(tibble(group1 = comp[1], group2 = comp[2],
                        estimate = NA, conf.low = NA, conf.high = NA,
                        p.value = NA))
        }
        
        tidy(t.test(grp1, grp2)) %>%
          mutate(group1 = comp[1], group2 = comp[2])
      })
    ) %>%
    select(gene, group1, group2, estimate, conf.low, conf.high, p.value)
  
  # BH correction per comparison
  # NOTE: Alternative is to remove group_by(group1, group2) and correct
  # across all tests simultaneously for more conservative FDR control
  t_test_results <- t_test_results %>%
    group_by(group1, group2) %>%
    mutate(
      padj = p.adjust(p.value, method = "BH"),
      significant = case_when(
        padj < 0.001 ~ "***",
        padj < 0.01  ~ "**",
        padj < 0.05  ~ "*",
        TRUE          ~ "ns"
      )
    ) %>%
    ungroup()
  
  # y positions for significance brackets
  t_test_results <- t_test_results %>%
    group_by(gene) %>%
    mutate(
      max_log2FC = max(df_FC$log2FoldChange[df_FC$gene == cur_group()$gene],
                       na.rm = TRUE),
      y.position = max_log2FC + 0.5 * row_number()
    ) %>%
    ungroup()
  
  # Mean summary
  mean_values <- df_FC %>%
    group_by(gene, condition) %>%
    summarise(
      mean_log2FC = mean(log2FoldChange, na.rm = TRUE),
      mean_FC     = 2^(mean_log2FC),
      .groups = "drop"
    )
  
  # Return everything you need downstream as a named list
  return(list(
    df_FC          = df_FC,
    t_test_results = t_test_results,
    mean_values    = mean_values
  ))
}

#Run Pipeline##########################################
# Run pipeline on all dataframes at once
###########################################
condition_levels <- c("control", "shutoff_16h", "FL_rescue")

comparisons <- list(
  c("control",     "shutoff_16h"),
  c("control",     "FL_rescue"),
  c("shutoff_16h", "FL_rescue")
)

# This runs the entire pipeline on every dataframe - adding leading_transcripts_4
# later just means adding one line to file_paths above, nothing else changes
results <- map(raw_data, ~ run_qpcr_pipeline(.x, comparisons, condition_levels))

names(results$leading_transcripts_1$df_FC)
names(results$leading_transcripts_1$mean_values)

#Pool Results###############################################################################
# Pool results across all datasets
################################################################################

# 1. Replicate-level data (for individual points on boxplots)
################################################################################
pooled_replicates <- map_dfr(results, ~ .x$df_FC, .id = "dataset") %>%
  filter(gene != "act1") %>%
  select(dataset, gene, condition, row_id, log2FoldChange)

clean_pooled_replicates <- pooled_replicates %>%
  filter(condition != "control") %>%
  dplyr::rename(log2FC = log2FoldChange) 

# 2. Mean log2FC + 95% CIs per gene/condition
#    log2Upper/log2Lower are log2 transformed from conf.high/conf.low
################################################################################
# 2. Mean log2FC + 95% CIs per gene/condition
################################################################################
pooled_means <- map_dfr(results, ~ .x$df_FC, .id = "dataset") %>%
  filter(gene != "act1") %>%
  group_by(dataset, gene, condition) %>%
  summarise(
    mean_log2FC = mean(log2FoldChange, na.rm = TRUE),
    mean_FC     = 2^(mean_log2FC),
    n           = n(),
    se          = sd(log2FoldChange, na.rm = TRUE) / sqrt(n),
    ci          = qt(0.975, df = n - 1) * se,
    log2upper   = mean_log2FC + ci,
    log2lower   = mean_log2FC - ci,
    .groups = "drop"
  ) %>%
  left_join(
    map_dfr(results, ~ .x$t_test_results, .id = "dataset") %>%
      filter(group1 == "control") %>%
      dplyr::rename(condition = group2) %>%
      select(dataset, gene, condition, padj, significant),
    by = c("dataset", "gene", "condition")
  ) 

clean_pooled_means <- pooled_means %>%
  filter(condition != "control") %>%
  mutate(condition = factor(condition, levels = c("control", "shutoff_16h", "FL_rescue")))

# 3a. Control comparisons (control vs shutoff_16h, control vs FL_rescue)
################################################################################
pooled_sig_vs_control <- map_dfr(results, ~ .x$t_test_results, .id = "dataset") %>%
  filter(group1 == "control") %>%
  mutate(significant = ifelse(significant == "ns", "", significant)) %>%
  left_join(
    pooled_means %>% select(dataset, gene, condition, log2upper),
    by = c("dataset", "gene", "group2" = "condition")
  ) %>%
  group_by(dataset, gene) %>%
  mutate(y.position = max(log2upper, na.rm = TRUE) + 0.5 * row_number()) %>%
  ungroup() %>%
  select(dataset, gene, group1, group2, padj, significant, y.position) 

pooled_sig_vs_control <- pooled_sig_vs_control %>%
  filter(gene != "act1") %>%
  mutate(condition = factor(group2, levels = c("control", "shutoff_16h", "FL_rescue")))

# 3b. Rescue vs shutoff comparison
################################################################################
pooled_sig_rescue_vs_shutoff <- map_dfr(results, ~ .x$t_test_results, .id = "dataset") %>%
  filter(group1 == "shutoff_16h" & group2 == "FL_rescue") %>%
  mutate(significant = ifelse(significant == "ns", "", significant)) %>%
  left_join(
    pooled_means %>% select(dataset, gene, condition, log2upper),
    by = c("dataset", "gene", "group2" = "condition")
  ) %>%
  group_by(dataset, gene) %>%
  mutate(y.position = max(log2upper, na.rm = TRUE) + 0.5 * row_number()) %>%
  ungroup() %>%
  select(dataset, gene, group1, group2, padj, significant, y.position)

pooled_sig_rescue_vs_shutoff <- pooled_sig_rescue_vs_shutoff %>%
  filter(gene != "act1")

#Plotting####
# Objects for plotting: 
clean_pooled_replicates
clean_pooled_means  
pooled_sig_vs_control
pooled_sig_rescue_vs_shutoff

plot_leading_rescue <- ggplot(clean_pooled_replicates, aes(x = gene, y = log2FC, fill = condition)) +
  geom_hline(yintercept = 0, linewidth = 1.0, color = "black") +
  
  # 2. Boxplot & Jitter (Using the 30-row raw data)
  geom_boxplot(aes(group = interaction(gene, condition)),
               position = position_dodge(width = 0.85),
               coef = 0, outlier.shape = NA, size = 0.8) +
  geom_jitter(aes(shape = condition),
              position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.85),
              size = 2.7, stroke = 1.2, fill = "white", color = "black") +
  
  # 3. Error Bars (Switching to the 10-row summary data)
  geom_errorbar(data = clean_pooled_means,
                aes(x = gene, ymin = log2lower, ymax = log2upper, group = condition),
                position = position_dodge(width = 0.85),
                width = 0.3, linewidth = 1, inherit.aes = FALSE) +
  
  # 4. Condition-specific Stars (Using the 10-row summary data)
  # LAYER A: The Numeric Mean (Black, bold, 2 sig figs)
  geom_text(data = clean_pooled_means,
            aes(x = gene, y = log2upper + 0.215 , label = sprintf("%.2f", mean_log2FC), group = condition),
            color = "black",
            position = position_dodge(0.85),
            size = 4.5, fontface = "bold.italic", vjust = 0, inherit.aes = FALSE) +
  
  # LAYER B: The Significance Stars (Black, bold)
  geom_text(data = pooled_sig_vs_control,
            aes(x = gene, y = y.position, label = significant, group = condition),
            color = "black",
            position = position_dodge(0.85),
            size = 8, fontface = "bold", vjust = 0, inherit.aes = FALSE) +
  
  # 5. Scales & Styling (From your template)
  scale_shape_manual(values = c("shutoff_16h" = 21, "FL_rescue" = 24),
                     labels = c("<b><i>nup211-so </i></b><b>(16h)</b>",
                                "<b>Ectopic Nup211 rescue</b>")) +
  scale_fill_manual(values = c("shutoff_16h" = "salmon", "FL_rescue" = "turquoise"),
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
    axis.line = element_line(linewidth = 1.0),
    axis.ticks = element_line(linewidth = 1.2),
    axis.ticks.length = unit(3, "mm"),
    axis.text.x = element_text(size = 18, face = "bold", colour = 'black', angle = 30, hjust = 1),
    axis.text.y = element_text(size = 16, face = "bold", colour = 'black'),
    axis.title.y = element_text(size = 18, face = "bold", colour = 'black'),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = ggtext::element_markdown(size = 14, face = "bold"),
    legend.position = "bottom",
    legend.key.width = unit(1, "cm"),
    legend.spacing.x = unit(0.5, "cm")
  )

# 6. Add the Significance Brackets (Using the 5-row label data)

pooled_sig_rescue_vs_shutoff <- pooled_sig_rescue_vs_shutoff %>%
  filter(significant != "") %>%  # remove non-significant brackets
  left_join(
    clean_pooled_means %>%
      group_by(gene) %>%
      summarise(max_log2upper = max(log2upper, na.rm = TRUE)),
    by = "gene"
  ) %>%
  mutate(y.position = max_log2upper + 1.2)  # adjust offset as needed

plot_leading_rescue_2 <- plot_leading_rescue + 
  geom_signif(
    data = pooled_sig_rescue_vs_shutoff,
    aes(
      x          = gene,
      xmin       = as.numeric(factor(gene, levels = c("alm1", "bub1", "mad1", "not1", "swi6"))) - 0.22,
      xmax       = as.numeric(factor(gene, levels = c("alm1", "bub1", "mad1", "not1", "swi6"))) + 0.22,
      annotations = significant,
      fontface   = "bold",
      y_position = y.position
    ),
    tip_length  = 0.01,
    size        = 1,
    textsize    = 8,
    vjust       = 0.2,
    inherit.aes = FALSE,
    manual      = TRUE
  )

plot_leading_rescue_2

#Supplemental tables###############################################################################
# Generate supplemental table (16h + rescue only, no 48h)
################################################################################
supp_table <- pooled_replicates %>%
  filter(condition != "control") %>%
  group_by(gene, condition) %>%
  mutate(rep = row_number(),
         condition = factor(condition, levels = condition_levels)) %>%
  select(-c(row_id, dataset)) %>%
  ungroup() %>%
  pivot_wider(
    names_from  = c(condition, rep),
    values_from = log2FoldChange,
    names_glue  = "{condition}_rep{rep}_log2FC"
  ) %>% 
  select(c(1,5:7,2:4)) %>%
  
  # Join mean + CI
  left_join(
    clean_pooled_means %>%
      select(gene, condition, mean_log2FC, log2upper, log2lower) %>%
      pivot_wider(
        names_from  = condition,
        values_from = c(mean_log2FC, log2upper, log2lower),
        names_glue  = "{.value}_{condition}"
      ),
    by = "gene"
  ) %>%
  
  # Join padj and annotations for 16h vs control and rescue vs control
  left_join(
    pooled_sig_vs_control %>%
      select(gene, group2, padj, significant) %>%
      pivot_wider(
        names_from  = group2,
        values_from = c(padj, significant),
        names_glue  = "{.value}_{group2}_v_control"
      ),
    by = "gene"
  ) %>%
  
  # Join padj and annotation for rescue vs 16h
  left_join(
    pooled_sig_rescue_vs_shutoff %>%
      select(gene, padj, significant) %>%
      dplyr::rename(
        padj_rescue_vs_16h       = padj,
        annotation_rescue_vs_16h = significant
      ),
    by = "gene"
  ) %>%
  
  # Clean up column names to match your existing format
  rename_with(~ gsub("significant", "annotation", .x)) %>%
  rename_with(~ gsub("shutoff_16h", "16h", .x)) %>%
  rename_with(~ gsub("FL_rescue", "rescue", .x)) %>%
  rename_with(~ gsub("log2upper", "log2upperCI", .x)) %>%
  rename_with(~ gsub("log2lower", "log2lowerCI", .x))

# Write to csv
# write.csv(supp_table, "Desktop/qPCR data/additional candidate testing/additional_supplemental_16h_rescue.csv", row.names = FALSE)

# chp1, mit1, swi6, bub1, mad1, for3, tea2, mus7, not1, alm1

# Read in the supplemental table info: 
# NOTE: RTqPCR-rescue.R now writes supplemental_16h_rescue.csv (48h data removed).

rescue_16h_table <- read_csv("~/Desktop/qPCR data/supplemental_files/supplemental_16h_rescue.csv") 

rescue_16h_table <- read_csv("~/Desktop/qPCR data/supplemental_files/supplemental_16_rescue_48h.csv") 
# Columns in A but not in B
setdiff(names(supp_table), names(rescue_16h_table))

# Columns in B but not in A
setdiff(names(rescue_16h_table), names(supp_table))

# Fix mismatched names:

rescue_16h_table <- rescue_16h_table %>%
  dplyr::rename(
    annotation_rescue_vs_16h = rescue_vs_16h_annotation,
    padj_rescue_vs_16h = rescue_vs_16h_padj
  )

setdiff(names(supp_table), names(rescue_16h_table))
setdiff(names(rescue_16h_table), names(supp_table))

supplemental_table_16h_rescue_RTqPCR <- bind_rows(
  rescue_16h_table,
  supp_table
)

supplemental_table_16h_rescue_RTqPCR <- supplemental_table_16h_rescue_RTqPCR %>% 
  select(!contains("48h"))

# Write to csv
# write.csv(supplemental_table_16h_rescue_RTqPCR, "Documents/RNA-seq manuscript/supplemental_table_16h_rescue_48h.csv", row.names = FALSE)

# Rename for ease of use:
working_supp <- supplemental_table_16h_rescue_RTqPCR

# Filter out the 48h data: 
working_supp <- working_supp %>% select(-contains("48h")) 

# Filter out data for pho1 and urg2: 
rows_to_remove <- c("pho1", "urg2","bub1","mad1","alm1")

working_supp <- working_supp %>%
  filter(!(gene %in% rows_to_remove))

# Pivot the replicate data points into long format again: 
rescue_16h_reps <- select(working_supp, c(1:7)) %>%
  # Pivot the replicate columns (rep1, rep2, rep3) to long format
  pivot_longer(
    cols = contains("_rep"), 
    names_to = c("condition", "replicate"), 
    names_pattern = "(.*)_rep(\\d+)_log2FC",
    values_to = "log2FC"
  )  %>%
  mutate(gene = factor(gene, levels = gene_levels))

rescue_16h_FC <- select(working_supp, c(1, 8:17)) %>%
  mutate(gene = factor(gene, levels = gene_levels)) %>%
  dplyr::rename(
    log2FC_16h       = mean_log2FC_16h,
    log2FC_rescue    = mean_log2FC_rescue,
    log2upper_16h    = log2upperCI_16h,
    log2upper_rescue = log2upperCI_rescue,
    log2lower_16h    = log2lowerCI_16h,
    log2lower_rescue = log2lowerCI_rescue,
    padj_16h         = annotation_16h_v_control, 
    padj_rescue      = padj_rescue_v_control,
    annot_16h        = annotation_16h_v_control, 
    annot_rescue     = annotation_rescue_v_control
  ) %>%
  pivot_longer(
    cols = contains("_16h") | contains("_rescue"),
    names_to = c(".value", "condition"),
    names_pattern = "(.*)_(16h|rescue)"
  ) %>% 
  group_by(gene, condition) %>%
  mutate(
    # Use 'log2upper' here so the label sits above the error bar, not just the mean
    y_position = max(log2upper, na.rm = TRUE) + 0.15
  ) %>%
  ungroup()

rescue_16h_FC <- rescue_16h_FC %>%
  mutate(
    # Create a clean label: if annot is 'ns', make it blank, otherwise keep the annot
    clean_label = case_when(
      annot == "ns" ~ "",      # If it says 'ns', make it empty
      is.na(annot)  ~ "",      # Handle potential NAs
      TRUE          ~ annot    # Otherwise, keep the original (e.g., "*", "**")
    )
  )

rescue_16h_FC <- rescue_16h_FC %>%
  mutate(
    # 1. Format the existing mean (log2FC) to 2 significant figures
    mean_label = sprintf("%.2f", log2FC),
    
    # 2. Position for the Number: Just above the error bar
    mean_y_pos = log2upper + 0.25,
    
    # 3. Position for the Star: Above the number
    star_y_pos = log2upper + 0.85
  )

# Create a separate dataframe for the Rescue vs 16h comparison
rescue_vs_16h_labels <- working_supp %>%
  # Select gene and the specific comparison columns
  select(gene, padj_rescue_vs_16h, annotation_rescue_vs_16h) %>%
  # Ensure genes match your specific factor levels
  mutate(gene = factor(gene, levels = gene_levels)) %>%
  # Find the highest value from the ALREADY created long dataframe to set bracket height
  # summarise only gives us the gene name and the max_y heights
  left_join(
    rescue_16h_FC %>% 
      group_by(gene) %>% 
      summarise(max_y = max(log2upper, na.rm = TRUE)),
    by = "gene"
  ) %>%
  mutate(
    # Set bracket position 1 unit above the highest error bar
    y_bracket = max_y + 1.6,
    # Blank out 'ns' labels
    clean_comp_label = case_when(
      annotation_rescue_vs_16h == "ns" ~ "",
      is.na(annotation_rescue_vs_16h) ~ "",
      TRUE ~ annotation_rescue_vs_16h
    )
  )

plot_rescue_16h <- ggplot(rescue_16h_reps, aes(x = gene, y = log2FC, fill = condition)) +
  geom_hline(yintercept = 0, linewidth = 1.0, color = "black") +
  
  # 2. Boxplot & Jitter (Using the 30-row raw data)
  geom_boxplot(aes(group = interaction(gene, condition)),
               position = position_dodge(width = 0.85),
               coef = 0, outlier.shape = NA, size = 0.8) +
  geom_jitter(aes(shape = condition),
              position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.85),
              size = 2.7, stroke = 1.2, fill = "white", color = "black") +
  
  # 3. Error Bars (Switching to the 10-row summary data)
  geom_errorbar(data = rescue_16h_FC,
                aes(x = gene, ymin = log2lower, ymax = log2upper, group = condition),
                position = position_dodge(width = 0.85),
                width = 0.3, linewidth = 1, inherit.aes = FALSE) +
  
  # 4. Condition-specific Stars (Using the 10-row summary data)
  # LAYER A: The Numeric Mean (Black, bold, 2 sig figs)
  # geom_text(data = rescue_16h_FC,
  #           aes(x = gene, y = mean_y_pos, label = mean_label, group = condition),
  #           color = "black",
  #           position = position_dodge(0.85),
  #           size = 4.5, fontface = "bold.italic", vjust = 0, inherit.aes = FALSE) +
  
  # LAYER B: The Significance Stars (Black, bold)
  geom_text(data = rescue_16h_FC,
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

rescue_vs_16h_labels_clean <- rescue_vs_16h_labels %>%
  filter(clean_comp_label != "")  # remove non-significant brackets

# 6. Add the Significance Brackets (Using the 5-row label data)
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

#CRISPR LEADING#####################################################################################
# **** Section 2: CRISPR strain RT-qPCR results ****
######################################################################################
######################################################################################

CRISPR_file_paths <- list(
  leading_transcripts_1 = "~/Desktop/qPCR data/additional candidate testing/CT_values_20260309_1.xlsx",
  leading_transcripts_2 = "~/Desktop/qPCR data/additional candidate testing/CT_values_20260310_2.xlsx",
  leading_transcripts_3 = "~/Desktop/qPCR data/additional candidate testing/CT_values_20260310_3.xlsx"
) 
# Read all files in one step - adding a new file is just one new line above
raw_data <- map(CRISPR_file_paths, read_excel)

# Groups: control, Nup211_863, Nup211_655
raw_data <- map(raw_data, ~ .x %>%
                  mutate(condition = case_when(
                    grepl("control",      sample) ~ "control",
                    grepl("C863",  sample) ~ "Nup211_863",
                    grepl("C655", sample) ~ "Nup211_655",
                    TRUE                          ~ NA_character_
                  ),
                  bio_rep = as.numeric(factor(gsub(".*_(\\d)$", "\\1", sample))))
)

###########################################
# Pipeline function - do this ONCE, use for all dataframes
###########################################
# NOTE: identical to run_qpcr_pipeline() above except for the condition names in the
# row_id case_when(). Kept separate deliberately; candidate for merging into one
# parameterised function during pipeline consolidation.
run_qpcr_pipeline_CRISPR <- function(df, comparisons, condition_levels) {
  
  # Groups: control, Nup211_863, Nup211_655
  
  # Technical replicate averaging
  tech_rep_avg <- df %>%
    group_by(sample, condition, bio_rep, gene) %>%
    summarise(mean_Ct = mean(Ct, na.rm = TRUE), .groups = "drop")
  
  # Format for downstream analysis
  ct_wide <- tech_rep_avg %>%
    mutate(row_id = paste(condition, bio_rep, sep = "_")) %>%
    select(row_id, gene, mean_Ct) %>%
    pivot_wider(names_from = gene, values_from = mean_Ct)
  
  # Delta Ct
  deltaCt_df <- ct_wide %>%
    pivot_longer(cols = -row_id, names_to = "gene", values_to = "Ct") %>%
    group_by(row_id) %>%
    mutate(act1_Ct = Ct[gene == "act1"],
           deltaCt = Ct - act1_Ct) %>%
    ungroup()
  
  # Mean delta Ct for controls
  mean_deltaCt_controls <- deltaCt_df %>%
    filter(grepl("control", row_id)) %>%
    group_by(gene) %>%
    summarise(mean_deltaCt_control = mean(deltaCt, na.rm = TRUE))
  
  # DeltaDeltaCt + fold changes
  df_FC <- deltaCt_df %>%
    left_join(mean_deltaCt_controls, by = "gene") %>%
    mutate(
      DeltaDeltaCt   = deltaCt - mean_deltaCt_control,
      FoldChange     = 2^(-DeltaDeltaCt),
      log2FoldChange = log2(FoldChange),
      condition = case_when(
        grepl("control",     row_id, ignore.case = TRUE) ~ "control",
        grepl("Nup211_863", row_id, ignore.case = TRUE) ~ "Nup211_863",
        grepl("Nup211_655",   row_id, ignore.case = TRUE) ~ "Nup211_655",
        TRUE ~ NA_character_
      ),
      condition = factor(condition, levels = condition_levels)
    )
  
  # Pairwise Welch's t-tests
  t_test_results <- df_FC %>%
    filter(!is.na(log2FoldChange)) %>%
    group_by(gene) %>%
    reframe(
      map_dfr(comparisons, function(comp) {
        grp1 <- log2FoldChange[condition == comp[1]]
        grp2 <- log2FoldChange[condition == comp[2]]
        
        if (length(grp1) < 2 | length(grp2) < 2) {
          return(tibble(group1 = comp[1], group2 = comp[2],
                        estimate = NA, conf.low = NA, conf.high = NA,
                        p.value = NA))
        }
        
        tidy(t.test(grp1, grp2)) %>%
          mutate(group1 = comp[1], group2 = comp[2])
      })
    ) %>%
    select(gene, group1, group2, estimate, conf.low, conf.high, p.value)
  
  # BH correction per comparison
  # NOTE: Alternative is to remove group_by(group1, group2) and correct
  # across all tests simultaneously for more conservative FDR control
  t_test_results <- t_test_results %>%
    group_by(group1, group2) %>%
    mutate(
      padj = p.adjust(p.value, method = "BH"),
      significant = case_when(
        padj < 0.001 ~ "***",
        padj < 0.01  ~ "**",
        padj < 0.05  ~ "*",
        TRUE          ~ "ns"
      )
    ) %>%
    ungroup()
  
  # y positions for significance brackets
  t_test_results <- t_test_results %>%
    group_by(gene) %>%
    mutate(
      max_log2FC = max(df_FC$log2FoldChange[df_FC$gene == cur_group()$gene],
                       na.rm = TRUE),
      y.position = max_log2FC + 0.5 * row_number()
    ) %>%
    ungroup()
  
  # Mean summary
  mean_values <- df_FC %>%
    group_by(gene, condition) %>%
    summarise(
      mean_log2FC = mean(log2FoldChange, na.rm = TRUE),
      mean_FC     = 2^(mean_log2FC),
      .groups = "drop"
    )
  
  # Return everything you need downstream as a named list
  return(list(
    df_FC          = df_FC,
    t_test_results = t_test_results,
    mean_values    = mean_values
  ))
}

###########################################
# Run pipeline on all dataframes at once
###########################################
condition_levels <- c("control", "Nup211_863", "Nup211_655")

comparisons <- list(
  c("control",     "Nup211_863"),
  c("control",     "Nup211_655"),
  c("Nup211_863", "Nup211_655")
)

# This runs the entire pipeline on every dataframe - adding leading_transcripts_4
# later just means adding one line to file_paths above, nothing else changes
CRISPR_results <- map(raw_data, ~ run_qpcr_pipeline_CRISPR(.x, comparisons, condition_levels))

names(CRISPR_results$leading_transcripts_1$df_FC)
names(CRISPR_results$leading_transcripts_1$mean_values)

################################################################################
# Pool results across all datasets
################################################################################

# 1. Replicate-level data (for individual points on boxplots)
################################################################################
pooled_replicates <- map_dfr(CRISPR_results, ~ .x$df_FC, .id = "dataset") %>%
  filter(gene != "act1") %>%
  select(dataset, gene, condition, row_id, log2FoldChange)

clean_pooled_replicates <- pooled_replicates %>%
  filter(condition != "control") %>%
  rename(log2FC = log2FoldChange)

# 2. Mean log2FC + 95% CIs per gene/condition
#    log2Upper/log2Lower are log2 transformed from conf.high/conf.low
################################################################################
# 2. Mean log2FC + 95% CIs per gene/condition
################################################################################
pooled_means <- map_dfr(CRISPR_results, ~ .x$df_FC, .id = "dataset") %>%
  filter(gene != "act1") %>%
  group_by(dataset, gene, condition) %>%
  summarise(
    mean_log2FC = mean(log2FoldChange, na.rm = TRUE),
    mean_FC     = 2^(mean_log2FC),
    n           = n(),
    se          = sd(log2FoldChange, na.rm = TRUE) / sqrt(n),
    ci          = qt(0.975, df = n - 1) * se,
    log2upper   = mean_log2FC + ci,
    log2lower   = mean_log2FC - ci,
    .groups = "drop"
  ) %>%
  left_join(
    map_dfr(CRISPR_results, ~ .x$t_test_results, .id = "dataset") %>%
      filter(group1 == "control") %>%
      rename(condition = group2) %>%
      select(dataset, gene, condition, padj, significant),
    by = c("dataset", "gene", "condition")
  ) 

clean_pooled_means <- pooled_means %>%
  filter(condition != "control") %>%
  mutate(condition = factor(condition, levels = c("control", "Nup211_863", "Nup211_655")))

# 3a. Control comparisons (control vs shutoff_16h, control vs FL_rescue)
################################################################################
pooled_sig_vs_control <- map_dfr(CRISPR_results, ~ .x$t_test_results, .id = "dataset") %>%
  filter(group1 == "control") %>%
  mutate(significant = ifelse(significant == "ns", "", significant)) %>%
  left_join(
    pooled_means %>% select(dataset, gene, condition, log2upper),
    by = c("dataset", "gene", "group2" = "condition")
  ) %>%
  group_by(dataset, gene) %>%
  mutate(y.position = max(log2upper, na.rm = TRUE) + 0.5 * row_number()) %>%
  ungroup() %>%
  select(dataset, gene, group1, group2, padj, significant, y.position) 

pooled_sig_vs_control <- pooled_sig_vs_control %>%
  filter(gene != "act1") %>%
  mutate(condition = factor(group2, levels = c("control", "Nup211_863", "Nup211_655")))

# 3b. Rescue vs shutoff comparison
################################################################################
pooled_sig_rescue_vs_shutoff <- map_dfr(CRISPR_results, ~ .x$t_test_results, .id = "dataset") %>%
  filter(group1 == "Nup211_863" & group2 == "Nup211_655") %>%
  mutate(significant = ifelse(significant == "ns", "", significant)) %>%
  left_join(
    pooled_means %>% select(dataset, gene, condition, log2upper),
    by = c("dataset", "gene", "group2" = "condition")
  ) %>%
  group_by(dataset, gene) %>%
  mutate(y.position = max(log2upper, na.rm = TRUE) + 0.5 * row_number()) %>%
  ungroup() %>%
  select(dataset, gene, group1, group2, padj, significant, y.position)

pooled_sig_rescue_vs_shutoff <- pooled_sig_rescue_vs_shutoff %>%
  filter(gene != "act1")

# Objects for plotting: 
clean_pooled_replicates
clean_pooled_means  
pooled_sig_vs_control
pooled_sig_rescue_vs_shutoff

#manually adjust positions like this : 
pooled_sig_vs_control <- pooled_sig_vs_control %>%
  mutate(y.position = if_else(gene == "swi6" & condition == "Nup211_655", 
                              1.00, 
                              y.position))

plot_leading_CRISPR <- ggplot(clean_pooled_replicates, aes(x = gene, y = log2FC, fill = condition)) +
  geom_hline(yintercept = 0, linewidth = 1.0, color = "black",linetype = 1) +
  
  # 2. Boxplot & Jitter (Using the 30-row raw data)
  geom_boxplot(aes(group = interaction(gene, condition)),
               position = position_dodge(width = 0.85),
               coef = 0, outlier.shape = NA, size = 0.8) +
  geom_jitter(aes(shape = condition),
              position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.85),
              size = 2.7, stroke = 1.2, fill = "white", color = "black") +
  
  # 3. Error Bars (Switching to the 10-row summary data)
  geom_errorbar(data = clean_pooled_means,
                aes(x = gene, ymin = log2lower, ymax = log2upper, group = condition),
                position = position_dodge(width = 0.85),
                width = 0.3, linewidth = 1, inherit.aes = FALSE) +
  
  # 4. Condition-specific Stars (Using the 10-row summary data)
  # LAYER A: The Numeric Mean (Black, bold, 2 sig figs)
  geom_text(data = clean_pooled_means,
            aes(x = gene, y = log2upper + 0.215 , label = sprintf("%.2f", mean_log2FC), group = condition),
            color = "black",
            position = position_dodge(0.85),
            size = 4.5, fontface = "bold.italic", vjust = 0, inherit.aes = FALSE) +
  
  # LAYER B: The Significance Stars (Black, bold)
  geom_text(data = pooled_sig_vs_control,
            aes(x = gene, y = y.position, label = significant, group = condition),
            color = "black",
            position = position_dodge(0.85),
            size = 8, fontface = "bold", vjust = 0, inherit.aes = FALSE) +
  
  # 5. Scales & Styling (From your template)
  scale_shape_manual(values = c("Nup211_863" = 21, "Nup211_655" = 24),
                     labels = c("<i>nup211<sub>1-863</sub></i>",
                                "<i>nup211<sub>1-655</sub></i>")) +
  scale_fill_manual(values = c("Nup211_863" = "thistle", "Nup211_655" = "sienna2"),
                    labels = c("<i>nup211<sub>1-863</sub></i>",
                               "<i>nup211<sub>1-655</sub></i>")) +
  scale_y_continuous(limits = c(-3, 5), breaks = seq(-3, 5, 1)) +
  labs(title = "RT-qPCR Results",
       subtitle =  "<b><i>nup211</i><sub><i>1-863</i></sub></b> <b>and</b> <b><i>nup211</i><sub><i>1-655</i></sub></b> <b>vs. wild-type</b>",
       x = NULL,
       y = bquote(bold(log[2] ~ "Fold Change")),
       shape = "Strain: ",
       fill = "Strain: ") +
  #theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_markdown(size = 16),
    axis.line = element_line(linewidth = 1.0),
    axis.ticks = element_line(linewidth = 1.2),
    axis.ticks.length = unit(3, "mm"),
    axis.text.x = element_text(size = 18, face = "bold", colour = 'black', angle = 30, hjust = 1),
    axis.text.y = element_text(size = 16, face = "bold", colour = 'black'),
    axis.title.y = element_text(size = 18, face = "bold", colour = 'black'),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = ggtext::element_markdown(size = 14, face = "bold"),
    legend.position = "bottom",
    legend.key.width = unit(1, "cm"),
    legend.spacing.x = unit(0.5, "cm")
  )

plot_leading_CRISPR

# 6. Add the Significance Brackets (Using the 5-row label data)

pooled_sig_rescue_vs_shutoff <- pooled_sig_rescue_vs_shutoff %>%
  filter(significant != "") %>%  # remove non-significant brackets
  left_join(
    clean_pooled_means %>%
      group_by(gene) %>%
      summarise(max_log2upper = max(log2upper, na.rm = TRUE)),
    by = "gene"
  ) %>%
  mutate(y.position = max_log2upper + 1.2)  # adjust offset as needed

plot_leading_CRISPR_2 <- plot_leading_CRISPR + 
  geom_signif(
    data = pooled_sig_rescue_vs_shutoff,
    aes(
      x          = gene,
      xmin       = as.numeric(factor(gene, levels = c("alm1", "bub1", "mad1", "not1", "swi6"))) - 0.22,
      xmax       = as.numeric(factor(gene, levels = c("alm1", "bub1", "mad1", "not1", "swi6"))) + 0.22,
      annotations = significant,
      fontface   = "bold",
      y_position = y.position
    ),
    tip_length  = 0.01,
    size        = 1,
    textsize    = 8,
    vjust       = 0.2,
    inherit.aes = FALSE,
    manual      = TRUE
  )

plot_leading_CRISPR_2

#CRISPR Supplemental###############################################################################
# Generate supplemental table (CRISPR strains vs 558 - wild type)
################################################################################
supp_leading_CRISPR <- pooled_replicates %>%
  filter(condition != "control") %>%
  group_by(gene, condition) %>%
  mutate(rep = row_number(),
         condition = factor(condition, levels = condition_levels)) %>%
  select(-c(row_id, dataset)) %>%
  ungroup() %>%
  pivot_wider(
    names_from  = c(condition, rep),
    values_from = log2FoldChange,
    names_glue  = "{condition}_rep{rep}_log2FC"
  ) %>% 
  select(c(1,5:7,2:4)) %>%
  
  # Join mean + CI
  left_join(
    clean_pooled_means %>%
      select(gene, condition, mean_log2FC, log2upper, log2lower) %>%
      pivot_wider(
        names_from  = condition,
        values_from = c(mean_log2FC, log2upper, log2lower),
        names_glue  = "{.value}_{condition}"
      ),
    by = "gene"
  ) %>%
  
  # Join padj and annotations for 16h vs control and rescue vs control
  left_join(
    pooled_sig_vs_control %>%
      select(gene, group2, padj, significant) %>%
      pivot_wider(
        names_from  = group2,
        values_from = c(padj, significant),
        names_glue  = "{.value}_{group2}_v_control"
      ),
    by = "gene"
  ) %>%
  
  # Join padj and annotation for rescue vs 16h
  left_join(
    pooled_sig_rescue_vs_shutoff %>%
      select(gene, padj, significant) %>%
      rename(
        padj_Nup211_863_vs_Nup211_655       = padj,
        annotation_Nup211_863_vs_Nup211_655 = significant
      ),
    by = "gene"
  ) %>%
  
  # Clean up column names to match your existing format
  rename_with(~ gsub("log2upper", "log2upperCI", .x)) %>%
  rename_with(~ gsub("log2lower", "log2lowerCI", .x))

existing_CRISPR_RTqPCR <- read_csv("~/Desktop/qPCR data/supplemental_files/supplemental_CRISPR_RTqPCR.csv")

# Columns in A but not in B
setdiff(names(existing_CRISPR_RTqPCR), names(supp_leading_CRISPR))

# Columns in B but not in A
setdiff(names(supp_leading_CRISPR), names(existing_CRISPR_RTqPCR))

# Fix mismatched names:

supp_leading_CRISPR_fixed <- supp_leading_CRISPR %>%
  rename(
    annotation_Nup211_863_v_control = significant_Nup211_863_v_control,
    annotation_Nup211_655_v_control = significant_Nup211_655_v_control
  )

setdiff(names(existing_CRISPR_RTqPCR), names(supp_leading_CRISPR_fixed))
setdiff(names(supp_leading_CRISPR_fixed), names(existing_CRISPR_RTqPCR))

supplemental_table_CRISPR_RTqPCR <- bind_rows(
  existing_CRISPR_RTqPCR,
  supp_leading_CRISPR_fixed
)

supplemental_table_CRISPR_RTqPCR

# Write to csv
# write.csv(supplemental_table_CRISPR_RTqPCR, "Documents/RNA-seq manuscript/supplemental_table_CRISPR_RTqPCR.csv", row.names = FALSE)

###############################################################################
#Merge with existing CRISPR Data                                              #
###############################################################################

CRISPR_PCR_table <- read_csv("~/Documents/RNA-seq manuscript/supplemental_table_CRISPR_RTqPCR.csv") 

rows_to_remove <- c("pho1", "urg2")

CRISPR_PCR_table <- CRISPR_PCR_table %>%
  filter(!(gene %in% rows_to_remove))

reps_long <- CRISPR_PCR_table %>%
  select(gene, 2:7) %>%
  pivot_longer(
    cols = -gene,
    names_to = c("condition", "replicate"),
    names_pattern = "Nup211_(\\d+)_?(rep\\d)_log2FC",
    values_to = "log2FC"
  )

# 2. Pivot the Metadata (padj and annotation)
meta_long <- CRISPR_PCR_table %>%
  select(gene, 14:17) %>%
  pivot_longer(
    cols = -gene,
    names_to = c(".value", "condition"),
    names_pattern = "(padj|annotation)_Nup211_(\\d+)_v_control"
  )

# 3. Join them together and keep the rest of the original table as a lookup
CRISPR_PCR_reps <- reps_long %>%
  left_join(meta_long, by = c("gene", "condition")) %>%
  mutate(condition = paste0("Nup211_", condition),
           gene = factor(gene, levels = gene_levels_CRISPR))

CRISPR_long_FC <- CRISPR_PCR_table %>%
  select(gene, 8:17) %>%
  mutate(gene = factor(gene, levels = gene_levels_CRISPR)) %>%
  # Rename columns to a consistent "Variable_Condition" format
  rename(
    log2FC_863     = mean_log2FC_Nup211_863,
    log2FC_655     = mean_log2FC_Nup211_655,
    log2upper_863  = log2upperCI_Nup211_863,
    log2upper_655  = log2upperCI_Nup211_655,
    log2lower_863  = log2lowerCI_Nup211_863,
    log2lower_655  = log2lowerCI_Nup211_655,
    padj_863       = padj_Nup211_863_v_control,
    padj_655       = padj_Nup211_655_v_control,
    annot_863      = annotation_Nup211_863_v_control,
    annot_655      = annotation_Nup211_655_v_control
  ) %>%
  # 2. Pivot into long format
  pivot_longer(
    cols = -gene,
    names_to = c(".value", "condition"),
    names_pattern = "(.*)_(863|655)"
  ) %>%
  # 3. Calculate plotting positions and clean labels
  mutate(
    clean_label = case_when(
      annot == "ns" ~ "",
      is.na(annot)  ~ "",
      TRUE          ~ annot
    ),
    mean_label = sprintf("%.2f", log2FC),
    mean_y_pos = log2upper + 0.25,
    star_y_pos = log2upper + 0.85
  )

# Create a separate dataframe for the 863 vs 655 comparison
comp_863_vs_655_labels <- CRISPR_PCR_table %>%
  # Select gene and the specific comparison columns (Indices 1, 18, 19)
  select(gene, 
         comp_padj = padj_Nup211_863_vs_Nup211_655, 
         comp_annot = annotation_Nup211_863_vs_Nup211_655) %>%
  
  # Ensure genes match your specific factor levels
  mutate(gene = factor(gene, levels = gene_levels_CRISPR)) %>%
  
  # Join with the max height from your previously created long dataframe
  left_join(
    CRISPR_long_FC %>% 
      group_by(gene) %>% 
      summarise(max_y = max(log2upper, na.rm = TRUE)),
    by = "gene"
  ) %>%
  
  # Set bracket and label logic
  mutate(
    # Set bracket position 1.6 units above the highest error bar (matching your template)
    y_bracket = max_y + 1.6,
    
    # Blank out 'ns' labels for the comparison
    clean_comp_label = case_when(
      comp_annot == "ns" ~ "",
      is.na(comp_annot)  ~ "",
      TRUE               ~ comp_annot
    )
  )

# Then to plot

# 1. Base Plot with Replicate Data

# Define the order once
target_order <- c("Nup211_863", "Nup211_655")

# Apply to all relevant dataframes
CRISPR_PCR_reps$condition <- factor(CRISPR_PCR_reps$condition, levels = target_order)

# The rows alternate: odd rows = Nup211_863, even rows = Nup211_655
CRISPR_long_FC$condition <- rep(c("Nup211_863", "Nup211_655"), times = nrow(CRISPR_long_FC)/2)
CRISPR_long_FC$condition <- factor(CRISPR_long_FC$condition, levels = target_order)

CRISPR_PCR_plot <- ggplot(CRISPR_PCR_reps, aes(x = gene, y = log2FC, fill = condition)) +
  geom_hline(yintercept = 0, linewidth = 1.0, color = "black") +
  
  # 2. Boxplot & Jitter (Individual Replicates)
  geom_boxplot(aes(group = interaction(gene, condition)),
               position = position_dodge(width = 0.85),
               coef = 0, outlier.shape = NA, size = 0.8) +
  geom_jitter(aes(shape = condition),
              position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.85),
              size = 2.7, stroke = 1.2, fill = "white", color = "black") +
  
  # 3. Error Bars (Using Summary Data)
  geom_errorbar(data = CRISPR_long_FC,
                aes(x = gene, ymin = log2lower, ymax = log2upper, group = condition),
                position = position_dodge(width = 0.85),
                width = 0.3, linewidth = 1, inherit.aes = FALSE) +
  
  # 4. Significance Stars & Means
  # geom_text(data = CRISPR_long_FC,
  #           aes(x = gene, y = mean_y_pos, label = mean_label, group = condition),
  #           color = "black", position = position_dodge(0.85),
  #           size = 4.5, fontface = "bold.italic", vjust = 0, inherit.aes = FALSE) +
  
  geom_text(data = CRISPR_long_FC,
            aes(x = gene, y = star_y_pos, label = clean_label, group = condition),
            color = "black", position = position_dodge(0.85),
            size = 8, fontface = "bold", vjust = 0, inherit.aes = FALSE) +
  
  # 5. Scales (Updated for 863/655)
  scale_shape_manual(
    values = c("Nup211_863" = 21, "Nup211_655" = 24),
    breaks = c("Nup211_863", "Nup211_655"), # Explicitly set order
    labels = c("Nup211<sub>1-863</sub>",
               "Nup211<sub>1-655</sub>")) +
  scale_fill_manual(
    values = c("Nup211_863" = "thistle", "Nup211_655" = "sienna2"),
    breaks = c("Nup211_863", "Nup211_655"), # Explicitly set order
    labels = c("Nup211<sub>1-863</sub>",
               "Nup211<sub>1-655</sub>")) +
  # 6. Styling
  labs(title = "CRISPR PCR Results",
       subtitle = bquote(bold("Nup211 truncations vs. control")),
       x = NULL, y = bquote(bold(log[2] ~ "Fold Change")),
       shape = "Truncation: ", fill = "Truncation: ") +
  scale_y_continuous(limits = c(-3, 4), breaks = seq(-3, 4, 1)) +
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

comp_863_vs_655_labels_clean <- comp_863_vs_655_labels %>%
  filter(clean_comp_label != "")  # remove non-significant brackets

# 7. Add Intergroup Brackets (863 vs 655)
CRISPR_PCR_plot_WIP <- CRISPR_PCR_plot + 
  geom_signif(
    data = comp_863_vs_655_labels_clean,
    aes(
      x = gene,
      xmin = as.numeric(gene) - 0.22, 
      xmax = as.numeric(gene) + 0.22, 
      annotations = clean_comp_label, 
      y_position = y_bracket 
    ),
    manual = TRUE, tip_length = 0.02, size = 1, textsize = 8, inherit.aes = FALSE
  )

CRISPR_PCR_plot_WIP

