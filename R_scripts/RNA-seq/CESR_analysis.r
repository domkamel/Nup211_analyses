############################################################### 
# Stress response gene analysis ------------------------
# Comparing nup211-so (16h) and (48h) vs CESR gene lists
###############################################################

library(readr)
library(dplyr)
library(tidyverse)
library(readxl)
library(ggtext)
library(patchwork)
library(ggplot2)


# Load data ---------------------------------------------------------------

# Helper to load CESR gene lists
load_cesr <- function(path) {
  read_delim(path, delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>%
    mutate(gene_name = ifelse(is.na(`Gene name`), `Systematic ID`, `Gene name`)) %>%
    select(c(1, 4, 3)) %>%
    rename(Gene_ID = `Systematic ID`, gene_product = `Product description`)
}

cesr_induced_df   <- load_cesr("Documents/RNA-seq manuscript/published_gene_lists/cesr_induced_pombase.tsv")
cesr_repressed_df <- load_cesr("Documents/RNA-seq manuscript/published_gene_lists/cesr_repressed_pombase.tsv")

# Pull just the gene names for downstream filtering
cesr_induced   <- cesr_induced_df   %>% pull(gene_name)
cesr_repressed <- cesr_repressed_df %>% pull(gene_name)

# DESeq2 results
res_16h <- read_csv("Rstudio/expression_analysis_16h/results/DE_spreadsheets/DESeq2_results_protein_coding.csv") %>%
  column_to_rownames(var = names(.)[1])

res_48h <- read_csv("Desktop/nup211_GEO_sumbission_20221014/processed data files/DEseq2Results.csv") %>%
  column_to_rownames(var = names(.)[1]) %>%
  filter(Type == "protein coding gene") %>%
  rename(gene_type = Type)


# Define DEG sets and backgrounds -----------------------------------------

padj_threshold   <- 1e-5
log2fc_threshold <- 1

# Backgrounds: all genes tested per timepoint
bg_16h <- res_16h$gene_name
bg_48h <- res_48h$gene_name

# DEG categories per timepoint
deg_16h_up   <- res_16h %>% filter(padj < padj_threshold, log2FoldChange >=  log2fc_threshold) %>% pull(gene_name)
deg_16h_down <- res_16h %>% filter(padj < padj_threshold, log2FoldChange <= -log2fc_threshold) %>% pull(gene_name)
deg_16h_ns   <- setdiff(bg_16h, c(deg_16h_up, deg_16h_down))

deg_48h_up   <- res_48h %>% filter(padj < padj_threshold, log2FoldChange >=  log2fc_threshold) %>% pull(gene_name)
deg_48h_down <- res_48h %>% filter(padj < padj_threshold, log2FoldChange <= -log2fc_threshold) %>% pull(gene_name)
deg_48h_ns   <- setdiff(bg_48h, c(deg_48h_up, deg_48h_down))


# Fisher's exact tests ----------------------------------------------------

run_fisher_2x2 <- function(query_genes, cesr_genes, background_genes) {
  query_genes <- intersect(query_genes, background_genes)
  cesr_genes  <- intersect(cesr_genes,  background_genes)
  
  a <- length(intersect(query_genes, cesr_genes))
  b <- length(setdiff(query_genes, cesr_genes))
  c <- length(setdiff(cesr_genes,  query_genes))
  d <- length(background_genes) - a - b - c
  
  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2))
  
  data.frame(
    odds_ratio = unname(ft$estimate),
    ci_low     = ft$conf.int[1],
    ci_high    = ft$conf.int[2],
    p_value    = ft$p.value
  )
}

# Build the comparison grid efficiently rather than 12 lines of bind_rows code lines
comparison_grid <- expand.grid(
  timepoint = c("16h", "48h"),
  deg       = c("Up", "Down", "NS"),
  cesr      = c("Induced", "Repressed"),
  stringsAsFactors = FALSE
)

# Lookup tables for the gene sets
deg_lookup <- list(
  "16h_Up" = deg_16h_up, "16h_Down" = deg_16h_down, "16h_NS" = deg_16h_ns,
  "48h_Up" = deg_48h_up, "48h_Down" = deg_48h_down, "48h_NS" = deg_48h_ns
)
cesr_lookup <- list("Induced" = cesr_induced, "Repressed" = cesr_repressed)
bg_lookup   <- list("16h" = bg_16h, "48h" = bg_48h)

results <- comparison_grid %>%
  rowwise() %>%
  mutate(
    fisher = list(run_fisher_2x2(
      deg_lookup[[paste(timepoint, deg, sep = "_")]],
      cesr_lookup[[cesr]],
      bg_lookup[[timepoint]]
    ))
  ) %>%
  unnest(fisher) %>%
  ungroup() %>%
  mutate(
    p_adj      = p.adjust(p_value, method = "BH"),
    log2OR     = log2(odds_ratio),
    log2_lower = log2(ci_low),
    log2_upper = log2(ci_high),
    label_sig  = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    ),
    timepoint = factor(timepoint, levels = c("16h", "48h")),
    deg       = factor(deg,       levels = c("Up", "Down", "NS")),
    cesr      = factor(cesr,      levels = c("Induced", "Repressed"))
  )

print(results)


# Proportion data for % CESR affected plot --------------------------------

build_prop_df <- function(res_df, timepoint_label) {
  res_df %>%
    mutate(
      cesr_status = case_when(
        gene_name %in% cesr_induced   ~ "Induced",
        gene_name %in% cesr_repressed ~ "Repressed"
      ),
      de_cat = case_when(
        padj < padj_threshold & log2FoldChange >=  log2fc_threshold ~ "Up",
        padj < padj_threshold & log2FoldChange <= -log2fc_threshold ~ "Down",
        TRUE                                                          ~ "NS"
      ),
      timepoint = timepoint_label
    ) %>%
    filter(!is.na(cesr_status)) %>%
    select(gene_name, cesr_status, de_cat, timepoint)
}

prop_df <- bind_rows(
  build_prop_df(res_16h, "16h"),
  build_prop_df(res_48h, "48h")
) %>%
  mutate(
    timepoint   = factor(timepoint,   levels = c("16h", "48h")),
    cesr_status = factor(cesr_status, levels = c("Induced", "Repressed")),
    de_cat      = factor(de_cat,      levels = c("Up", "Down", "NS"))
  ) %>%
  dplyr::count(timepoint, cesr_status, de_cat) %>%
  group_by(timepoint, cesr_status) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()



# Define theme and aesthetics for plotting -----------------------------------


cesr_fill_scale <- scale_fill_manual(
  name   = "***nup211-so***",
  values = c("16h" = "#1b9e77", "48h" = "#D55E00")
)

cesr_theme <- theme_grey(base_size = 14) +
  theme(
    plot.title        = element_text(size = 18, face = "bold"),
    plot.subtitle     = element_text(size = 16, face = "bold"),
    axis.line         = element_blank(),
    axis.ticks        = element_line(size = 1.2),
    axis.ticks.length = unit(3, "mm"),
    panel.grid.major.y = element_line(color = "white", linewidth = 0.5, linetype = 1),
    panel.grid.minor.y = element_line(color = "white", linewidth = 0.5, linetype = 1),
    panel.grid.major.x = element_blank(),
    axis.text.x       = element_text(size = 18, face = "bold", colour = "black",
                                     angle = 0, hjust = 0.5),
    axis.text.y       = element_text(size = 16, face = "bold", colour = "black"),
    axis.title.y      = element_markdown(family = "Arial", face = "bold", size = 20),
    axis.title.x      = element_text(size = 18, face = "bold", colour = "black",
                                     margin = margin(t = 20)),
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 2),
    strip.background  = element_rect(fill = "black", colour = "black"),
    strip.text        = element_text(size = 16, face = "bold", colour = "white"),
    panel.spacing     = unit(2, "lines"),
    legend.title      = element_markdown(size = 16),
    legend.text       = element_markdown(size = 16, face = "bold"),
    legend.position   = "bottom",
    legend.key.width  = unit(1, "cm"),
    legend.spacing.x  = unit(0.5, "cm")
  )


# % CESR affected plot ----------------------------------------------------
cesr.labs <- c("Induced CESR genes", "Repressed CESR genes")
names(cesr.labs) <- c("Induced", "Repressed")

prop_plot <- ggplot(prop_df, aes(x = de_cat, y = proportion, fill = timepoint)) +
  geom_col(position = position_dodge(width = 0.78),
           width = 0.65, color = "black", linewidth = 0.8) +
  geom_text(aes(label = ifelse(proportion >= 0.002,
                               scales::percent(proportion, accuracy = 0.1),
                               "")),
            position = position_dodge(width = 0.78),
            vjust = -0.5, size = 4.5, fontface = "bold") +
  facet_wrap(~ cesr_status, labeller = labeller(cesr_status = cesr.labs)) +
  cesr_fill_scale +
  scale_y_continuous(labels = scales::percent,
                     expand = expansion(mult = c(0, 0.12))) +
  labs(
    #title = "Distribution of induced and repressed CESR genes",
    x     = "Expression change",
    y     = "Percentage (%)" # % of CESR genes or Percent affected
  ) +
  cesr_theme


# Fisher's Odds Ratio plot  -----------------------------------------------
cesr.labs <- c("Induced CESR genes", "Repressed CESR genes")
names(cesr.labs) <- c("Induced", "Repressed")

or_plot <- ggplot(results, aes(x = deg, y = log2OR, fill = timepoint)) +
  geom_hline(yintercept = 0, linewidth = 1.0, color = "black") +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.65, color = "black", linewidth = 0.8) +
  geom_errorbar(aes(ymin = log2_lower, ymax = log2_upper),
                position = position_dodge(width = 0.7),
                width = 0.25, linewidth = 1) +
  geom_text(aes(
    label = ifelse(label_sig == "ns", "", label_sig),
    y = ifelse(log2OR >= 0, log2_upper + 0.35, log2_lower - 0.5)
  ),
  position = position_dodge(width = 0.7),
  color = "black", size = 8, fontface = "bold") +
  facet_wrap(~ cesr, labeller = labeller(cesr = cesr.labs)) +
  cesr_fill_scale +
  scale_y_continuous(breaks = -4:3) +
  coord_cartesian(ylim = c(-4.5, 2.8)) +
  labs(
    #title = "CESR enrichment analysis",
    x     = "Expression change",
    y     = "log<sub>2</sub> (Odds Ratio)"
  ) +
  cesr_theme


# Combining plots for figure ----------------------------------------------

combined <- prop_plot / or_plot +
  plot_layout(axis_titles = "collect",
              guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "right",
    plot.tag = element_text(size = 20, face = "bold")
  )

print(combined)


# Exporting Results -------------------------------------------------------

library(flextable)
library(officer)

# ---- Helper to format p-values with superscript exponents ------------------
format_pval <- function(x) {
  sci      <- formatC(x, format = "e", digits = 2)
  mantissa <- sub("e.*", "", sci)
  exponent <- as.integer(sub(".*e", "", sci))
  list(mantissa = mantissa, exponent = as.character(exponent))
}

# ---- Build a clean results table from the analysis dataframe ---------------
# Includes percentages by merging in prop_df
deg_sizes <- tibble(
  timepoint = c("16h","16h","16h","48h","48h","48h"),
  deg       = c("Up","Down","NS","Up","Down","NS"),
  DEG_n     = c(length(deg_16h_up), length(deg_16h_down), length(deg_16h_ns),
                length(deg_48h_up), length(deg_48h_down), length(deg_48h_ns))
) %>%
  mutate(
    timepoint = factor(timepoint, levels = c("16h", "48h")),
    deg       = factor(deg,       levels = c("Up", "Down", "NS"))
  )

results_table <- results %>%
  left_join(
    prop_df %>%
      rename(deg = de_cat, cesr = cesr_status) %>%
      select(timepoint, cesr, deg, n, proportion),
    by = c("timepoint", "cesr", "deg")
  ) %>%
  left_join(deg_sizes, by = c("timepoint", "deg")) %>%
  arrange(timepoint, cesr, deg) %>%        # sort factors
  transmute(
    Timepoint    = as.character(timepoint),
    CESR_set     = as.character(cesr),
    DEG          = as.character(deg),
    DEG_n        = DEG_n,                 # total genes in this DEG category
    Overlap_n    = n,                     # CESR genes in this category
    Percent      = scales::percent(proportion, accuracy = 0.1),
    Odds_Ratio   = round(odds_ratio, 2),
    CI_low       = round(ci_low, 2),
    CI_high      = round(ci_high, 2),
    log2OR       = round(log2OR, 2),
    log2_lower   = round(log2_lower, 2),
    log2_upper   = round(log2_upper, 2),
    p_value      = signif(p_value, 3),
    p_adj        = signif(p_adj, 3),
    Significance = label_sig
  )

print(results_table)

# ---- Save CSV --------------------------------------------------------------
# write.csv(results_table,
#           "Documents/RNA-seq manuscript/Transcriptome_CESR_analysis_2x3.csv",
#           row.names = FALSE)

# ---- Build flextable for Word export ---------------------------------------
pval_parts <- format_pval(results_table$p_value)
padj_parts <- format_pval(results_table$p_adj)

ft <- flextable(results_table) %>%
  set_header_labels(
    Timepoint    = "Timepoint",
    CESR_set     = "CESR set",
    DEG          = "DEG category",
    DEG_n        = "# of DEGs",
    CESR_n       = "# of CESR genes",
    Percent      = "% of CESR set",
    Odds_Ratio   = "Odds Ratio",
    CI_low       = "CI low",
    CI_high      = "CI high",
    log2OR       = "log\u2082 OR",
    log2_lower   = "log\u2082 CI low",
    log2_upper   = "log\u2082 CI high",
    p_value      = "p-value",
    p_adj        = "BH-adjusted p",
    Significance = "Significance"
  ) %>%
  bold(part = "header") %>%
  merge_v(j = c("Timepoint", "CESR_set")) %>%
  valign(j = c("Timepoint", "CESR_set"), valign = "top") %>%
  bold(j = c("Timepoint", "CESR_set")) %>%
  autofit()

# Format p-values with superscript exponents
for (i in seq_len(nrow(results_table))) {
  ft <- compose(
    ft, i = i, j = "p_value",
    value = as_paragraph(
      as_chunk(paste0(pval_parts$mantissa[i], " \u00d7 10")),
      as_chunk(pval_parts$exponent[i],
               props = fp_text(vertical.align = "superscript", font.size = 12))
    )
  )
  ft <- compose(
    ft, i = i, j = "p_adj",
    value = as_paragraph(
      as_chunk(paste0(padj_parts$mantissa[i], " \u00d7 10")),
      as_chunk(padj_parts$exponent[i],
               props = fp_text(vertical.align = "superscript", font.size = 12))
    )
  )
}

# make fit to landscape
landscape_section <- prop_section(
  page_size = page_size(orient = "landscape"),
  type = "continuous"
)

#  Save Word doc _---------------------------------------------------------

# fit to landscape width
ft <- ft %>% fit_to_width(max_width = 9.5)

doc <- read_docx() %>%
  body_add_flextable(ft) %>%
  body_end_block_section(block_section(landscape_section))
# print(doc, target = "Documents/RNA-seq manuscript/Transcriptome_CESR_analysis_2x3.docx")


# Contingency plots similar to what Weigang mentioned ---------------------

# 16h contingency table
tab_16h <- results_table %>%
  filter(Timepoint == "16h") %>%
  mutate(cell = paste0(Overlap_n, " (", Percent, ")")) %>%
  select(DEG, CESR_set, cell) %>%
  pivot_wider(names_from = CESR_set, values_from = cell) %>%
  rename(`Induced CESR` = Induced, `Repressed CESR` = Repressed) %>%
  arrange(factor(DEG, levels = c("Up", "Down", "NS")))

# 48h contingency table
tab_48h <- results_table %>%
  filter(Timepoint == "48h") %>%
  mutate(cell = paste0(Overlap_n, " (", Percent, ")")) %>%
  select(DEG, CESR_set, cell) %>%
  pivot_wider(names_from = CESR_set, values_from = cell) %>%
  rename(`Induced CESR` = Induced, `Repressed CESR` = Repressed) %>%
  arrange(factor(DEG, levels = c("Up", "Down", "NS")))

print(tab_16h)
print(tab_48h)
