# TurboID summary tables for thesis: 

install.packages("flexlsx")
library(flexlsx)
library(readxl)
library(readr)
library(tidyr)
library(dplyr)
library(stringr)
library(writexl)
library("openxlsx2")


TurboID_Summary <- read_excel("~/Documents/TurboID Mass Spec Analysis/TurboID_summary.xlsx")

names(TurboID_Summary) <- c("protein",
                         "description",
                         "dis3_mean_abundance",
                         "nup211_mean_abundance",
                         "untagged_mean_abundance",
                         "dis3_peptides",
                         "nup211_peptides",
                         "untagged_peptides",
                         "detected_dis3",
                         "detected_nup211",
                         "detected_untagged")



TurboID_Coverage <- read_csv("~/Documents/TurboID Mass Spec Analysis/protein_coverage_results.csv")                       


TurboID_Coverage <- TurboID_Coverage %>% select(c(1,2,5))


wide_TurboID_Coverage <- TurboID_Coverage %>%
  group_by(protein) %>%
  # Pivot to wide format for easy comparison
  pivot_wider(
    id_cols = c(protein),
    names_from = experiment,
    values_from = c(coverage_percent),
    values_fill = list(coverage_percent = 0)
  )

names(wide_TurboID_Coverage) <- c("protein",
                            "untagged_coverage",
                            "dis3_coverage",
                            "nup211_coverage")


TurboID_comprehensive <- TurboID_Summary %>%
  left_join(wide_TurboID_Coverage, by = c("protein")) %>%
  arrange(desc(nup211_mean_abundance))


TurboID_comprehensive <- TurboID_comprehensive %>% select(c(1,2,4,7,14,3,6,13,5,8,12,10,9,11))


TurboID_comprehensive <- TurboID_comprehensive %>%
  arrange(desc(nup211_mean_abundance))


library(scales)
TurboID_comprehensive <- TurboID_comprehensive %>%
  mutate(across(c(3,6,9), ~ .x * 100))

library(dplyr)
library(flextable)

ft <- flextable(TurboID_comprehensive)

ft <- colformat_num(
  ft,
  j = c("dis3_mean_abundance",
        "nup211_mean_abundance",
        "untagged_mean_abundance",
        "untagged_coverage",
        "dis3_coverage",
        "nup211_coverage"),
  digits = 2
)


# Rename visible column labels
ft <- set_header_labels(
  ft,
  protein = "Protein",
  description = "Description",
  nup211_mean_abundance = "Mean relative abundance (%)",
  nup211_peptides = "# of peptides",
  nup211_coverage = "Protein coverage (%)",
  dis3_mean_abundance = "Mean relative abundance (%)",
  dis3_peptides = "# peptides",
  dis3_coverage = "Protein coverage (%)",
  untagged_mean_abundance = "Mean relative abundance (%)",
  untagged_peptides = "# peptides",
  untagged_coverage = "Protein coverage (%)",
  detected_nup211 = "Nup211-TurboID",
  detected_dis3 = "Dis3-TurboID",
  detected_untagged = "Untagged control"
)

# Add grouped headers
ft <- add_header_row(
  ft,
  values = c("", "", "Nup211-TurboID", "Dis3-TurboID","Untagged control", "Detected in:"),
  colwidths = c(1, 1, 3, 3, 3, 3)
)

# Formatting tweaks
ft <- align(ft, align = "center", part = "all")
ft <- align(ft, j = c("protein", "description"), align = "left")
ft <- autofit(ft)
band_rows <- seq(2, nrow(summary_table), by = 2)

# Optional: nicer borders
ft <- theme_booktabs(ft)

ft <- align(
  ft,
  i = 1,
  align = "center",
  part = "header"
)



ft <- bg(
  ft,
  i = seq(1, nrow(Summary_Full), by = 2),
  bg = "#F2F2F2",
  part = "body"
)

library(officer)
library(flextable)
vline <- fp_border(color = "black", width = 0.75)

ft2 <- flextable::border(
  ft,
  j = c("description","nup211_mean_abundance","dis3_mean_abundance","untagged_mean_abundance",
        "detected_nup211"),
  border.left = vline,
  part = "all"
)

ft3 <- colformat_num(
  ft2,
  digits = 2
)

str(Summary_Full$nup211_coverage)

# save_as_docx(ft3, path = "Nup211_TurboID_Thesis_Table.docx")


library("flexlsx")
wb <- openxlsx2::wb_workbook()$add_worksheet("TurboID_Summary")

# add the flextable ft to the workbook, sheet "TurboID_Summary"
# offset the table to cell 'C2'
wb <- wb_add_flextable(wb, "TurboID_Summary", ft3, dims = "C2")

# save the workbook to a temporary xlsx file
tmpfile <- tempfile(fileext = ".xlsx")
wb$save(tmpfile)

#wb$save("~/Documents/TurboID_thesis.xlsx")



# Then separate comprehensive sets


nup211_all <- TurboID_comprehensive %>%
  filter(detected_nup211)

dis3_all <- TurboID_comprehensive %>%
  filter(detected_dis3)

untagged_all <- TurboID_comprehensive %>%
  filter(detected_untagged)


nup211_specific <- TurboID_comprehensive %>%
  filter(detected_nup211 & !detected_untagged)

dis3_specific <- TurboID_comprehensive %>%
  filter(detected_dis3 & !detected_untagged)




nup211_only <- TurboID_comprehensive %>%
  filter(detected_nup211, !detected_dis3, !detected_untagged)

nup211_dis3 <- TurboID_comprehensive %>%
  filter(detected_nup211, detected_dis3)


nup211_dis3_only <- TurboID_comprehensive %>%
  filter(detected_nup211, detected_dis3, !detected_untagged)


dis3_only <- TurboID_comprehensive %>%
  filter(!detected_nup211, detected_dis3, !detected_untagged)

untagged_only <- TurboID_comprehensive %>%
  filter(detected_untagged, !detected_nup211, !detected_dis3)

# untagged only, nup211 only, and dis3 only datasets







##################################################################################

library(flextable)
library(dplyr)
library(officer)

make_turboid_flextable <- function(df, band_rows = NULL) {
  
  ft <- flextable(df)
  
  # Numeric formatting
  num_cols <- intersect(
    c(
      "dis3_mean_abundance",
      "nup211_mean_abundance",
      "untagged_mean_abundance",
      "untagged_coverage",
      "dis3_coverage",
      "nup211_coverage"
    ),
    colnames(df)
  )
  
  ft <- colformat_num(ft, j = num_cols, digits = 2)
  
  # Header labels (only apply to columns that exist)
  header_map <- list(
    protein = "Protein",
    description = "Description",
    nup211_mean_abundance = "Mean relative abundance (%)",
    nup211_peptides = "# of peptides",
    nup211_coverage = "Protein coverage (%)",
    dis3_mean_abundance = "Mean relative abundance (%)",
    dis3_peptides = "# peptides",
    dis3_coverage = "Protein coverage (%)",
    untagged_mean_abundance = "Mean relative abundance (%)",
    untagged_peptides = "# peptides",
    untagged_coverage = "Protein coverage (%)",
    detected_nup211 = "Nup211-TurboID",
    detected_dis3 = "Dis3-TurboID",
    detected_untagged = "Untagged control"
  )
  
  header_map <- header_map[names(header_map) %in% colnames(df)]
  ft <- do.call(set_header_labels, c(list(ft), header_map))
  
  # Grouped header (only if full table)
  if (all(c("nup211_mean_abundance","dis3_mean_abundance","untagged_mean_abundance") %in% colnames(df))) {
    ft <- add_header_row(
      ft,
      values = c("", "", "Nup211-TurboID", "Dis3-TurboID", "Untagged control", "Detected in:"),
      colwidths = c(1, 1, 3, 3, 3, 3)
    )
  }
  
  # Alignment & layout
  ft <- align(ft, align = "center", part = "all")
  ft <- align(ft, j = intersect(c("protein", "description"), colnames(df)), align = "left")
  ft <- align(ft, i = 1, align = "center", part = "header")
  ft <- autofit(ft)
  ft <- theme_booktabs(ft)
  
  # Row banding
  if (!is.null(band_rows)) {
    bg_cols <- colnames(df)   # include protein too
    ft <- bg(
      ft,
      i = band_rows,
      j = bg_cols,
      bg = "#F2F2F2",
      part = "body"
    )
  }
  
  # Vertical separators
  vline <- fp_border(color = "black", width = 0.75)
  
  sep_cols <- intersect(
    c(
      "description",
      "nup211_mean_abundance",
      "dis3_mean_abundance",
      "untagged_mean_abundance",
      "detected_nup211"
    ),
    colnames(df)
  )
  
  ft <- flextable::border(
    ft,
    j = sep_cols,
    border.left = vline,
    part = "all"
  )
  
  ft
}

ft1 <- make_turboid_flextable(
  TurboID_comprehensive,
  band_rows = seq(1, nrow(TurboID_comprehensive), by = 2)
)


ft2 <- make_turboid_flextable(
  untagged_only,
  band_rows = seq(1, nrow(untagged_only), by = 2)
)

ft3 <- make_turboid_flextable(
  dis3_only,
  band_rows = seq(1, nrow(dis3_only), by = 2)
)

ft4 <- make_turboid_flextable(
  nup211_only,
  band_rows = seq(1, nrow(nup211_only), by = 2)
)

ft5 <- make_turboid_flextable(
  nup211_all,
  band_rows = seq(1, nrow(nup211_all), by = 2)
)

ft6 <- make_turboid_flextable(
  nup211_dis3_only,
  band_rows = seq(1, nrow(nup211_dis3_only), by = 2)
)

ft7 <- make_turboid_flextable(
  untagged_all,
  band_rows = seq(1, nrow(untagged_all), by = 2)
)





fts <- list(
  specific_Nup211 = ft4,
  specific_Dis3 = ft3,
  specific_untagged = ft2,
  all_Nup211 = ft5,
  Nup211_Dis3 = ft6,
  all_untagged = ft7
)

wb <- wb_workbook()

for (nm in names(fts)) {
  wb <- wb$add_worksheet(nm)
  wb <- wb_add_flextable(wb, nm, fts[[nm]], dims = "B2")
  
  wb$set_col_widths(
    sheet = nm,
    cols = 1:50,
    widths = "auto"
  )
}

# wb$save("~/Documents/TurboID_thesis_tables.xlsx")



