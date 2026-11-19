library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggrepel)

acute_res_df <- readxl::read_xlsx("~/shared/PRJ_TONGTA/prj_atac_Lyn/dat/TOBIAS_bindetect/Acute_vs_Control/bindetect_results.xlsx")
res <- readxl::read_xlsx("~/shared/PRJ_TONGTA/prj_atac_Lyn/dat/TOBIAS_bindetect/Chronic_vs_Control/bindetect_results.xlsx")
res <- readxl::read_xlsx("~/shared/PRJ_TONGTA/prj_atac_Lyn/dat/TOBIAS_bindetect/Chronic_vs_Acute/bindetect_results.xlsx")

comparison = "Acute_Control"

acute_res_df <-
  acute_res_df %>%
  mutate(
    case_mean_score = .data[["Acute_mean_score"]],
    case_bound = .data[["Acute_bound"]],
    control_mean_score = .data[["Control_mean_score"]],
    control_bound = .data[["Control_bound"]],
    change = .data[[paste0(comparison, "_change")]],
    pvalue = as.double(.data[[paste0(comparison, "_pvalue")]]),
    highlighted = .data[[paste0(comparison, "_highlighted")]],
    padj = p.adjust(pvalue, method = "BH")
  ) %>%
  arrange(desc(abs(change))) %>%
  dplyr::select(-Acute_mean_score, -Acute_bound, -Control_mean_score, -Control_bound,
                -Acute_Control_change, -Acute_Control_pvalue, -Acute_Control_highlighted)

acute_res_df <- res
chronic_res_df <- res
chronic_acute_df <- res

writexl::write_xlsx(list(acute_res_df, chronic_res_df, chronic_acute_df),
                    "~/shared/PRJ_TONGTA/prj_atac_Lyn/res/tobias_res/TOBIAS_results_20261005.xlsx")


filter_res <-
  res %>%
  mutate(
    change = .data[[paste0(comparison, "_change")]],
    pvalue = as.double(.data[[paste0(comparison, "_pvalue")]]),
    padj = p.adjust(pvalue, method = "BH"),
    neg_log10_padj = -log10(padj),
    sig = if_else(padj < 0.05 & change >= 0.1, "up",
                  if_else(padj < 0.05 & change <= -0.1, "down", "ns")),
    gene_labels = if_else(sig == "ns", NA_character_, motif_id)
  )

case_label = "Acute"
ctrl_label = "Control"

volcano_plot <-
  filter_res %>%
  ggplot(., aes(x = change, y = neg_log10_padj, color = sig)) +
  geom_point(alpha = 0.6) +
  geom_text_repel(
    aes(label = gene_labels), size = 3
  ) +
  geom_vline(xintercept = c(-0.1, 0.1), color = "darkred",
             linewidth = 0.4, linetype = "dashed") +
  scale_color_manual(
    values = c(down = "#1b7837", ns = "#898781", up = "#762a83"),
    labels = c(down = ctrl_label, ns = "Not significant", up = case_label),
    name = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    #axis.line = element_line(color = "#c3c2b7"),
    legend.position = "top"
  ) +
  labs(x = paste("Differential binding score:", paste(case_label, "VS", ctrl_label)),
       y = "-log10(FDR)",
       title = paste("TOBIAS BINDetect:", paste(case_label, "VS", ctrl_label)))


acute_ctrl_plot <- volcano_plot
chronic_ctrl_plot <- volcano_plot
chronic_acute_plot <- volcano_plot


pdf("~/shared/PRJ_TONGTA/prj_atac_Lyn/res/Figures_Tables/volcano_plots_differential_TOBIAS_binding_score_20261002.pdf", width = 7, height = 6)
plot(acute_ctrl_plot)
plot(chronic_ctrl_plot)
plot(chronic_acute_plot)
dev.off()

#
# Comparison with Chromvar results
#

acute_res_chromvar <- readxl::read_xlsx("~/shared/PRJ_TONGTA/prj_atac_Lyn/res/Figures_Tables/differential_motif_results_20260921.xlsx",
                                        sheet = 1)
chronic_res_chromvar <- readxl::read_xlsx("~/shared/PRJ_TONGTA/prj_atac_Lyn/res/Figures_Tables/differential_motif_results_20260921.xlsx",
                                        sheet = 2)
chronic_acute_res_chromvar <- readxl::read_xlsx("~/shared/PRJ_TONGTA/prj_atac_Lyn/res/Figures_Tables/differential_motif_results_20260921.xlsx",
                                        sheet = 3)


sig_chromvar <-
  chronic_res_chromvar %>%
  filter(abs(logFC) >= 0.5 & adj.P.Val <= 0.1) %>%
  mutate(ChromVAR = "TRUE") %>%
  dplyr::rename("motif_ID" = "ID",
                "gene_symbol" = "symbol")

sig_tobias <- acute_res_df %>%
  mutate(TOBIAS_top_10perc = ifelse(top_10_quantile == "TRUE", "TRUE", "FALSE")) %>%
  mutate(TOBIAS_highlighted = ifelse(Acute_Control_highlighted == TRUE, "TRUE", "FALSE")) %>%
  filter(TOBIAS_top_10perc == "TRUE" | TOBIAS_highlighted == "TRUE") %>%
  dplyr::rename("motif_ID" = "name",
                "gene_symbol" = "motif_id")


sig_acute_res <- full_join(sig_chromvar, sig_tobias) %>%
  dplyr::select(motif_ID, gene_symbol, ChromVAR, TOBIAS_top_10perc, TOBIAS_highlighted) %>%
  mutate(ChromVAR = ifelse(is.na(ChromVAR), "FALSE", ChromVAR),
         TOBIAS_top_10perc = ifelse(is.na(TOBIAS_top_10perc), "FALSE", TOBIAS_top_10perc),
         TOBIAS_highlighted = ifelse(is.na(TOBIAS_highlighted), "FALSE", TOBIAS_highlighted)) %>%
  arrange(desc(ChromVAR), desc(TOBIAS_top_10perc), desc(TOBIAS_highlighted))

sig_chronic_res <- full_join(sig_chromvar, sig_tobias) %>%
  dplyr::select(motif_ID, gene_symbol, ChromVAR, TOBIAS_top_10perc, TOBIAS_highlighted) %>%
  mutate(ChromVAR = ifelse(is.na(ChromVAR), "FALSE", ChromVAR),
         TOBIAS_top_10perc = ifelse(is.na(TOBIAS_top_10perc), "FALSE", TOBIAS_top_10perc),
         TOBIAS_highlighted = ifelse(is.na(TOBIAS_highlighted), "FALSE", TOBIAS_highlighted)) %>%
  arrange(desc(ChromVAR), desc(TOBIAS_top_10perc), desc(TOBIAS_highlighted))


sig_chronic_acute_res <- full_join(sig_chromvar, sig_tobias) %>%
  dplyr::select(motif_ID, gene_symbol, ChromVAR, TOBIAS_top_10perc, TOBIAS_highlighted) %>%
  mutate(ChromVAR = ifelse(is.na(ChromVAR), "FALSE", ChromVAR),
         TOBIAS_top_10perc = ifelse(is.na(TOBIAS_top_10perc), "FALSE", TOBIAS_top_10perc),
         TOBIAS_highlighted = ifelse(is.na(TOBIAS_highlighted), "FALSE", TOBIAS_highlighted)) %>%
  arrange(desc(ChromVAR), desc(TOBIAS_top_10perc), desc(TOBIAS_highlighted))


writexl::write_xlsx(list("Acute_vs_Control" = sig_acute_res,
                         "Chronic_vs_Control" = sig_chronic_res,
                         "Chronic_vs_Acute" = sig_chronic_acute_res),
                    "~/shared/PRJ_TONGTA/prj_atac_Lyn/res/tobias_res/ChromVAR_TOBIAS_results_20261005.xlsx")




























