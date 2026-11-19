#
# Heatmap
#

library(tidyverse)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggthemes)
library(grid)

# JASPAR motif family
lines <- readLines("~/shared/PRJ_TONGTA/prj_atac_Lyn/dat/JASPAR2026_CORE_non-redundant_pfms_transfac.txt")

record_id <- cumsum(c(FALSE, head(lines, -1) == "//"))
records <- split(lines, record_id)

parse_record <- function(rec) {
  get_field <- function(tag) {
    hit <- grep(paste0("^", tag, "\\s"), rec, value = TRUE)
    if (length(hit) == 0) return(NA_character_)
    trimws(sub(paste0("^", tag, "\\s+"), "", hit[1]))
  }
  get_cc <- function(key) {
    hit <- grep(paste0("^CC\\s+", key, ":"), rec, value = TRUE)
    if (length(hit) == 0) return(NA_character_)
    trimws(sub(paste0("^CC\\s+", key, ":"), "", hit[1]))
  }
  
  data.frame(
    motif_id = get_field("AC"),
    gene_name = get_field("ID"),
    family = get_cc("tf_family"),
    tf_class = get_cc("tf_class"),
    stringsAsFactors = FALSE
  )
}

jaspar_annotation <- do.call(rbind, lapply(records, parse_record))
jaspar_annotation <- jaspar_annotation[!is.na(jaspar_annotation), ]

# TOBIAS results
acute_res_df <- readxl::read_xlsx("~/shared/PRJ_TONGTA/prj_atac_Lyn/dat/TOBIAS_bindetect/Acute_vs_Control/bindetect_results.xlsx")
chronic_res_df <- readxl::read_xlsx("~/shared/PRJ_TONGTA/prj_atac_Lyn/dat/TOBIAS_bindetect/Chronic_vs_Control/bindetect_results.xlsx")
#chronic_acute_df <- readxl::read_xlsx("~/shared/PRJ_TONGTA/prj_atac_Lyn/dat/TOBIAS_bindetect/Chronic_vs_Acute/bindetect_results.xlsx")


acute_res_df <-
  acute_res_df %>%
  mutate(Acute_Control_padj = p.adjust(Acute_Control_pvalue, method = "BH")) %>%
  dplyr::select(name, motif_id, Acute_Control_change, Acute_Control_padj)

chronic_res_df <-
  chronic_res_df %>%
  mutate(Chronic_Control_padj = p.adjust(Chronic_Control_pvalue, method = "BH")) %>%
  dplyr::select(name, motif_id, Chronic_Control_change, Chronic_Control_padj)

# 
# chronic_acute_df <-
#   chronic_acute_df %>%
#   mutate(Chronic_Acute_padj = p.adjust(Chronic_Acute_pvalue, method = "BH"))

# Significant motifs
acute_motifs <-
  acute_res_df %>%
  filter(Acute_Control_padj < 0.05 & abs(Acute_Control_change) >= 0.1) %>%
  arrange(desc(abs(Acute_Control_change))) %>%
  slice_head(n = 30) %>%
  .$name

chronic_motifs <-
  chronic_res_df %>%
  filter(Chronic_Control_padj < 0.05 & abs(Chronic_Control_change) >= 0.1) %>%
  arrange(desc(abs(Chronic_Control_change))) %>%
  slice_head(n = 30) %>%
  .$name

combined <- full_join(acute_res_df, chronic_res_df) %>%
  mutate(
    label = paste0(motif_id, " (", name, ")")
  ) %>%
  filter(name %in% c(acute_motifs, chronic_motifs))


top_motifs <- combined %>%
  mutate(Acute_minus_Chronic = Acute_Control_change - Chronic_Control_change) %>%
  dplyr::rename("Acute" = "Acute_Control_change",
                "Chronic" = "Chronic_Control_change")

# Matrices for plotting
top_motifs <-
  top_motifs %>%
  dplyr::rename("gene_name" = "motif_id",
                "motif_id" = "name") %>%
  left_join(jaspar_annotation %>%
              dplyr::select(motif_id, family), by = "motif_id") %>%
  mutate(
    direction = if_else(Acute >= 0, "Up", "Down"),
    direction = factor(direction, levels = c("Up", "Down"))
  ) %>%
  arrange(family) %>%
  group_by(family) %>%
  mutate(mean_diff = mean(Acute_minus_Chronic)) %>%
  ungroup() %>%
  arrange(desc(mean_diff))

mat_diff <- top_motifs %>%
  dplyr::select(label, Acute_minus_Chronic) %>%
  column_to_rownames("label") %>%
  as.matrix()

mat_main <- top_motifs %>%
  dplyr::select(label, Acute, Chronic) %>%
  column_to_rownames("label") %>%
  as.matrix()

lim_main <- max(abs(mat_main), na.rm = TRUE)
lim_diff <- max(abs(mat_diff), na.rm = TRUE)


# Heatmap setting and plotting
top_motifs$family <- as.character(top_motifs$family)
top_motifs$family[is.na(top_motifs$family)] <- "Other"

families <- unique(top_motifs$family)
tableau20 <- ggthemes::tableau_color_pal("Tableau 20")(20)
fam_colors <- setNames(tableau20[1:length(families)], families)

row_anno <- rowAnnotation(
  Family = top_motifs$family,
  col = list(Family = fam_colors),
  show_legend = FALSE,
  show_annotation_name = FALSE,
  annotation_name_side = "top",
  annotation_name_rot = 0,
  annotation_name_gp = gpar(fontsize = 10)
)

slice_order <-
  top_motifs %>%
  distinct(mean_diff, family) %>%
  arrange(mean_diff, family)

cluster_anno <- rowAnnotation(
  cluster_label = anno_block(
    labels = slice_order$family,
    labels_gp = gpar(fontsize = 9),
    labels_rot = 0,
    labels_just = "right",
    labels_offset = unit(1, "npc"),
    gp = gpar(col = NA, fill = NA)
  ),
  width = unit(7.5, "cm"),
  show_annotation_name = FALSE
)

row_split_df <- data.frame(
  Mean_Diff = top_motifs$mean_diff
)

acute_lab <-
  top_motifs %>%
  filter(motif_id %in% acute_motifs) %>%
  .$label
chronic_lab <-
  top_motifs %>%
  filter(motif_id %in% chronic_motifs) %>%
  .$label

ht_main <-
  Heatmap(
    mat_main,
    name = "change",
    cluster_columns = FALSE,
    row_split = row_split_df,
    row_title = NULL,
    cluster_rows = TRUE,
    cluster_row_slices = FALSE,
    col = colorRamp2(c(-lim_main, 0, lim_main), c("#4393c3", "#f7f7f7", "#d6604d")),
    show_row_names = FALSE,
    column_split = c("Acute/Control", "Chronic/Control"),
    column_labels = c("", ""),
    cell_fun = function(j, i, x, y, width, height, fill) {
      val <- mat_main[i, j]
      if(!is.na(val)) {
        lab <- sprintf("%.3f", val)
        if (val > 0) lab <- paste0("+", lab)
        
        motif <- rownames(mat_main)[i]
        col_name <- colnames(mat_main)[j]
        
        mark <- (col_name == "Acute" && motif %in% acute_lab) ||
                (col_name == "Chronic" && motif %in% chronic_lab)
        if (mark) lab <- paste0(lab, "*")
        
        grid.text(lab, x, y, gp = gpar(fontsize = 7))
      }
    }
  )


ht_diff <-
  Heatmap(
    mat_diff,
    name = "Acute - Chronic",
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    col = colorRamp2(c(-lim_diff, 0, lim_diff), c("#542788", "#f7f7f7", "#b35806")),
    show_row_names = TRUE,
    row_names_side = "right",
    row_names_gp = gpar(fontsize = 9, fontface = 'italic'),
    width = unit(0.5, "cm"),
    column_title = "",
    column_labels = ""
  )

cluster_anno + row_anno + ht_main + ht_diff

pdf("~/shared/PRJ_TONGTA/prj_atac_Lyn/res/Figures_Tables/TOBIAS_heatmap_tf_change_top30_motifs.pdf",
    width = 10, height = 10)
draw(cluster_anno + row_anno + ht_main + ht_diff,
     heatmap_legend_side = "right",
     annotation_legend_side = "right")
dev.off()

jpeg("~/shared/PRJ_TONGTA/prj_atac_Lyn/res/Figures_Tables/TOBIAS_heatmap_tf_change_top30_motifs.jpg",
     width = 10, height = 10, units = "in", res = 300)
draw(cluster_anno + row_anno + ht_main + ht_diff,
     heatmap_legend_side = "right",
     annotation_legend_side = "right")
dev.off()


# ChromVAR results
acute_res_df <- 
  readxl::read_xlsx("~/shared/PRJ_TONGTA/prj_atac_Lyn/res/Figures_Tables/differential_motif_results_20260921.xlsx",
                    sheet = 1) %>%
  dplyr::select(ID, symbol, logFC, adj.P.Val) %>%
  dplyr::rename("motif_id" = "ID",
                "gene_name" = "symbol",
                "Acute_logFC" = "logFC",
                "Acute_padj" = "adj.P.Val")

chronic_res_df <- 
  readxl::read_xlsx("~/shared/PRJ_TONGTA/prj_atac_Lyn/res/Figures_Tables/differential_motif_results_20260921.xlsx",
                    sheet = 2) %>%
  dplyr::select(ID, symbol, logFC, adj.P.Val) %>%
  dplyr::rename("motif_id" = "ID",
                "gene_name" = "symbol",
                "Chronic_logFC" = "logFC",
                "Chronic_padj" = "adj.P.Val")

# Significant motifs
acute_motifs <-
  acute_res_df %>%
  filter(Acute_padj < 0.1 & abs(Acute_logFC) >= 0.5) %>%
  .$motif_id

chronic_motifs <-
  chronic_res_df %>%
  filter(Chronic_padj < 0.1 & abs(Chronic_logFC) >= 0.5) %>%
  .$motif_id

combined <- full_join(acute_res_df, chronic_res_df) %>%
  mutate(
    label = paste0(gene_name, " (", motif_id, ")")
  ) %>%
  filter(motif_id %in% c(acute_motifs, chronic_motifs))


top_motifs <- combined %>%
  mutate(Acute_minus_Chronic = Acute_logFC - Chronic_logFC) %>%
  dplyr::rename("Acute" = "Acute_logFC",
                "Chronic" = "Chronic_logFC")

# Matrices for plotting
top_motifs <-
  top_motifs %>%
  left_join(jaspar_annotation %>%
              dplyr::select(motif_id, family), by = "motif_id") %>%
  mutate(
    direction = if_else(Acute >= 0, "Up", "Down"),
    direction = factor(direction, levels = c("Up", "Down"))
  ) %>%
  arrange(family) %>%
  group_by(family) %>%
  mutate(mean_diff = mean(Acute_minus_Chronic)) %>%
  ungroup() %>%
  arrange(mean_diff)

mat_main <- top_motifs %>%
  dplyr::select(label, Acute, Chronic) %>%
  column_to_rownames("label") %>%
  as.matrix()

mat_diff <- top_motifs %>%
  dplyr::select(label, Acute_minus_Chronic) %>%
  column_to_rownames("label") %>%
  as.matrix()

lim_main <- max(abs(mat_main), na.rm = TRUE)
lim_diff <- max(abs(mat_diff), na.rm = TRUE)


# Heatmap setting and plotting
top_motifs$family <- as.character(top_motifs$family)
top_motifs$family[is.na(top_motifs$family)] <- "Other"

families <- unique(top_motifs$family)
tableau20 <- ggthemes::tableau_color_pal("Tableau 20")(20)
fam_colors <- setNames(tableau20[1:length(families)], families)

row_anno <- rowAnnotation(
  Family = top_motifs$family,
  col = list(Family = fam_colors),
  show_legend = FALSE,
  show_annotation_name = FALSE,
  annotation_name_side = "top",
  annotation_name_rot = 0,
  annotation_name_gp = gpar(fontsize = 10)
)

slice_order <-
  top_motifs %>%
  distinct(mean_diff, family) %>%
  arrange(mean_diff, family)

cluster_anno <- rowAnnotation(
  cluster_label = anno_block(
    labels = slice_order$family,
    labels_gp = gpar(fontsize = 9),
    labels_rot = 0,
    labels_just = "right",
    labels_offset = unit(1, "npc"),
    gp = gpar(col = NA, fill = NA)
  ),
  width = unit(7.5, "cm"),
  show_annotation_name = FALSE
)

row_split_df <- data.frame(
  Mean_Diff = top_motifs$mean_diff
)

ht_main <-
  Heatmap(
    mat_main,
    name = "logFC",
    cluster_columns = FALSE,
    row_split = row_split_df,
    row_title = NULL,
    cluster_rows = TRUE,
    cluster_row_slices = FALSE,
    col = colorRamp2(c(-lim_main, 0, lim_main), c("#4393c3", "#f7f7f7", "#d6604d")),
    show_row_names = FALSE,
    column_split = c("Acute/Control", "Chronic/Control"),
    column_labels = c("", ""),
    cell_fun = function(j, i, x, y, width, height, fill) {
      val <- mat_main[i, j]
      if(!is.na(val)) {
        if (val > 0) {
          val = paste0("+", as.character(sprintf("%.3f", val)))
        } else {
          val = as.character(sprintf("%.3f", val))
        }
        grid.text(val, x, y, gp = gpar(fontsize = 7))
      }
    }
  )


ht_diff <-
  Heatmap(
    mat_diff,
    name = "Acute - Chronic",
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    col = colorRamp2(c(-lim_diff, 0, lim_diff), c("#542788", "#f7f7f7", "#b35806")),
    show_row_names = TRUE,
    row_names_side = "right",
    row_names_gp = gpar(fontsize = 9, fontface = 'italic'),
    width = unit(0.5, "cm"),
    column_title = "",
    column_labels = ""
  )

cluster_anno + row_anno + ht_main + ht_diff

pdf("~/shared/PRJ_TONGTA/prj_atac_Lyn/res/Figures_Tables/ChromVAR_heatmap_tf_change.pdf",
    width = 10, height = 10)
draw(cluster_anno + row_anno + ht_main + ht_diff,
     heatmap_legend_side = "right",
     annotation_legend_side = "right")
dev.off()

jpeg("~/shared/PRJ_TONGTA/prj_atac_Lyn/res/Figures_Tables/ChromVAR_heatmap_tf_change.jpg",
     width = 10, height = 10, units = "in", res = 300)
draw(cluster_anno + row_anno + ht_main + ht_diff,
     heatmap_legend_side = "right",
     annotation_legend_side = "right")
dev.off()












