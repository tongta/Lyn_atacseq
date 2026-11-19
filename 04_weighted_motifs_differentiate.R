library(tidyverse)
library(chromVAR)
library(weightedMotifAccess)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
library(motifmatchr)
library(JASPAR2024)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TFBSTools)
library(RSQLite)


peaks <- import("~/shared/PRJ_TONGTA/prj_atac_Lyn/res/macs3_peaks/consensus_peaks.bed", format = "BED")
keep <- as.character(seqnames(peaks)) %in% standardChromosomes(peaks) &
  as.character(seqnames(peaks)) != "chrM"
peaks <- peaks[keep]
peaks <- sortSeqlevels(peaks)
peaks <- sort(peaks)
seqlevels(peaks) <- seqlevelsInUse(peaks)

# Count summary
peakcount_summary <- read.delim("~/shared/PRJ_TONGTA/prj_atac_Lyn/res/featureCounts_res/peak_counts.txt.summary",
                                row.names = 1)
total_reads <-
  colSums(peakcount_summary) %>%
  as.data.frame() %>%
  rownames_to_column("samples") %>%
  dplyr::rename("depth" = 2) %>%
  mutate(samples = sub(".*\\.shifted_bam\\.", "", samples) %>%
           sub("\\.shifted.*", "", .)) %>%
  mutate(samples = sub("\\.", "-", samples))

# Metadata
coldata <-
  names(bam_files) %>%
  as.tibble() %>%
  dplyr::rename("samples" = 1) %>%
  mutate("condition" = c(rep("acute", 3),
                         rep("control", 3),
                         rep("acute", 3),
                         rep("chronic", 3),
                         rep("control", 3),
                         rep("acute", 3),
                         rep("chronic", 3),
                         rep("control", 3))) %>%
  mutate("batch" = c(rep("batch1", 3), rep("batch2", 18), rep("batch1", 3)) %>%
           factor(., levels = c("batch1", "batch2"))) %>%
  left_join(., total_reads) %>%
  column_to_rownames("samples") %>%
  as.data.frame()

# summarizedExperiment object for chromVAR
bam_files <- list.files("~/shared/PRJ_TONGTA/prj_atac_Lyn/res/shifted_bam",
                        pattern = "\\.shifted\\.sorted\\.bam$", full.names = TRUE)

names(bam_files) <- sub(".*/", "", bam_files) %>% sub("\\.shifted.*", "", .)


se <- getWeightedCounts(
  atac = bam_files,
  ranges = peaks,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  resize = TRUE,
  fragWeight = TRUE,
  peakWeight = TRUE
)
colData(se) <- DataFrame(coldata[colnames(se), ])

# GC-bias correction and filtering
seqlevelsStyle(se) <- "UCSC"
se <- addGCBias(se, genome = BSgenome.Hsapiens.UCSC.hg38)

se <- filterSamples(se, min_depth = 1500, min_in_peaks = 0.15, shiny = FALSE)
se <- filterPeaks(se, non_overlapping = TRUE)

# Match motifs to peaks
opts <- list(species = 9606, collection = "CORE")

jaspar_db <- JASPAR2024()
jaspar_connect <- dbConnect(SQLite(), db(jaspar_db))
motifs <- getMatrixSet(jaspar_connect, opts)

dbDisconnect(jaspar_connect)

motif_ix <- matchMotifs(motifs, se, genome = BSgenome.Hsapiens.UCSC.hg38)
motif_ix <- matchMotifs(motifs, se, genome = BSgenome.Hsapiens.UCSC.hg38, out = "positions")
motif_ix <- matchMotifs(motifs, se, genome = BSgenome.Hsapiens.UCSC.hg38, out = "positions", p.cutoff = 0.01)

# Deviations
dev <- computeDeviations(object = se, annotations = motif_ix)



#
# Comparisons with batch covariate
#

library(limma)

# Acute vs control
dev_acute <- dev[, colData(dev)$condition != "chronic"]
colData(dev_acute)$condition <- factor(colData(dev_acute)$condition, levels = c("control", "acute"))

z_score_acute <- deviationScores(dev_acute)

design_acute <- model.matrix(~ condition + batch, data = colData(dev_acute))
fit_acute <- eBayes(lmFit(z_score_acute, design_acute))
res_acute <- topTable(fit_acute, coef = "conditionacute", number = Inf)

# Chronic vs control
dev_chronic <- dev[, colData(dev)$condition != "acute" & colData(dev)$batch == "batch2"]
colData(dev_chronic)$condition <- factor(colData(dev_chronic)$condition, levels = c("control", "chronic"))

z_score_chronic <- deviationScores(dev_chronic)

design_chronic <- model.matrix(~ condition, data = colData(dev_chronic))
fit_chronic <- eBayes(lmFit(z_score_chronic, design_chronic))
res_chronic <- topTable(fit_chronic, coef = "conditionchronic", number = Inf)

# Chronic vs Acute
dev_chronic_vs_acute <- dev[, colData(dev)$condition != "control" & colData(dev)$batch == "batch2"]
colData(dev_chronic_vs_acute)$condition <- factor(colData(dev_chronic_vs_acute)$condition, levels = c("acute", "chronic"))

z_score_chronic_vs_acute <- deviationScores(dev_chronic_vs_acute)

design_chronic_vs_acute <- model.matrix(~ condition, data = colData(dev_chronic_vs_acute))
fit_chronic_vs_acute <- eBayes(lmFit(z_score_chronic_vs_acute, design_chronic_vs_acute))
res_chronic_vs_acute <- topTable(fit_chronic_vs_acute, coef = "conditionchronic", number = Inf)



plot(res_acute$logFC, -log10(res_acute$P.Value),
     xlab = "logFC", ylab = "-log10(p)",
     col = ifelse(res_acute$adj.P.Val < 0.1, "darkred", "grey"))

plot(res_chronic$logFC, -log10(res_chronic$P.Value),
     xlab = "logFC", ylab = "-log10(p)",
     col = ifelse(res_chronic$adj.P.Val < 0.1, "darkred", "grey"))

plot(res_chronic_vs_acute$logFC, -log10(res_chronic_vs_acute$P.Value),
     xlab = "logFC", ylab = "-log10(p)",
     col = ifelse(res_chronic_vs_acute$adj.P.Val < 0.1, "darkred", "grey"))



# Volcano Plots for sig motifs
library(ggplot2)
library(ggrepel)

motif_map <- data.frame(
  ID = names(motifs),
  symbol = vapply(motifs, TFBSTools::name, character(1)),
  stringsAsFactors = FALSE
)

plot_volcano <- function(res, motif_map, title, fdr_cutoff = 0.05, n_labels = 50,
                         case_label = "UP", ctrl_label = "DOWN") {
  df <- res
  df$ID <- rownames(df)
  df <- merge(df, motif_map, by = "ID")
  
  df$status <- "ns"
  df$status[df$adj.P.Val < fdr_cutoff & df$logFC > 0] <- "up"
  df$status[df$adj.P.Val < fdr_cutoff & df$logFC < 0] <- "down"
  df$status <- factor(df$status, levels = c("down", "ns", "up"))
  
  # horizontal reference line at the actual raw-p boundary implied by this
  # dataset's FDR cutoff (not a fixed p=0.05 line, which wouldn't match the coloring)
  sig_p <- df$P.Value[df$adj.P.Val < fdr_cutoff]
  hline_y <- if (length(sig_p)) -log10(max(sig_p)) else NA
  
  top_labels <- df[df$status != "ns", ]
  top_labels <- top_labels[order(top_labels$adj.P.Val), ]
  top_labels <- head(top_labels, n_labels)
  
  ggplot(df, aes(x = logFC, y = -log10(P.Value), color = status)) +
    geom_point(size = 1.6, alpha = 0.8) +
    geom_text_repel(
      data = top_labels, aes(label = symbol),
      size = 3, color = "#0b0b0b", max.overlaps = 20,
      segment.color = "#898781", segment.size = 0.3
    ) +
    scale_color_manual(
      values = c(down = "#1b7837", ns = "#898781", up = "#762a83"),
      labels = c(down = ctrl_label, ns = "Not significant", up = case_label),
      name = NULL
    ) +
    geom_vline(xintercept = 0, color = "#c3c2b7", linewidth = 0.4) +
    { if (!is.na(hline_y)) geom_hline(yintercept = hline_y, color = "darkred",
                                      linewidth = 0.4, linetype = "dashed") } +
    labs(title = title, x = "log2 fold change", y = expression(-log[10](p))) +
    theme_classic(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      #axis.line = element_line(color = "#c3c2b7"),
      legend.position = "top"
    )
}

acute_plot <- plot_volcano(res_acute, motif_map, "Acute vs Control", fdr_cutoff = 0.1,
                           case_label = "Acute", ctrl_label = "Control")
chronic_plot <- plot_volcano(res_chronic, motif_map, "Chronic vs Control", fdr_cutoff = 0.1,
                             case_label = "Chronic", ctrl_label = "Control")
chronic_vs_acute_plot <- plot_volcano(res_chronic_vs_acute, motif_map, "Chronic vs Acute", fdr_cutoff = 0.1,
                                      case_label = "Chronic", ctrl_label = "Acute")


pdf("~/shared/PRJ_TONGTA/prj_atac_Lyn/res/volcano_plots_differential_motifs_20260916.pdf", width = 7, height = 6)
plot(acute_plot)
plot(chronic_plot)
plot(chronic_vs_acute_plot)
dev.off()

# Map with gene symbol
res_acute <-
  res_acute %>%
  rownames_to_column("ID") %>%
  left_join(motif_map)

res_chronic <-
  res_chronic %>%
  rownames_to_column("ID") %>%
  left_join(motif_map)

res_chronic_vs_acute <-
  res_chronic_vs_acute %>%
  rownames_to_column("ID") %>%
  left_join(motif_map)



writexl::write_xlsx(list(res_acute, res_chronic, res_chronic_vs_acute),
                    "~/shared/PRJ_TONGTA/prj_atac_Lyn/res/Figures_Tables/differential_motif_results_20260921.xlsx",
                    format_headers = TRUE)








