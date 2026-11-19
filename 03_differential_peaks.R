library(DESeq2)
library(tidyverse)
# library(edgeR)

# Prepare data and metadata
count_mat <- read_delim("~/shared/PRJ_TONGTA/prj_atac_Lyn/res/featureCounts_res/peak_counts.txt", skip = 1)
count_mat <- count_mat[ , c(1, 7:ncol(count_mat))]

count_mat <-
  count_mat %>%
  column_to_rownames("Geneid") %>%
  as.matrix()

# clean up column names
colnames(count_mat) <-
  str_remove(colnames(count_mat), "\\./shifted_bam/") %>%
  str_remove(., "\\.shifted.sorted.bam")

coldata <-
  colnames(count_mat) %>%
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
  column_to_rownames("samples") %>%
  as.data.frame()

# Removing batch effect
# adj_count_mat <- sva::ComBat_seq(count_mat, batch = coldata$batch, group = NULL)

# Acute vs Control analysis
acute_coldata <- coldata[coldata$condition != "chronic", ]
acute_coldata$condition <- factor(acute_coldata$condition, levels = c("control", "acute"))

dds <- DESeqDataSetFromMatrix(
  countData = count_mat[, c(rownames(acute_coldata))] %>% .[rowSums(.) > 0, ],
  colData = acute_coldata,
  design = ~ batch + condition
)

acute_dds <- DESeq(dds)
DESeq2::plotMA(results(acute_dds))

res_acute_vs_ctrl <- results(acute_dds)
res_acute_vs_ctrl <-
  res_acute_vs_ctrl %>%
  as.data.frame() %>%
  rownames_to_column("peak_id") %>%
  as_tibble() %>%
  filter(!is.na(padj))


# Chronic vs Control analysis
chronic_coldata <- coldata[coldata$condition != "acute", ]
chronic_coldata$condition <- factor(chronic_coldata$condition, levels = c("control", "chronic"))

dds <- DESeqDataSetFromMatrix(
  countData = count_mat[, c(rownames(chronic_coldata))] %>% .[rowSums(.) > 0, ],
  colData = chronic_coldata,
  design = ~ batch + condition
)
#dds <- estimateSizeFactors(dds, type = "poscounts")
chronic_dds <- DESeq(dds)
DESeq2::plotMA(results(chronic_dds))
res_chronic_vs_ctrl <- results(chronic_dds)

res_chronic_vs_ctrl <-
  res_chronic_vs_ctrl %>%
  as.data.frame() %>%
  rownames_to_column("peak_id") %>%
  as_tibble() %>%
  filter(!is.na(padj))

#
# Annotate peaks to genes
#
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Mapped peak ids to genomic coordinates
peak_coords <- read_delim("~/shared/PRJ_TONGTA/prj_atac_Lyn/res/macs3_peaks/consensus_peaks.saf") %>%
  dplyr::rename("peak_id" = "GeneID")

# Acute vs Control
res_acute_vs_ctrl <-
  left_join(res_acute_vs_ctrl, peak_coords)

peak_gr <- GRanges(seqnames = res_acute_vs_ctrl$Chr,
                   ranges = IRanges(start = res_acute_vs_ctrl$Start, end = res_acute_vs_ctrl$End))

anno_peak <- annotatePeak(peak_gr, TxDb = txdb, annoDb = "org.Hs.eg.db", tssRegion = c(-3000, 3000))
anno_peak_df <-
  as.data.frame(anno_peak) %>%
  dplyr::rename("Chr" = "seqnames",
                "Start" = "start",
                "End" = "end")

res_acute_vs_ctrl <-
  res_acute_vs_ctrl %>%
  left_join(., anno_peak_df)


# Chronic vs Control
res_chronic_vs_ctrl <-
  left_join(res_chronic_vs_ctrl, peak_coords)

peak_gr <- GRanges(seqnames = res_chronic_vs_ctrl$Chr,
                   ranges = IRanges(start = res_chronic_vs_ctrl$Start, end = res_chronic_vs_ctrl$End))

anno_peak <- annotatePeak(peak_gr, TxDb = txdb, annoDb = "org.Hs.eg.db", tssRegion = c(-3000, 3000))
anno_peak_df <-
  as.data.frame(anno_peak) %>%
  dplyr::rename("Chr" = "seqnames",
                "Start" = "start",
                "End" = "end")

res_chronic_vs_ctrl <-
  res_chronic_vs_ctrl %>%
  left_join(., anno_peak_df)


# Volcano Plots for sig peaks
library(ggplot2)
library(ggrepel)

plot_volcano <- function(res_df, title, padj_threshold = 0.05, fc_threshold = 1, top_n_label = 30,
                         case_label = "UP", control_label = "Down") {
  
  g <- res_df %>%
    ungroup() %>%
    mutate(sig = case_when(
      padj < padj_threshold & log2FoldChange >= fc_threshold ~ "UP",
      padj < padj_threshold & log2FoldChange <= -fc_threshold ~ "DOWN",
      TRUE ~ "Not Sig"
    )) %>%
    mutate(sig = factor(sig, levels = c("UP", "DOWN", "Not Sig"))) %>%
    arrange(sig == "Not Sig", padj, desc(abs(log2FoldChange))) %>%
    group_by(sig) %>%
    mutate(sig_sequence = 1:n()) %>%
    mutate(label = ifelse(sig != "Not Sig" & sig_sequence <= top_n_label, SYMBOL, NA)) %>%
    ggplot(., aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
    geom_point(size = 1.6, alpha = 0.75) +
    geom_text_repel(aes(label = label),
                    size = 2.5, color = "#0b0b0b", fontface = "italic",
                    max.overlaps = Inf, segment.size = 0.3, segment.color = "#898781",
                    min.segment.length = 0, box.padding = 0.4
    ) +
    scale_color_manual(
      name = "Significance",
      values = c("UP" = "#762a83", "DOWN" = "#1b7837", "Not Sig" = "#c3c2b7"),
      labels = c("UP" = case_label, "DOWN" = control_label, "Not Sig" = "Not Sig")) +
    geom_vline(xintercept = fc_threshold, color = "darkred", linewidth = 0.4, linetype = "dashed") +
    geom_vline(xintercept = -fc_threshold, color = "darkred", linewidth = 0.4, linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_threshold), color = "darkred", linewidth = 0.4, linetype = "dashed") +
    labs(title = title, x = "log2 Fold Change", y = "-log10(adjusted P-value)", color = NULL) +
    theme_classic(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      #axis.line = element_line(color = "#c3c2b7"),
      legend.position = "top"
    )
  
  return(g)
}

acute_plot <- plot_volcano(res_acute_vs_ctrl, "Acute vs Control", 0.05, 0.3, 30,
                           case_label = "Higher in Acute", control_label = "Higher in Control")
chronic_plot <- plot_volcano(res_chronic_vs_ctrl, "Chronic vs Control", 0.05, 0.3, 30,
                             case_label = "Higher in Chronic", control_label = "Higher in Control")

pdf("~/shared/PRJ_TONGTA/prj_atac_Lyn/res/volcano_plots_differential_peaks_20260909.pdf", width = 7, height = 6)
plot(acute_plot)
plot(chronic_plot)
dev.off()

ggsave("~/shared/PRJ_TONGTA/prj_atac_Lyn/res/volcano_plots_differential_peaks_acute_vs_control_20260909.jpg",
       plot = acute_plot, width = 7, height = 6, dpi = 300)
ggsave("~/shared/PRJ_TONGTA/prj_atac_Lyn/res/volcano_plots_differential_peaks_chronic_vs_control_20260909.jpg",
       plot = chronic_plot, width = 7, height = 6, dpi = 300)


deg_res <-
  list(
  "Acute_vs_Control" = res_acute_vs_ctrl,
  "Chronic_vs_Control" = res_chronic_vs_ctrl
)

saveRDS(deg_res, "~/shared/PRJ_TONGTA/prj_atac_Lyn/res/Figures_Tables/deg_res_20260908.rds")



#
# ClusterProfilter - GO enrichment
#
library(clusterProfiler)
library(org.Hs.eg.db)
library(ChIPseeker)

deg_res <- readRDS("~/shared/PRJ_TONGTA/prj_atac_Lyn/res/Figures_Tables/deg_res_20260908.rds")
res_acute_vs_ctrl <- deg_res$Acute_vs_Control
res_chronic_vs_ctrl <- deg_res$Chronic_vs_Control

acute_universe <- res_acute_vs_ctrl$geneId[!is.na(res_acute_vs_ctrl$geneId)]
chronic_universe <- res_chronic_vs_ctrl$geneId[!is.na(res_chronic_vs_ctrl$geneId)]

acute_go_up <- enrichGO(gene = res_acute_vs_ctrl %>%
                          filter(!is.na(geneId)) %>%
                          filter(padj <= 0.05 & log2FoldChange > 0.3) %>% .$geneId,
                        universe = acute_universe,
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENTREZID",
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.1) %>%
  .@result %>%
  filter(pvalue <= 0.05) %>%
  mutate(group = "Acute - UP")

acute_go_down <- enrichGO(gene = res_acute_vs_ctrl %>%
                            filter(!is.na(geneId)) %>%
                            filter(padj <= 0.05 & log2FoldChange < -0.3) %>% .$geneId,
                        universe = acute_universe,
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENTREZID",
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.1) %>%
  .@result %>%
  filter(pvalue <= 0.05) %>%
  mutate(group = "Acute - DOWN")


chronic_go_up <- enrichGO(gene = res_chronic_vs_ctrl %>%
                            filter(!is.na(geneId)) %>%
                            filter(padj <= 0.05 & log2FoldChange > 0) %>% .$geneId,
                        universe = chronic_universe,
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENTREZID",
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.1) %>%
  .@result %>%
  filter(pvalue <= 0.05) %>%
  mutate(group = "Chronic - UP")

chronic_go_down <- enrichGO(gene = res_chronic_vs_ctrl %>%
                              filter(!is.na(geneId)) %>%
                              filter(padj <= 0.05 & log2FoldChange < 0) %>% .$geneId,
                        universe = acute_universe,
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENTREZID",
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.1) %>%
  .@result %>%
  filter(pvalue <= 0.05) %>%
  mutate(group = "Chronic - DOWN")


writexl::write_xlsx(
  list(
    "Acute_vs_Control" = res_acute_vs_ctrl %>% filter(padj <= 0.05 & abs(log2FoldChange) > 0.3),
    "Chronic_vs_Control" = res_chronic_vs_ctrl %>% filter(padj <= 0.05 & abs(log2FoldChange) > 0.3),
    "Acute_BP_UP" = acute_go_up,
    "Acute_BP_DOWN" = acute_go_down,
    "Chronic_BP_UP" = chronic_go_up,
    "Chronic_BP_DOWN" = chronic_go_down
  ), path = "~/shared/PRJ_TONGTA/prj_atac_Lyn/res/Figures_Tables/differential_peaks_results.xlsx"
)


top_BP_dotplot <-
  full_join(acute_go_up, full_join(acute_go_down, full_join(chronic_go_up, chronic_go_down))) %>%
  mutate(GeneRatio = str_split(GeneRatio, "/")) %>%
  mutate(GeneRatio = lapply(GeneRatio, as.integer)) %>%
  rowwise() %>%
  mutate(GeneRatio = GeneRatio[[1]]/GeneRatio[[2]]) %>%
  arrange(p.adjust, desc(FoldEnrichment)) %>%
  mutate(is_sig = p.adjust < 0.05) %>%
  mutate(group = factor(group, levels = c("Acute - UP", "Acute - DOWN", "Chronic - UP", "Chronic - DOWN"))) %>%
  group_by(group) %>%
  slice_head(n = 10) %>%
  arrange(FoldEnrichment) %>%
  mutate(Description = factor(Description, levels = Description)) %>%
  ggplot(., aes(x = FoldEnrichment, y = Description)) +
  geom_point(aes(size = Count, fill = -log10(pvalue), stroke = ifelse(is_sig, 1, 0)),
             shape = 21, color = "#0b0b0b") +
  scale_fill_gradient(low = "#ece7f2", high = "#0570b0", name = "-log10(P.adjusted)") +
  scale_size(range = c(2, 8), name = "Gene Count") +
  labs(x = "Fold Enrichment", y = NULL, title = "GO Biological Process") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "#e1e0d9", linewidth = 0.3),
    axis.text.y = element_text(color = "#0b0b0b"),
    axis.text.x = element_text(color = "#52514e"),
    plot.title  = element_text(color = "#0b0b0b"),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold")
  ) +
  facet_grid(group ~ ., scales = "free_y")


pdf("~/shared/PRJ_TONGTA/prj_atac_Lyn/res/Figures_Tables/top_10_BP_enrichment.pdf", width = 15, height = 8)
plot(top_BP_dotplot)
dev.off()

ggsave("~/shared/PRJ_TONGTA/prj_atac_Lyn/res/Figures_Tables/top_10_BP_enrichment.jpg",
       plot = top_BP_dotplot, width = 15, height = 8, dpi = 300)






