library(Gviz)
library(rtracklayer)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

gene_symbol <- c("TP53", "HOXB9", "HOXC9", "HOXC10", "HOXD10", "HOXD11",
                 "IRF2", "IRF3", "IRF4", "IRF5", "IRF6", "IRF7", "IRF8", "IRF9")
gene_symbol <-
  res_acute %>% filter(adj.P.Val <= 0.1) %>% dplyr::select(ID, symbol) %>%
  full_join(., res_chronic %>% filter(adj.P.Val <= 0.1) %>% dplyr::select(ID, symbol)) %>%
  full_join(., res_chronic_vs_acute %>% filter(adj.P.Val <= 0.1) %>% dplyr::select(ID, symbol)) %>%
  .$symbol %>% sort() %>% unique()
gene_symbol <- str_split(gene_symbol, "::") %>% unlist()

entrez_id <- mapIds(org.Hs.eg.db, keys = gene_symbol, column = "ENTREZID", keytype = "SYMBOL")

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
all_genes <- genes(txdb, single.strand.genes.only = FALSE)
gene_grl <- all_genes[names(all_genes) %in% entrez_id]

gene_span <- unlist(range(gene_grl))
gene_span$entrez_id <- names(gene_span)
gene_span$symbol <- names(entrez_id)[match(gene_span$entrez_id, entrez_id)]

plot_locus <- function(symbol, gr, control_bw_path, acute_bw_path, chronic_bw_path,
                       se, motif_pos, motif_name) {
  chr <- as.character(seqnames(gr))
  gene_start <- start(gr)
  gene_end <- end(gr)
  
  # Searching for motifs in the gene area
  flank = 50000
  start <- gene_start - flank
  end <- gene_end + flank
  region <- GRanges(chr, IRanges(start, end))
  
  locus_motifs <- subsetByOverlaps(motif_pos[[motif_name]], region)
  
  if (length(locus_motifs) == 0) {
    flank = 100000
    start <- gene_start - flank
    end <- gene_end + flank
    region <- GRanges(chr, IRanges(start, end))
    locus_motifs <- subsetByOverlaps(motif_pos[[motif_name]], region)
    if (length(locus_motifs) == 0) {
      return(invisible(NULL))
    }
  }
  
  motif_start <- min(start(locus_motifs))
  motif_end <- max(end(locus_motifs))
  
  if (motif_start < gene_start & motif_end > gene_end) {
    start <- motif_start - 5000
    end <- motif_end + 5000
  } else if (motif_start < gene_start & motif_end < gene_end) {
    start <- motif_start - 5000
    end <- gene_end
  } else if (motif_start > gene_start & motif_end > gene_end) {
    start <- gene_start
    end <- motif_end + 5000
  } else {
    start <- gene_start - 5000
    end <- gene_end + 5000
  }
  
  region <- GRanges(chr, IRanges(start, end))
  
  locus_motifs <- subsetByOverlaps(motif_pos[[motif_name]], region)
  motif_track <- AnnotationTrack(
    locus_motifs, genome = "hg38", chromosome = chr,
    name = motif_name, fill = "#377eb8", col = "#377eb8", stacking = "dense"
  )
  
  
  itrack <- IdeogramTrack(genome = "hg38", chromosome = chr)
  gtrack <- GenomeAxisTrack()
  
  grtrack <- GeneRegionTrack(
    txdb, chromosome = chr, start = start, end = end,
    name = symbol, geneSymbol = TRUE, showId = FALSE
  )
  gene(grtrack) <- mapIds(org.Hs.eg.db, keys = gene(grtrack), column = "SYMBOL", keytype = "ENTREZID")
  grtrack <- grtrack[gene(grtrack) == symbol]
  
  control_bw <- import.bw(control_bw_path, which = region)
  control_bw <- keepSeqlevels(control_bw, chr, pruning.mode = "coarse")
  
  acute_bw <- import.bw(acute_bw_path, which = region)
  acute_bw <- keepSeqlevels(acute_bw, chr, pruning.mode = "coarse")
  
  chronic_bw <- import.bw(chronic_bw_path, which = region)
  chronic_bw <- keepSeqlevels(chronic_bw, chr, pruning.mode = "coarse")
  
  y_max = max(c(score(control_bw), score(acute_bw), score(chronic_bw)))
  control_track <- DataTrack(control_bw, genome = "hg38", type = "histogram",
                             name = "Control (mean RPGC)", fill = "#4daf4a", col.histogram = "#4daf4a",
                             ylim = c(0, y_max))
  acute_track <- DataTrack(acute_bw, genome = "hg38", type = "histogram",
                           name = "Acute (mean RPGC)", fill = "#984ea3", col.histogram = "#984ea3",
                           ylim = c(0, y_max))
  chronic_track <- DataTrack(chronic_bw, genome = "hg38", type = "histogram",
                             name = "Chronic (mean RPGC)", fill = "#ff7f00", col.histogram = "#ff7f00",
                             ylim = c(0, y_max))
  
  locus_peaks <- subsetByOverlaps(rowRanges(se), region)
  peak_track <- AnnotationTrack(
    locus_peaks, genome = "hg38", chromosome = chr, name = "Peaks",
    fill = "#e41a1c", col = "#e41a1c", stacking = "dense"
  )
  
  plotTracks(
    list(itrack, gtrack, grtrack, control_track, acute_track, chronic_track, peak_track, motif_track),
    from = start, to = end, chromosome = chr,
    sizes = c(1, 1, 4, 2, 2, 2, 1, 1),
    main = symbol
  )
}


plot_locus(
  symbol = "TP53",
  gr = gene_span[gene_span$symbol == "TP53"],
  control_bw_path = "~/shared/PRJ_TONGTA/prj_atac_Lyn/res/bigwig_rpgc/control_chronic_avg.bw",
  acute_bw_path = "~/shared/PRJ_TONGTA/prj_atac_Lyn/res/bigwig_rpgc/acute_batch2_avg.bw",
  chronic_bw_path = "~/shared/PRJ_TONGTA/prj_atac_Lyn/res/bigwig_rpgc/chronic_avg.bw",
  se = se, motif_pos = motif_ix, motif_name = "MA0106.3"
)


gene_span <- gene_span[gene_span@seqnames %in% standardChromosomes(gene_span), ]
motif_list <- 
  res_acute %>% filter(symbol %in% gene_symbol) %>%
  full_join(., res_chronic %>% filter(symbol %in% gene_symbol)) %>%
  full_join(., res_chronic_vs_acute %>% filter(symbol %in% gene_symbol)) %>%
  dplyr::select(ID, symbol) %>%
  unique() %>%
  arrange(symbol) %>%
  dplyr::rename("SYMBOL" = "symbol")

plot_list_table <-
  gene_symbol %>%
  as_tibble() %>%
  dplyr::rename("SYMBOL" = 1) %>%
  left_join(., motif_list) %>%
  rowwise() %>%
  mutate(gr = list(gene_span[gene_span$symbol == SYMBOL])) %>%
  filter(!is.na(ID))


pdf("~/shared/PRJ_TONGTA/prj_atac_Lyn/res/Figures_Tables/locus_motif_plots_rpgc_norm.pdf",
    width = 12, height = 10)
for (i in 1:nrow(plot_list_table)) {
  print(plot_list_table$SYMBOL[[i]])
  plot_locus(
    symbol = plot_list_table$SYMBOL[[i]],
    gr = plot_list_table$gr[[i]],
    control_bw_path = "~/shared/PRJ_TONGTA/prj_atac_Lyn/res/bigwig_rpgc/control_chronic_avg.bw",
    acute_bw_path = "~/shared/PRJ_TONGTA/prj_atac_Lyn/res/bigwig_rpgc/acute_batch2_avg.bw",
    chronic_bw_path = "~/shared/PRJ_TONGTA/prj_atac_Lyn/res/bigwig_rpgc/chronic_avg.bw",
    se = se,
    motif_pos = motif_ix,
    motif_name = plot_list_table$ID[[i]]
  )
}
dev.off()












