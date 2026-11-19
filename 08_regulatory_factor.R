library(data.table)
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)

bind_dir <- "/home/tongta/shared/PRJ_TONGTA/prj_atac_Lyn/dat/TOBIAS_bindetect"
res <- fread(file.path(bind_dir, "Chronic_vs_Control/bindetect_results.txt"))

res[ ,group := fcase(
  grepl("^(FOS|JUN|ATF3|BATF)", motif_id), "AP1",
  grepl("^(MYOD1|MYF5|MYF6|MYOG)", motif_id), "MRF",
  motif_id == "KLF15", "KLF15",
  default = NA_character_
)]

res %>%
  filter(!is.na(group)) %>%
  dplyr::select(group, name, motif_id, cluster, total_tfbs, Chronic_Control_change, Chronic_Control_pvalue) %>%
  arrange(group, Chronic_Control_change)

read_sites <- function(prefix, cond) {
  f <- file.path(bind_dir, "Chronic_vs_Control", prefix, "beds", paste0(prefix, "_", cond, "_bound.bed"))
  if (!file.exists(f)) return(NULL)
  
  d <- tryCatch(fread(f, select = 1:3, col.names = c("chr", "start", "end")),
                error = function(e) NULL)
  
  if (is.null(d) || nrow(d) == 0) return(NULL)
  d[, start := start + 1L]
  d[, `:=`(motif = prefix, cond = cond)]
}

ap1_prefixes <- res[group == "AP1", output_prefix]
ap1 <- rbindlist(c(lapply(ap1_prefixes, read_sites, cond = "Chronic"),
                   lapply(ap1_prefixes, read_sites, cond = "Control")))

ap1_gr <- GRanges(ap1$chr, IRanges(ap1$start, ap1$end),
                  motif = ap1$motif, cond = ap1$cond)

# The genes of interest
gtf <- import("/home/tongta/shared/ref/GRCh38.p14/gene_annotation/gencode.v50.chr_patch_hapl_scaff.annotation.gtf")
targets <- c("MYOD", "MYF5", "MYF6", "MYOG", "KLF15")
genes <- gtf[gtf$type == "gene" & toupper(gtf$gene_name) %in% targets]


# Table of AP-1 sites within 100 kb of each gene
win <- 100000
ov <- findOverlaps(ap1_gr, genes, maxgap = win)
q <- ap1_gr[queryHits(ov)]
g <- genes[subjectHits(ov)]

tss <- ifelse(as.character(strand(g)) == "-", end(g), start(g))
mid <- (start(q) + end(q))/2
sign <- ifelse(as.character(strand(g)) == "-", -1, 1)

tab <- data.table(
  gene = g$gene_name,
  AP1_motif = q$motif,
  condition = q$cond,
  site = paste0(seqnames(q), ":", start(q), "-", end(q)),
  dist_gene = distance(q, g), # 0 = site lies inside the gene
  dist_TSS = round((mid - tss) * sign) # negative = upstream of the TSS
) %>%
  as.tibble() %>%
  arrange(gene, dist_TSS)

acute_ctrl_tab <- tab
chronic_ctrl_tab <- tab

ctrl_tab <-
  acute_ctrl_tab %>%
  filter(condition == "Control") %>%
  intersect(chronic_ctrl_tab)

writexl::write_xlsx(list(acute_ctrl_tab, chronic_ctrl_tab, ctrl_tab), "~/shared/PRJ_TONGTA/prj_atac_Lyn/res/tobias_res/AP1_sites_near_genes.xlsx",
                    col_names = TRUE, format_headers = TRUE)


