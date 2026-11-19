library(tidyverse)
library(chromVAR)
library(motifmatchr)
library(JASPAR2024)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
library(TFBSTools)
library(RSQLite)

# Peak coordinates
peaks <- import("~/shared/PRJ_TONGTA/prj_atac_Lyn/res/macs3_peaks/consensus_peaks.bed", format = "BED")
keep <- as.character(seqnames(peaks)) %in% standardChromosomes(peaks) &
  as.character(seqnames(peaks)) != "chrM"
peaks <- peaks[keep]


# Peak counts
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

count_mat <- count_mat[keep, ]

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

# Building the summarizedExperiment for chromVAR
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
  left_join(., total_reads) %>%
  column_to_rownames("samples") %>%
  as.data.frame()

se <- SummarizedExperiment(
  assays = list(counts = count_mat),
  rowRanges = peaks,
  colData = coldata
)


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

# Deviations
dev <- computeDeviations(object = se, annotations = motif_ix)

#
# Acute vs Control test
#
dev_acute <- dev[, colData(dev)$condition != "chronic"]
colData(dev_acute)$condition <- factor(colData(dev_acute)$condition, levels = c("control", "acute"))

var_acute <- computeVariability(dev_acute)
plotVariability(var_acute, use_plotly = FALSE)

# Differential analysis
diff_acute_dev <- differentialVariability(dev_acute, groups = "condition")
diff_acute_dev <- diff_acute_dev[order(diff_acute_dev$p_value), ]



#
# Chronic vs Control test
#
dev_chronic <- dev[, colData(dev)$condition != "acute" & colData(dev)$batch == "batch2"]
colData(dev_chronic)$condition <- factor(colData(dev_chronic)$condition, levels = c("control", "chronic"))

var_chronic <- computeVariability(dev_chronic)
plotVariability(var_chronic, use_plotly = FALSE)

# Differential analysis
diff_chronic_dev <- differentialVariability(dev_chronic, groups = "condition")
diff_chronic_dev <- diff_chronic_dev[order(diff_chronic_dev$p_value), ]




#
# Chronic vs Acute test
#
dev_chronic_vs_acute <- dev[, colData(dev)$condition != "control" & colData(dev)$batch == "batch2"]
colData(dev_chronic_vs_acute)$condition <- factor(colData(dev_chronic_vs_acute)$condition, levels = c("acute", "chronic"))

var_chronic_vs_acute <- computeVariability(dev_chronic_vs_acute)
plotVariability(var_chronic_vs_acute, use_plotly = FALSE)

# Differential analysis
diff_chronic_acute_dev <- differentialVariability(dev_chronic_vs_acute, groups = "condition")
diff_chronic_acute_dev <- diff_chronic_acute_dev[order(diff_chronic_acute_dev$p_value), ]



#
# With batch covariate
#

library(limma)

# Acute vs control
z_score_acute <- deviationScores(dev_acute)

design_acute <- model.matrix(~ condition + batch, data = colData(dev_acute))
fit_acute <- eBayes(lmFit(z_score_acute, design_acute))
res_acute <- topTable(fit_acute, coef = "conditionacute", number = Inf)

# Chronic vs control
z_score_chronic <- deviationScores(dev_chronic)

design_chronic <- model.matrix(~ condition, data = colData(dev_chronic))
fit_chronic <- eBayes(lmFit(z_score_chronic, design_chronic))
res_chronic <- topTable(fit_chronic, coef = "conditionchronic", number = Inf)

# Chronic vs Acute
z_score_chronic_vs_acute <- deviationScores(dev_chronic_vs_acute)

design_chronic_vs_acute <- model.matrix(~ condition, data = colData(dev_chronic_vs_acute))
fit_chronic_vs_acute <- eBayes(lmFit(z_score_chronic_vs_acute, design_chronic_vs_acute))
res_chronic_vs_acute <- topTable(fit_chronic_vs_acute, coef = "conditionchronic", number = Inf)

