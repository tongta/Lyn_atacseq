
rownames(se) <- paste0(seqnames(se), ":", start(se), "-", end(se))
weighted_count_mat <- assay(se)

library(limma)

# Acute vs control
coldata_acute <- as.data.frame(colData(dev_acute))

design_acute <- model.matrix(~ condition + batch, data = coldata_acute)
fit_acute <- eBayes(lmFit(weighted_count_mat[rownames(coldata_acute), ], design_acute))
res_peaks_acute <- topTable(fit_acute, coef = "conditionacute", number = Inf)
