library(motifmatchr)
library(TFBSTools)
library(universalmotif)
library(RSQLite)

jaspar_db <- JASPAR2024()
sq24 <- dbConnect(SQLite(), db(jaspar_db))

opts <- list(species = 9606, collection = "CORE")
pfm_list <- getMatrixSet(sq24, opts)

motifs <- convert_motifs(pfm_list, class = "universalmotif-universalmotif")
write_meme(motifs, "~/shared/PRJ_TONGTA/prj_atac_Lyn/dat/jaspar_human_core.meme")


dbDisconnect(sq24)
