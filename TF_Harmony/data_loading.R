# data_loading.R
# Loads preprocessed RDS files created by preprocess_data.R.
# If the RDS files don't exist, falls back to reading raw TSVs and processing inline.
# Run `Rscript preprocess_data.R` once to generate the RDS files.

if (file.exists("Data/dt.rds")) {
  ## --- Fast path: load preprocessed RDS files ---

  ## Color scheme for network graphs overlaying activator/repressor data with TF activity
  testcolors <- c("red", "gray", "blue", "gray", "darkred", "darkblue")
  names(testcolors) <- c("Activator", "Minimally Active", "Repressor", "Unknown", "Upregulated", "Downregulated")

  ## Pre-calculated harmony data with family annotations, exp17 removed, names cleaned
  dt <- readRDS("Data/dt.rds")
  ## TF family lookup table
  tfswithfamilies <- readRDS("Data/tfswithfamilies.rds")
  ## Pairwise harmony with fewer columns (used in Pairwise Analyses)
  newdt <- readRDS("Data/newdt.rds")

  ## DEGs in narrow format, with color column, HSF.A4A fix, and exp17 removed
  narpv <- readRDS("Data/narpv.rds")
  ## Unique gene IDs from DEGs, used for Target Regulation selectize
  allgeneids <- readRDS("Data/allgeneids.rds")
  ## Ordered unique TF names from DEGs
  idoptions <- readRDS("Data/idoptions.rds")
  ## TF ID to TF name mapping, derived from narpv
  ids <- readRDS("Data/ids.rds")

  ## Known transcription effector domains
  allteds <- readRDS("Data/allteds.rds")
  ## Pre-built vertices for network graphs with shapes assigned
  vertices <- readRDS("Data/vertices.rds")

  ## PWMs for motif sorting and PFMs for motifStack
  pwms <- readRDS("Data/pwms_processed.rds")
  pfms <- readRDS("Data/pfms.rds")

  ## Cell type expression from Benfey lab
  cte <- readRDS("Data/cte.rds")
  ## Just-in-time datasets for roots and shoots
  jitr <- readRDS("Data/jitr.rds")
  jits <- readRDS("Data/jits.rds")

  ## Pre-computed phylogenetic dendrogram from MSA (pruned per-render in server)
  phylo_dend <- readRDS("Data/phylo_dend.rds")

} else {
  ## --- Fallback: process from raw files (slow, first run only) ---
  message("RDS files not found in Data/. Loading from raw TSVs (slow). Run preprocess_data.R to fix this.")

  ## Setting up color scheme that is used in a couple of the network graphs
  testcolors <- c("red", "gray", "blue", "gray", "darkred", "darkblue")
  names(testcolors) <- c("Activator", "Minimally Active", "Repressor", "Unknown", "Upregulated", "Downregulated")

  ## This msa is used for the protein sequence dendrogram
  htmsa <- Biostrings::readAAMultipleAlignment("Data/msa.fa")

  ## Pre-compute the phylogenetic tree from the MSA
  phylo_alignment <- msa::msaConvert(htmsa, type = "seqinr::alignment")
  phylo_distance <- seqinr::dist.alignment(phylo_alignment)
  phylo_tree <- ape::bionj(phylo_distance)
  phylo_dend <- as.dendrogram.phylo(phylo_tree)

  ## Known transcription effector domains
  allteds <- fread("Data/allteds.tsv")

  ## Reading in the pre-calculated harmony
  dt <- fread("Data/harmonytable_nobatch_nebs.tsv")

  setnames(dt,
    c("TF", "Intersect_Concordant", "Intersect_Discordant",
      "PValue_Concordant", "PValue_Discordant",
      "Correlation_Concordant", "Correlation_Discordant",
      "Harmony_Concordant", "Harmony_Discordant"),
    c("TF1", "Concordant_Intersect", "Discordant_Intersect",
      "Concordant_PValue", "Discordant_PValue",
      "Concordant_Correlation", "Discordant_Correlation",
      "Concordant_Harmony", "Discordant_Harmony"))

  dt <- dt[!(TF1 == TF2)]
  dt[, TF1 := sub("-.*", "", TF1)]
  dt[, TF2 := sub("-.*", "", TF2)]

  tfswithfamilies <- fread("Data/tfidswithfamilies.tsv")

  setkey(dt, TF1)
  setkey(tfswithfamilies, Name)
  dt <- dt[tfswithfamilies[, .(Name, TF1_Family = Family)]]
  setkey(dt, TF2)
  dt <- dt[tfswithfamilies[, .(Name, TF2_Family = Family)]]

  pwms <- readRDS("Data/pwms.RDS")
  pfms <- convert_motifs(pwms, class = "motifStack-pfm")

  narpv <- fread("Data/newnarpv_nebs_nobatch.tsv")
  narpv[, Color := ifelse(sign(log2FoldChange) < 0, "Blue", "Red")]
  narpv[TF == "HSF.A4A", TF_ID := "AT4G18880"]

  exp17nogo <- narpv[EXP == "1-17", unique(TF)]
  narpv <- narpv[EXP != "1-17"]
  dt <- dt[!(TF1 %in% exp17nogo)]
  dt <- dt[!(TF2 %in% exp17nogo)]

  allgeneids <- unique(narpv$rn)

  gis <- narpv[, .(unique(rn))]
  vertices <- merge(gis, allteds, by.x = "V1", by.y = "Locus", all.x = T)
  vertices <- merge(vertices, unique(narpv[, .(TF_ID, TF)]), by.x = "V1", by.y = "TF_ID", all.x = T)

  cte <- fread("Data/cellTypeExpression.tsv", key = "Gene ID")

  jitr <- fread("Data/JITGenes_root.csv", skip = 1,
    select = c("Gene", "FDR adjusted p-value", "First Response (Just-in-time bin)"),
    col.names = c("GeneID", "pvalue", "JIT"))
  setkey(jitr, GeneID)

  jits <- fread("Data/JITGenes_shoot.csv", skip = 1,
    select = c("AtID", "FDR adjusted p-value", "First Response (Just-in-time bin)"),
    col.names = c("GeneID", "pvalue", "JIT"))
  setkey(jits, GeneID)

  idoptions <- narpv[order(TF), unique(TF)]
  ids <- unique(narpv[, .(Ids = TF_ID, TF)])

  newdt <- dt[, .(
    TF1 = sub("AT.*?_", "", TF1),
    TF2 = sub("AT.*?_", "", TF2),
    Concordant = Concordant_Harmony,
    Discordant = Discordant_Harmony
  )][order(-Concordant)]

  fulltfids <- fread("Data/all_ath_tf_ids.tsv")

  vertices <- merge(vertices, fulltfids, by.x = "V1", by.y = "Gene_ID", all.x = T)
  vertices[!(is.na(Family)), shape := "triangle"]
  vertices[(is.na(Family)), shape := "circle"]
  vertices[!(is.na(TF)), shape := "triangle"]
  vertices <- unique(vertices)
}
