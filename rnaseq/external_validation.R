##################################
# Anna Schwager
# CNRS UMR9018, Institut Gustave Roussy
# 2026
##################################

library(GEOquery)
library(dplyr)
library(readxl)
library(tibble)
library(ggplot2)
library(ggpubr)
library(DESeq2)
library(GSEABase)
library(singscore)
library(hgu133plus2.db)

################### Paths ################################

setwd("/Users/annaschwager/Documents/projects/MCL/revision_analysis/public_validation/")
output_dir <- paste0(getwd(),"/output/")

################### Load chr19 MCL signature ###########################

mcl_signature <- read_excel(
  "/Users/annaschwager/Documents/projects/MCL/DAG_MCL_deseq2/results/chr19_upregulated_genes_with_distance_and_function.xlsx",
  sheet = 1
) %>%
  pull(gene) %>%
  unique()


################### Load our RNA-seq data ######################################
my_samples <- read.csv2(
  "/Users/annaschwager/Documents/projects/MCL/revision_analysis/rnaseq/input/counts.csv")

counts_ours <- my_samples %>%
  column_to_rownames("X") %>%
  as.matrix()

storage.mode(counts_ours) <- "integer"

keep_samples_ours <- grep(
  "^MCL_|^naive_blood_",
  colnames(counts_ours),
  value = TRUE
)

counts_ours <- counts_ours[
  ,
  keep_samples_ours,
  drop = FALSE
]

################### Load GSE271664 ########################################
### 51 MCL, raw HTSeq counts 

getGEOSuppFiles(
  "GSE271664",
  makeDirectory = TRUE,
  baseDir = ".",
  filter_regex = "HTSeq_counts"
)

counts_271664 <- read.csv(
  "GSE271664/GSE271664_HTSeq_counts.csv.gz",
  check.names = FALSE
)

counts_271664_mat <- counts_271664 %>%
  column_to_rownames("Gene") %>%
  as.matrix()

storage.mode(counts_271664_mat) <- "integer"

################### Load GSE305144 ########################################
### 6 primary MCL, TPM 

getGEOSuppFiles(
  "GSE305144",
  makeDirectory = TRUE,
  baseDir = ".",
  filter_regex = "primaryMCL_TPM"
)

tpm_305144 <- read.table(
  "GSE305144/GSE305144_primaryMCL_TPM.txt.gz",
  header = TRUE,
  sep = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

tpm_305144_mat <- tpm_305144 %>%
  filter(
    !is.na(Gene_Name),
    Gene_Name != ""
  ) %>%
  dplyr::select(
    Gene_Name,
    `25021R-01-01`,
    `25021R-01-02`,
    `25021R-01-03`,
    `25021R-01-04`,
    `25021R-01-05`,
    `25021R-01-06`
  ) %>%
  group_by(Gene_Name) %>%
  summarise(
    across(
      everything(),
      ~ median(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  column_to_rownames("Gene_Name") %>%
  as.matrix()

storage.mode(tpm_305144_mat) <- "numeric"
dim(tpm_305144_mat)

################### Load GSE141335 ########################################
### Primary MCL, TPM

getGEOSuppFiles(
  "GSE141335",
  makeDirectory = TRUE,
  baseDir = ".",
  filter_regex = "PT_RNAseq"
)

tpm_141335_part1 <- read_excel(
  "GSE141335/GSE141335_PT_RNAseq_part1.xlsx",
  sheet = 1
)

tpm_141335_part2 <- read_excel(
  "GSE141335/GSE141335_PT_RNAseq_part2.xlsx",
  sheet = 1
)

colnames(tpm_141335_part1)[1] <- "GeneName"
colnames(tpm_141335_part2)[1] <- "GeneName"

# remove BRD4-treated samples
samples_part1_untreated <- setdiff(
  colnames(tpm_141335_part1),
  "GeneName"
)

samples_part1_untreated <- samples_part1_untreated[
  !grepl(
    "BRD4|NVP2|DMSO",
    samples_part1_untreated,
    ignore.case = TRUE
  )
]

# remove DMSO and drug-treated samples
samples_part2_untreated <- setdiff(
  colnames(tpm_141335_part2),
  "GeneName"
)

samples_part2_untreated <- samples_part2_untreated[
  !grepl(
    "BRD4|NVP2|DMSO",
    samples_part2_untreated,
    ignore.case = TRUE
  )
]

part1_baseline <- tpm_141335_part1 %>%
  dplyr::select(
    GeneName,
    all_of(samples_part1_untreated)
  )

part2_baseline <- tpm_141335_part2 %>%
  dplyr::select(
    GeneName,
    all_of(samples_part2_untreated)
  )

# Gene rows occur in different orders in the two workbooks,
# therefore merge explicitly by GeneName

tpm_141335 <- full_join(
  part1_baseline,
  part2_baseline,
  by = "GeneName"
)

### Convert to gene x sample matrix 
tpm_141335_mat <- tpm_141335 %>%
  filter(
    !is.na(GeneName),
    GeneName != ""
  ) %>%
  group_by(GeneName) %>%
  summarise(
    across(
      everything(),
      ~ median(as.numeric(.x), na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  column_to_rownames("GeneName") %>%
  as.matrix()

storage.mode(tpm_141335_mat) <- "numeric"

################### Create common gene universe #############################
common_genes <- Reduce(
  intersect,
  list(
    rownames(counts_ours),
    rownames(counts_271664_mat),
    rownames(tpm_305144_mat),
    rownames(tpm_141335_mat)
  )
)

### Restrict all datasets
counts_ours_common <- counts_ours[
  common_genes,
  ,
  drop = FALSE
]

counts_271664_common <- counts_271664_mat[
  common_genes,
  ,
  drop = FALSE
]

tpm_305144_common <- tpm_305144_mat[
  common_genes,
  ,
  drop = FALSE
]

tpm_141335_common <- tpm_141335_mat[
  common_genes,
  ,
  drop = FALSE
]

################### Expression filtering #############################
# Our RNA-seq:
# >=10 counts in at least 3 samples

keep_ours <- rowSums(
  counts_ours_common >= 10
) >= 3

# GSE271664:
# >=10 counts in at least 3 samples

keep_271664 <- rowSums(
  counts_271664_common >= 10
) >= 3

# GSE305144:
# TPM >=1 in at least 2/6 MCL samples

keep_305144 <- rowSums(
  tpm_305144_common >= 1
) >= 2

# GSE141335:
# TPM >=1 in at least 3 MCL samples

keep_141335 <- rowSums(
  tpm_141335_common >= 1
) >= 3

################### Final common ranking universe ###################
keep_final <- (
  keep_ours &
    keep_271664 &
    keep_305144 &
    keep_141335
)

################### Final matrices ###################################
counts_ours_final <- counts_ours_common[
  keep_final,
  ,
  drop = FALSE
]

counts_271664_final <- counts_271664_common[
  keep_final,
  ,
  drop = FALSE
]

tpm_305144_final <- tpm_305144_common[
  keep_final,
  ,
  drop = FALSE
]

tpm_141335_final <- tpm_141335_common[
  keep_final,
  ,
  drop = FALSE
]

################### Normalize count datasets #########################
### Our cohort 
dds_ours <- DESeqDataSetFromMatrix(
  countData = counts_ours_final,
  colData = data.frame(
    row.names = colnames(counts_ours_final)
  ),
  design = ~ 1
)

dds_ours <- estimateSizeFactors(
  dds_ours
)

vsd_ours <- vst(
  dds_ours,
  blind = TRUE
)

expr_ours <- assay(
  vsd_ours
)

### GSE271664 
dds_271664 <- DESeqDataSetFromMatrix(
  countData = counts_271664_final,
  colData = data.frame(
    row.names = colnames(counts_271664_final)
  ),
  design = ~ 1
)

dds_271664 <- estimateSizeFactors(
  dds_271664
)

vsd_271664 <- vst(
  dds_271664,
  blind = TRUE
)

expr_271664 <- assay(
  vsd_271664
)

################### Prepare TPM datasets #############################
# log2 transformation - does not change ranks, but produces conventional expression values

expr_305144 <- log2(
  tpm_305144_final + 1
)

expr_141335 <- log2(
  tpm_141335_final + 1
)

################### Common chr19 signature ###########################
mcl_signature_common <- intersect(
  mcl_signature,
  rownames(expr_ours)
)

################### Calculate singscores in MCL samples ########################
signature_set <- GeneSet(
  mcl_signature_common,
  geneIdType = SymbolIdentifier(),
  setName = "Chr19_MCL_signature"
)

rank_ours <- rankGenes(
  expr_ours
)

rank_271664 <- rankGenes(
  expr_271664
)

rank_305144 <- rankGenes(
  expr_305144
)

rank_141335 <- rankGenes(
  expr_141335
)

scores_ours <- simpleScore(
  rank_ours,
  upSet = signature_set,
  centerScore = TRUE
)

scores_271664 <- simpleScore(
  rank_271664,
  upSet = signature_set,
  centerScore = TRUE
)

scores_305144 <- simpleScore(
  rank_305144,
  upSet = signature_set,
  centerScore = TRUE
)

scores_141335 <- simpleScore(
  rank_141335,
  upSet = signature_set,
  centerScore = TRUE
)

score_ours <- scores_ours %>%
  rownames_to_column("sample") %>%
  transmute(
    sample,
    score = TotalScore,
    group = case_when(
      grepl(
        "^naive_blood_",
        sample
      ) ~ "Naive B cells",
      
      grepl(
        "^MCL_",
        sample
      ) ~ "Our MCL",
      
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(group)
  )

score_271664 <- scores_271664 %>%
  rownames_to_column("sample") %>%
  transmute(
    sample,
    score = TotalScore,
    group = "GSE271664 MCL"
  )

score_305144 <- scores_305144 %>%
  rownames_to_column("sample") %>%
  transmute(
    sample,
    score = TotalScore,
    group = "GSE305144 MCL"
  )

score_141335 <- scores_141335 %>%
  rownames_to_column("sample") %>%
  transmute(
    sample,
    score = TotalScore,
    group = "GSE141335 MCL"
  )

score_combined <- bind_rows(
  score_ours,
  score_271664,
  score_305144,
  score_141335
)

score_combined$group <- factor(
  score_combined$group,
  levels = c(
    "Naive B cells",
    "Our MCL",
    "GSE271664 MCL",
    "GSE305144 MCL",
    "GSE141335 MCL"
  )
)

################### PLOT #############################################

comparisons_all <- list(
  c("Naive B cells", "Our MCL"),
  c("Naive B cells", "GSE271664 MCL"),
  c("Naive B cells", "GSE305144 MCL"),
  c("Naive B cells", "GSE141335 MCL")
)

cols_combined <- c(
  "Naive B cells" = "#8ecae6",
  "Our MCL"       = "#CE4441",
  "GSE271664 MCL" = "#CE4441",
  "GSE305144 MCL" = "#CE4441",
  "GSE141335 MCL" = "#CE4441"
)

p_combined <- ggplot(
  score_combined,
  aes(
    x = group,
    y = score,
    fill = group
  )
) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    size = 2,
    alpha = 0.7
  ) +
  scale_fill_manual(
    values = cols_combined
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  labs(
    x = NULL,
    y = "Chr19 signature singscore"
  ) +
  stat_compare_means(
    comparisons = comparisons_all,
    method = "wilcox.test",
    label = "p.format"
  )

p_combined

ggsave(
  filename = paste0(
    output_dir,
    "chr19_signature_singscore_validation.pdf"
  ),
  plot = p_combined,
  width = 10,
  height = 14,
  units = "cm"
)


################### GSE132929: CROSS-LYMPHOMA ANALYSIS ###############
gse_132929 <- getGEO(
  "GSE132929",
  GSEMatrix = TRUE,
  getGPL = FALSE
)

eset_132929 <- gse_132929[[1]]

expr_132929 <- exprs(eset_132929)
pheno_132929 <- pData(eset_132929)

### Prepare diagnosis metadata 
meta_132929 <- pheno_132929 %>%
  rownames_to_column("sample") %>%
  transmute(
    sample,
    diagnosis = sub(
      "^diagnosis: ",
      "",
      characteristics_ch1
    )
  )

### Map probes to gene symbols
probe_annot_132929 <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = rownames(expr_132929),
  columns = "SYMBOL",
  keytype = "PROBEID"
) %>%
  filter(
    !is.na(SYMBOL),
    SYMBOL != ""
  ) %>%
  distinct(
    PROBEID,
    SYMBOL
  )

### Collapse probes to gene level 
expr_132929_gene_mat <- as.data.frame(expr_132929) %>%
  rownames_to_column("PROBEID") %>%
  inner_join(
    probe_annot_132929,
    by = "PROBEID"
  ) %>%
  dplyr::select(
    SYMBOL,
    all_of(colnames(expr_132929))
  ) %>%
  group_by(SYMBOL) %>%
  summarise(
    across(
      everything(),
      ~ median(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  column_to_rownames("SYMBOL") %>%
  as.matrix()

storage.mode(expr_132929_gene_mat) <- "numeric"

### Signature coverage 

signature_132929 <- intersect(
  mcl_signature,
  rownames(expr_132929_gene_mat)
)

### Calculate singscore 

signature_set_132929 <- GeneSet(
  signature_132929,
  geneIdType = SymbolIdentifier(),
  setName = "Chr19_MCL_signature"
)

rank_132929 <- rankGenes(
  expr_132929_gene_mat
)

scores_132929 <- simpleScore(
  rank_132929,
  upSet = signature_set_132929,
  centerScore = TRUE
)

score_132929 <- scores_132929 %>%
  rownames_to_column("sample") %>%
  transmute(
    sample,
    score = TotalScore
  ) %>%
  left_join(
    meta_132929,
    by = "sample"
  )

# Exclude Double Hit Lymphoma because n = 1
score_132929_plot <- score_132929 %>%
  filter(
    diagnosis != "Double Hit Lymphoma"
  )

################### Plot #############################################
score_132929_plot <- score_132929_plot %>%
  mutate(
    diagnosis_plot = case_when(
      diagnosis == "Burkitt's Lymphoma" ~ "Burkitt",
      
      diagnosis == "Diffuse Large B-cell Lymphoma" ~ "DLBCL",
      
      diagnosis == "Follicular Lymphoma" ~ "FL",
      
      diagnosis == "Marginal Zone Lymphoma" ~ "MZL",
      
      diagnosis ==
        "High-grade B-cell Lymphoma Not Otherwise Specified" ~
        "HGBCL-NOS",
      
      diagnosis == "Mantle Cell Lymphoma" ~ "MCL",
      
      TRUE ~ diagnosis
    )
  )


################### Plot #############################################
score_132929_plot$diagnosis_plot <- factor(
  score_132929_plot$diagnosis_plot,
  levels = c(
    "MCL",
    "FL",
    "MZL",
    "HGBCL-NOS",
    "DLBCL",
    "Burkitt"
  )
)

cols_132929 <- c(
  "MCL"       = "#CE4441",
  "FL"        = "#BD6EA1",
  "MZL"       = "#DE813F",
  "HGBCL-NOS" = "#294F62",
  "DLBCL"     = "#8A9B73",
  "Burkitt"   = "#EFBC5C"
)

comparisons_132929 <- list(
  c("MCL", "FL"),
  c("MCL", "MZL"),
  c("MCL", "HGBCL-NOS"),
  c("MCL", "DLBCL"),
  c("MCL", "Burkitt")
)


p_132929 <- ggplot(
  score_132929_plot,
  aes(
    x = diagnosis_plot,
    y = score,
    fill = diagnosis_plot
  )
) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    size = 1.5,
    alpha = 0.6
  ) +
  scale_fill_manual(
    values = cols_132929
  ) +
  stat_compare_means(
    comparisons = comparisons_132929,
    method = "wilcox.test",
    label = "p.format"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  labs(
    x = NULL,
    y = "Chr19 MCL signature score"
  )

p_132929

ggsave(
  filename = paste0(
    output_dir,
    "GSE132929_chr19_signature_cross_lymphoma.pdf"
  ),
  plot = p_132929,
  width = 10,
  height = 12,
  units = "cm"
)


