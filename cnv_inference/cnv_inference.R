##################################
# Anna Schwager
# CNRS UMR9018, Institut Gustave Roussy 
# 2026 
##################################

library(CaSpER)
library(dplyr)
library(tibble)
library(data.table)
library(stringr)
library(ggplot2)
library(patchwork)
library(ggpubr)


################### Loading the inputs ##############################
setwd("/Users/annaschwager/Documents/projects/MCL/revision_analysis/rnaseq/input")
output_dir = "/Users/annaschwager/Documents/projects/MCL/revision_analysis/rnaseq/output/"
counts = read.csv2("counts.csv", header = T)

################### Prepare inputs ##############################
counts_casper <- counts %>%
  dplyr::select(
    X,
    starts_with("GRANTA_"),
    starts_with("MCL_"),
    starts_with("naive_blood_")
  ) %>%
  tibble::column_to_rownames("X")

counts_casper$GRANTA_merged <-
  counts_casper$GRANTA_1 +
  counts_casper$GRANTA_2 +
  counts_casper$GRANTA_3

counts_casper <- as.matrix(counts_casper)
storage.mode(counts_casper) <- "numeric"

control_samples <- grep(
  "^naive_blood_",
  colnames(counts_casper),
  value = TRUE
)

tumor_samples <- c(
  grep("^MCL_", colnames(counts_casper), value = TRUE),
  grep("^GRANTA_", colnames(counts_casper), value = TRUE)
)

loh <- readBAFExtractOutput(
  path = "baf_output",
  sequencing.type = "bulk",
  suffix = "snp"
)

names(loh) <- sub("\\.snp$", "", names(loh))

loh.name.mapping <- data.frame(
  loh.name = names(loh),
  sample.name = names(loh),
  stringsAsFactors = FALSE
)

################### Prepare annotation ##############################
data("hg38_cytoband", package = "CaSpER")

annotation <- generateAnnotation(
  id_type = "hgnc_symbol",
  genes = rownames(counts_casper),
  ishg19 = FALSE,
  centromere = centromere_hg38
)

#saveRDS(annotation,file = "annotation.rds")

annotation <- readRDS("annotation.rds")
common_genes <- intersect(annotation$Gene, rownames(counts_casper))

counts_casper <- counts_casper[
  common_genes,
  ,
  drop = FALSE
]

annotation <- annotation[
  match(common_genes, annotation$Gene),
  ,
  drop = FALSE
]

counts_casper <- counts_casper[
  annotation$Gene,
  ,
  drop = FALSE
]

sample_order <- colnames(counts_casper)
loh <- loh[sample_order]

################### Create and run Casper ##############################
hist(log2(rowMeans(counts_casper) + 1),
     breaks = 100,
     main = "Mean expression",
     xlab = "log2(mean counts + 1)")
#expression cutoff should be 0.1-1

casper_object <- CreateCasperObject(
  raw.data = counts_casper,
  annotation = annotation,
  control.sample.ids = control_samples,
  cytoband = cytoband_hg38,
  loh.name.mapping = loh.name.mapping,
  cnv.scale = 3, 
  loh.scale = 3,
  method = "iterative",
  loh = loh,
  project = "MCL_CaSpER",
  matrix.type = "raw",
  sequencing.type = "bulk",
  expr.cutoff = 0.2,
  log.transformed = FALSE,
  genomeVersion = "hg38",
  filter = "median"
)

casper_object <- runCaSpER(casper_object)

saveRDS(
  casper_object,
  file = "MCL_CaSpER_object.rds"
)

################### Plotting ##############################
casper_object <- readRDS("MCL_CaSpER_object.rds")

cna_calls <- extractLargeScaleEvents(casper_object, thr = 0.85)

sample_order <- c(
  grep("^naive_blood_", rownames(cna_calls), value = TRUE),
  grep("^MCL_our_",     rownames(cna_calls), value = TRUE),
  grep("^MCL_DAG_",     rownames(cna_calls), value = TRUE),
  "GRANTA_merged"
)

pdf(paste(output_dir,"CaSpER_arm_level_CNA_heatmap.pdf"), width = 11, height = 3.7)

pheatmap::pheatmap(
  cna_calls[sample_order, , drop = FALSE],
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = c("#3B6FB6", "white", "#C84545"),
  breaks = c(-1.5, -0.5, 0.5, 1.5),
  border_color = "black"
)

dev.off()

pdf(paste(output_dir,"CaSpER_arm_level_CNA_heatmap_no_zeros.pdf"), width = 7, height = 3.7)

plot_mat <- cna_calls[sample_order, , drop = FALSE]
plot_mat <- plot_mat[, colSums(plot_mat != 0) > 0, drop = FALSE]

pheatmap::pheatmap(
  plot_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = c("#3B6FB6", "white", "#C84545"),
  breaks = c(-1.5, -0.5, 0.5, 1.5),
  border_color = "black"
)

dev.off()


############# RNA-derived BAF analysis across chromosome 19p  #################
# Outputs:
#   1. Rolling median |BAF - 0.5| across chromosome 19p
#   2. Sample-level median |BAF - 0.5| around the probe,
#      with pairwise Wilcoxon P-values

### Load files 
baf_dir <- "baf_output"

snp_files <- list.files(
  path = baf_dir,
  pattern = "\\.snp$",
  full.names = TRUE
)

read_baf_file <- function(file) {
  
  x <- data.table::fread(
    file,
    header = FALSE,
    col.names = c(
      "chr",
      "position",
      "ref_allele",
      "alt_allele",
      "alt_count",
      "coverage"
    )
  )
  
  x %>%
    mutate(
      sample = sub(
        "\\.snp$",
        "",
        basename(file)
      ),
      
      chr = sub(
        "^chr",
        "",
        as.character(chr)
      ),
      
      position = as.numeric(position),
      alt_count = as.numeric(alt_count),
      coverage = as.numeric(coverage),
      
      ref_count = coverage - alt_count,
      BAF = alt_count / coverage
    ) %>%
    filter(
      coverage > 0,
      alt_count >= 0,
      ref_count >= 0,
      is.finite(BAF),
      between(BAF, 0, 1)
    )
}

baf <- bind_rows(
  lapply(
    snp_files,
    read_baf_file
  )
)

print(summary(baf$BAF))

baf <- baf %>%
  mutate(
    display_group = case_when(
      str_detect(sample, "^naive_blood_") ~ "Naive B cells",
      
      str_detect(sample, "^MCL_DAG_") |
        str_detect(sample, "^MCL_our_") ~ "Primary MCL",
      
      str_detect(sample, "^GRANTA_") ~ "GRANTA",
      
      TRUE ~ NA_character_
    ),
    
    display_group = factor(
      display_group,
      levels = c(
        "Naive B cells",
        "Primary MCL",
        "GRANTA"
      )
    )
  ) %>%
  filter(
    !is.na(display_group),
    
    # Exclude the merged GRANTA dataset because it duplicates
    # the individual GRANTA libraries
    sample != "GRANTA_merged"
  )

sample_annotation <- baf %>%
  distinct(
    display_group,
    sample
  ) %>%
  arrange(
    display_group,
    sample
  )

print(sample_annotation)


### Set parameters
# Probe coordinates, hg38
probe_chr   <- "19"
probe_start <- 478637
probe_end   <- 702132

# Chromosome 19p interval, hg38
region_start <- 1
region_end   <- 24400000

# Downsample all retained loci to the same read depth
target_depth <- 30

# Rolling-profile parameters
window_n <- 301
step_n   <- 25

grid_step <- 200000
max_interpolation_gap <- 1000000

# Sample-level summary interval: probe ±1 Mb
summary_flank <- 1000000

summary_start <- max(
  region_start,
  probe_start - summary_flank
)

summary_end <- min(
  region_end,
  probe_end + summary_flank
)

downsampling_seed <- 123

profile_colours <- c(
  "Naive B cells" = "#56B9D5",
  "Primary MCL"   = "#E75B55",
  "GRANTA"        = "#665099"
)

output_pdf <- "chr19p_BAF_profile_and_probe_summary.pdf"

### Retain chromosome 19p loci with depth >= target depth
baf_19p <- baf %>%
  filter(
    chr == probe_chr,
    position >= region_start,
    position <= region_end,
    coverage >= target_depth
  ) %>%
  arrange(
    sample,
    position
  )

### QC
sample_qc <- baf_19p %>%
  group_by(
    display_group,
    sample
  ) %>%
  summarise(
    n_snps = n(),
    median_depth = median(
      coverage,
      na.rm = TRUE
    ),
    first_position = min(
      position,
      na.rm = TRUE
    ),
    last_position = max(
      position,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

print(sample_qc)


probe_qc <- baf_19p %>%
  filter(
    position >= probe_start,
    position <= probe_end
  ) %>%
  count(
    display_group,
    sample,
    name = "n_probe_snps"
  )

probe_qc <- sample_annotation %>%
  left_join(
    probe_qc,
    by = c(
      "display_group",
      "sample"
    )
  ) %>%
  mutate(
    n_probe_snps = coalesce(
      n_probe_snps,
      0L
    )
  )

print(probe_qc)


### Equalize depth (downsample every SNP to target_depth)
set.seed(downsampling_seed)

alt_count_downsampled <- mapply(
  FUN = function(alt, ref) {
    
    stats::rhyper(
      nn = 1,
      m = alt,
      n = ref,
      k = target_depth
    )
  },
  
  alt = baf_19p$alt_count,
  ref = baf_19p$ref_count
)

baf_19p_ds <- baf_19p %>%
  mutate(
    alt_count_ds = as.numeric(
      alt_count_downsampled
    ),
    
    ref_count_ds =
      target_depth - alt_count_ds,
    
    BAF_ds =
      alt_count_ds / target_depth,
    
    baf_deviation =
      abs(BAF_ds - 0.5)
  )

### Edge-adjusted rolling-median function
calculate_rolling_baf <- function(
    df,
    window_n = 301,
    step_n = 25
) {
  
  df <- df %>%
    arrange(position)
  
  n_total <- nrow(df)
  
  if (n_total < window_n) {
    return(tibble())
  }
  
  half_window <- floor(
    window_n / 2
  )
  
  focal_indices <- sort(
    unique(
      c(
        1L,
        
        seq(
          from = 1L,
          to = n_total,
          by = step_n
        ),
        
        n_total
      )
    )
  )
  
  bind_rows(
    lapply(
      focal_indices,
      function(i) {
        
        if (i <= half_window) {
          
          # Forward window at the telomeric edge
          first_index <- 1L
          last_index  <- window_n
          
        } else if (
          i > n_total - half_window
        ) {
          
          # Backward window near the centromeric edge
          first_index <-
            n_total - window_n + 1L
          
          last_index <- n_total
          
        } else {
          
          # Centred window elsewhere
          first_index <-
            i - half_window
          
          last_index <-
            i + half_window
        }
        
        idx <- first_index:last_index
        
        tibble(
          # Plot the estimate at the focal SNP position
          position = df$position[i],
          
          median_baf_deviation = median(
            df$baf_deviation[idx],
            na.rm = TRUE
          ),
          
          window_start = min(
            df$position[idx],
            na.rm = TRUE
          ),
          
          window_end = max(
            df$position[idx],
            na.rm = TRUE
          ),
          
          n_snps = length(idx)
        )
      }
    )
  )
}

rolling_baf <- baf_19p_ds %>%
  group_by(
    sample,
    display_group
  ) %>%
  group_modify(
    ~ calculate_rolling_baf(
      df = .x,
      window_n = window_n,
      step_n = step_n
    )
  ) %>%
  ungroup()


### Add edge anchors so profiles span the displayed 19p arm
rolling_baf_extended <- rolling_baf %>%
  group_by(
    sample,
    display_group
  ) %>%
  group_modify(
    ~ bind_rows(
      tibble(
        position = region_start,
        
        median_baf_deviation =
          first(.x$median_baf_deviation),
        
        window_start =
          first(.x$window_start),
        
        window_end =
          first(.x$window_end),
        
        n_snps =
          first(.x$n_snps)
      ),
      
      .x,
      
      tibble(
        position = region_end,
        
        median_baf_deviation =
          last(.x$median_baf_deviation),
        
        window_start =
          last(.x$window_start),
        
        window_end =
          last(.x$window_end),
        
        n_snps =
          last(.x$n_snps)
      )
    )
  ) %>%
  ungroup() %>%
  arrange(
    sample,
    position
  )


### Interpolate sample profiles onto a common genomic grid
profile_grid <- seq(
  from = region_start,
  to = region_end,
  by = grid_step
)

profile_grid <- sort(
  unique(
    c(
      profile_grid,
      probe_start,
      probe_end,
      region_end
    )
  )
)

interpolate_profile <- function(
    df,
    grid,
    max_gap
) {
  
  df <- df %>%
    arrange(position) %>%
    distinct(
      position,
      .keep_all = TRUE
    )
  
  if (nrow(df) < 2) {
    return(tibble())
  }
  
  interpolation <- stats::approx(
    x = df$position,
    y = df$median_baf_deviation,
    xout = grid,
    rule = 1,
    ties = "ordered"
  )
  
  nearest_distance <- vapply(
    grid,
    function(grid_position) {
      
      min(
        abs(
          df$position - grid_position
        )
      )
    },
    numeric(1)
  )
  
  profile_values <- interpolation$y
  
  # Prevent interpolation through long regions lacking estimates
  profile_values[
    nearest_distance > max_gap
  ] <- NA_real_
  
  tibble(
    position = interpolation$x,
    median_baf_deviation =
      profile_values
  )
}

rolling_baf_grid <- rolling_baf_extended %>%
  group_by(
    sample,
    display_group
  ) %>%
  group_modify(
    ~ interpolate_profile(
      df = .x,
      grid = profile_grid,
      max_gap = max_interpolation_gap
    )
  ) %>%
  ungroup()

### Calculate cohort-median rolling profiles
cohort_baf_profiles <- rolling_baf_grid %>%
  group_by(
    display_group,
    position
  ) %>%
  summarise(
    n_samples = sum(
      !is.na(median_baf_deviation)
    ),
    
    cohort_median = if (
      n_samples > 0
    ) {
      median(
        median_baf_deviation,
        na.rm = TRUE
      )
    } else {
      NA_real_
    },
    
    .groups = "drop"
  )

### Calculate one probe-region value per sample

probe_sample_summary <- baf_19p_ds %>%
  filter(
    position >= summary_start,
    position <= summary_end
  ) %>%
  group_by(
    display_group,
    sample
  ) %>%
  summarise(
    n_snps = n(),
    
    median_baf_deviation = median(
      baf_deviation,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    display_group = factor(
      display_group,
      levels = c(
        "Naive B cells",
        "Primary MCL",
        "GRANTA"
      )
    )
  )

print(
  probe_sample_summary %>%
    arrange(
      display_group,
      sample
    )
)

wilcox_comparisons <- list(
  c(
    "Naive B cells",
    "Primary MCL"
  ),
  c(
    "Naive B cells",
    "GRANTA"
  ),
  c(
    "Primary MCL",
    "GRANTA"
  )
)

### Final plots
p_baf_19p <- ggplot() +
  
  annotate(
    geom = "rect",
    xmin = probe_start / 1e6,
    xmax = probe_end / 1e6,
    ymin = -Inf,
    ymax = Inf,
    fill = "grey60",
    alpha = 0.25
  ) +
  
  # Individual samples
  geom_line(
    data = rolling_baf_grid,
    mapping = aes(
      x = position / 1e6,
      y = median_baf_deviation,
      group = sample,
      colour = display_group
    ),
    linewidth = 0.30,
    alpha = 0.10,
    na.rm = TRUE
  ) +
  
  # Cohort medians
  geom_line(
    data = cohort_baf_profiles,
    mapping = aes(
      x = position / 1e6,
      y = cohort_median,
      colour = display_group
    ),
    linewidth = 1.20,
    na.rm = TRUE
  ) +
  
  scale_colour_manual(
    values = profile_colours,
    name = NULL
  ) +
  
  scale_x_continuous(
    limits = c(
      region_start,
      region_end
    ) / 1e6,
    
    breaks = seq(
      0,
      24,
      by = 2
    ),
    
    expand = expansion(
      mult = c(0, 0.005)
    )
  ) +
  
  scale_y_continuous(
    limits = c(0, 0.5),
    
    breaks = seq(
      0,
      0.5,
      by = 0.1
    ),
    
    expand = expansion(
      mult = c(0, 0.02)
    )
  ) +
  
  theme_classic(
    base_size = 11
  ) +
  
  theme(
    legend.position = c(
      0.81,
      0.84
    ),
    
    legend.background =
      element_blank(),
    
    plot.title = element_text(
      size = 13
    ),
    
    plot.subtitle = element_text(
      size = 9.5
    )
  ) +
  
  labs(
    title =
      "RNA-derived allelic imbalance across chromosome 19p",
    
    subtitle = paste0(
      "Rolling median over ",
      window_n,
      " consecutive SNPs after downsampling each locus to "
    ),
    
    x = "Chromosome 19 position (Mb)",
    
    y = expression(
      "Median " * "|" * BAF - 0.5 * "|"
    )
  )


wilcox_label_positions <- c(
  0.38,
  0.43,
  0.48
)

p_probe_summary <- ggplot(
  probe_sample_summary,
  aes(
    x = display_group,
    y = median_baf_deviation,
    colour = display_group
  )
) +
  
  geom_point(
    position = position_jitter(
      width = 0.08,
      height = 0,
      seed = 123
    ),
    size = 2.7
  ) +
  
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.50,
    linewidth = 0.75
  ) +
  
  ggpubr::stat_compare_means(
    comparisons = wilcox_comparisons,
    method = "wilcox.test",
    method.args = list(
      exact = FALSE
    ),
    label = "p.format",
    label.y = wilcox_label_positions,
    tip.length = 0.01,
    size = 3.4
  ) +
  
  scale_colour_manual(
    values = profile_colours,
    guide = "none"
  ) +
  
  scale_y_continuous(
    limits = c(0, 0.5),
    
    breaks = seq(
      0,
      0.5,
      by = 0.1
    ),
    
    expand = expansion(
      mult = c(0, 0.02)
    )
  ) +
  
  theme_classic(
    base_size = 11
  ) +
  
  theme(
    axis.title.x =
      element_blank(),
    
    plot.title = element_text(
      size = 13
    ),
    
    plot.subtitle = element_text(
      size = 9.5
    )
  ) +
  
  labs(
    title =
      "Allelic imbalance around the chromosome 19 probe",
    
    subtitle = paste0(
      "Interval: ",
      round(summary_start / 1e6, 2),
      "–",
      round(summary_end / 1e6, 2),
      " Mb"
    ),
    
    y = expression(
      "Median " * "|" * BAF - 0.5 * "|"
    )
  )

final_baf_plot <- (
  p_baf_19p | p_probe_summary
) +
  patchwork::plot_layout(
    widths = c(5, 1)
  )

print(final_baf_plot)

ggsave(
  filename = paste(output_dir, output_pdf),
  plot = final_baf_plot,
  width = 28,
  height = 8,
  units = "cm",
  device = cairo_pdf
)


