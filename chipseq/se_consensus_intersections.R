##################################
# Anna Schwager
# CNRS UMR9018, Institut Gustave Roussy 
# 2024 
##################################
library(ggvenn)
library(ChIPpeakAnno)
library(plyr)
library(EnsDb.Hsapiens.v86)
library(biomaRt)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(clusterProfiler)

############### Consensus and condition-specific SEs ##################
### Loading the SEs as GRanges ###
mcl <- lapply(paste("mcl/", list.files("mcl"), sep = ""), toGRanges, format = "BED", header = FALSE)
naive <- lapply(paste("naive/", list.files("naive"), sep = ""), toGRanges, format = "BED", header = FALSE)
germinal <- lapply(paste("germinal/", list.files("germinal"), sep = ""), toGRanges, format = "BED", header = FALSE)

## Finding consensus SEs with at least 80% region overlap within condition
mcl_overlap <- findOverlapsOfPeaks(mcl[[1]], mcl[[2]], mcl[[3]], mcl[[4]], mcl[[5]], minoverlap = 0.8)
mcl_ol_gr <- mcl_overlap$mergedPeaks #1121 peaks
write_csv(as.data.frame(mcl_ol_gr), "results/mcl_consensus_se.csv")

naive_overlap <- findOverlapsOfPeaks(naive[[1]], naive[[2]], naive[[3]], minoverlap = 0.8)
naive_ol_gr <- naive_overlap$mergedPeaks #752 peaks
write_csv(as.data.frame(naive_ol_gr), "results/naive_consensus_se.csv")

germinal_overlap <- findOverlapsOfPeaks(germinal[[1]], germinal[[2]], minoverlap = 0.8)
germinal_ol_gr <- germinal_overlap$mergedPeaks # 520 peaks
write_csv(as.data.frame(germinal_ol_gr), "results/germinal_consensus_se.csv")

## Median sizes of SEs ##
dfm <- as.data.frame(mcl_ol_gr)
dfn <- as.data.frame(naive_ol_gr)
dfg <- as.data.frame(germinal_ol_gr)
sizes <- data.frame(condition = c(rep("mcl", length(dfm$width)),
                                  rep("naive", length(dfn$width)),
                                  rep("germinal", length(dfg$width))),
                                  size = c(dfm$width, dfn$width, dfg$width))
stats <- ddply(sizes, "condition", summarise, size.median=median(size))
stats

#condition size.median
#1  germinal     26195.5
#2       mcl     40871.0
#3     naive     17727.0

median_sizes <- ggplot(sizes, aes(x=size, fill = condition)) +
                geom_density(alpha=.2) +
                geom_vline(data=stats, aes(xintercept=size.median,  colour=condition),
                linetype="dashed") + 
                theme_bw() +
                labs(title="Sizes of SEs + median") +
                xlim(0, 200000)

## Intersecting the consensus SEs between conditions
overlapp_all <- findOverlapsOfPeaks(mcl_ol_gr, naive_ol_gr, germinal_ol_gr, minoverlap = 0.8)
venn <- makeVennDiagram(overlapp_all)

unique_mcl <- overlapp_all$peaklist$mcl_ol_gr
common_for_all <- overlapp_all$peaklist$`mcl_ol_gr///naive_ol_gr///germinal_ol_gr`
common_mcl_naive <- overlapp_all$peaklist$`mcl_ol_gr///naive_ol_gr`
common_mcl_germinal <- overlapp_all$peaklist$`mcl_ol_gr///germinal_ol_gr`
write_csv(as.data.frame(common_for_all), "results/common_for_all.csv")
write_csv(as.data.frame(common_mcl_naive), "results/common_mcl_naive.csv")
write_csv(as.data.frame(common_mcl_germinal), "results/common_mcl_germinal.csv")
write_csv(as.data.frame(unique_mcl), "results/unique_mcl.csv")

export.bed(unique_mcl, "results/unique_se_mcl.bed")
export.bed(common_for_all, "results/common_for_all.bed")
export.bed(common_mcl_naive, "results/common_mcl_naive.bed")
export.bed(common_mcl_germinal, "results/common_mcl_germinal.bed")

########################## Annotating SEs #####################
### Annotating SEs ###
## preparing annotation files
anno_ensdb <- toGRanges(EnsDb.Hsapiens.v86)
## annotating peaksets of interest 
common_for_all_annot <- annotatePeakInBatch(common_for_all,
                                            output = "nearestLocation",
                                            PeakLocForDistance = "middle",
                                            FeatureLocForDistance = "TSS",
                                            AnnotationData = anno_ensdb)
unique_mcl_annot <- annotatePeakInBatch(unique_mcl,
                                  output = "nearestLocation",
                                  PeakLocForDistance = "middle",
                                  FeatureLocForDistance = "TSS",
                                  AnnotationData = anno_ensdb)
## adding gene symbols 
common_for_all_annot_df <- as.data.frame(common_for_all_annot)
symbols <- addGeneIDs(common_for_all_annot_df$feature, 
                      orgAnn = "org.Hs.eg.db", 
                      feature_id_type = "ensembl_gene_id",
                      IDs2Add = "symbol")
colnames(symbols) <- c("feature", "symbol")
common_for_all_annot_df <- left_join(common_for_all_annot_df, symbols)

unique_mcl_annot_df <- as.data.frame(unique_mcl_annot)
symbols <- addGeneIDs(unique_mcl_annot$feature, 
                      orgAnn = "org.Hs.eg.db", 
                      feature_id_type = "ensembl_gene_id",
                      IDs2Add = "symbol")
colnames(symbols) <- c("feature", "symbol")
unique_mcl_annot_df <- left_join(unique_mcl_annot_df, symbols)

write_csv(common_for_all_annot_df, "results/common_for_all_annot.csv")
write_csv(unique_mcl_annot_df, "results/unique_mcl_annot.csv")


############# Enrichment analysis for SE-associated genes ###############
## GO over-representation analysis ## 
common_for_all_go <- getEnrichedGO(common_for_all_annot, 
                             orgAnn = "org.Hs.eg.db", 
                             feature_id_type = "ensembl_gene_id")

unique_mcl_go <- getEnrichedGO(unique_mcl_annot, 
                                   orgAnn = "org.Hs.eg.db", 
                                   feature_id_type = "ensembl_gene_id")

p_common_for_all_go <- enrichmentPlot(common_for_all_go, n = 10)
p_unique_mcl_go <- enrichmentPlot(unique_mcl_go, n = 10)

## adding gene symbols
unique_mcl_bp <- as.data.frame(unique_mcl_go$bp)
symbols <- addGeneIDs(unique_mcl_bp$EntrezID, 
                      orgAnn = "org.Hs.eg.db", 
                      feature_id_type = "entrez_id",
                      IDs2Add = "symbol")
colnames(symbols) <- c("EntrezID", "symbol")
unique_mcl_bp <- left_join(unique_mcl_bp, symbols)

unique_mcl_mf <- as.data.frame(unique_mcl_go$mf)
symbols <- addGeneIDs(unique_mcl_mf$EntrezID, 
                      orgAnn = "org.Hs.eg.db", 
                      feature_id_type = "entrez_id",
                      IDs2Add = "symbol")
colnames(symbols) <- c("EntrezID", "symbol")
unique_mcl_mf <- left_join(unique_mcl_mf, symbols)

unique_mcl_cc <- as.data.frame(unique_mcl_go$cc)
symbols <- addGeneIDs(unique_mcl_mf$EntrezID, 
                      orgAnn = "org.Hs.eg.db", 
                      feature_id_type = "entrez_id",
                      IDs2Add = "symbol")
colnames(symbols) <- c("EntrezID", "symbol")
unique_mcl_cc <- left_join(unique_mcl_cc, symbols)

write.csv2(unique_mcl_bp, "results/go_bp_unique_se_mcl.csv")
write.csv2(unique_mcl_cc, "results/go_cc_unique_se_mcl.csv")
write.csv2(unique_mcl_mf, "results/go_mf_unique_se_mcl.csv")

common_for_all_bp <- as.data.frame(common_for_all_go$bp)
symbols <- addGeneIDs(common_for_all_bp$EntrezID, 
                      orgAnn = "org.Hs.eg.db", 
                      feature_id_type = "entrez_id",
                      IDs2Add = "symbol")
colnames(symbols) <- c("EntrezID", "symbol")
common_for_all_bp <- left_join(common_for_all_bp, symbols)

common_for_all_cc <- as.data.frame(common_for_all_go$cc)
symbols <- addGeneIDs(common_for_all_cc$EntrezID, 
                      orgAnn = "org.Hs.eg.db", 
                      feature_id_type = "entrez_id",
                      IDs2Add = "symbol")
colnames(symbols) <- c("EntrezID", "symbol")
common_for_all_cc <- left_join(common_for_all_cc, symbols)

common_for_all_mf <- as.data.frame(common_for_all_go$mf)
symbols <- addGeneIDs(common_for_all_mf$EntrezID, 
                      orgAnn = "org.Hs.eg.db", 
                      feature_id_type = "entrez_id",
                      IDs2Add = "symbol")
colnames(symbols) <- c("EntrezID", "symbol")
common_for_all_mf <- left_join(common_for_all_mf, symbols)

write.csv2(common_for_all_bp, "results/go_bp_common_se.csv")
write.csv2(common_for_all_cc, "results/go_cc_common_se.csv")
write.csv2(common_for_all_mf, "results/go_mf_common_se.csv")

## Over-representation analysis for the disease gene network
entrez_unique_mcl <- addGeneIDs(unique_mcl_annot_df$feature, 
                      orgAnn = "org.Hs.eg.db", 
                      feature_id_type = "ensembl_gene_id",
                      IDs2Add = "entrez_id")

dgn <- enrichDGN(entrez_unique_mcl$entrez_id,
                 pvalueCutoff = 0.5,
                 minGSSize = 5,
                 qvalueCutoff = 0.5,
                 readable = TRUE)
dotplot(dgn, showCategory=10) + ggtitle("DisGeNET enrichment, unique SEs MCL")
write.csv2(dgn, "results/disgenet_enrichment_unique_mcl.csv")

## Over-representation analysis for the network of cancer gene
ncg <- enrichNCG(entrez$entrez_id,
                 pvalueCutoff = 0.5,
                 minGSSize = 5,
                 qvalueCutoff = 0.5,
                 readable = TRUE)
dotplot(ncg, showCategory=10) + ggtitle("Network of Cancer Gene enrichment, unique SEs MCL")
write.csv2(ncg, "results/ncg_enrichment_unique_mcl.csv")


## KEGG enrichment analysis
kk <- enrichKEGG(entrez$entrez_id,
                 organism     = 'hsa',
                 pvalueCutoff = 0.5,
                 qvalueCutoff = 0.6,
                 minGSSize = 5)
dotplot(kk, showCategory=10) + ggtitle("KEGG enrichment, unique SEs MCL")
write.csv2(kk, "results/kegg_enrichment_unique_mcl.csv")

############# SEs per chromosome ###############
table(unique_mcl_annot_df$seqnames)

chromsizes <- read.csv2("./chromSizes.csv", header = TRUE)
chromsizes$normsizes <- chromsizes$size/min(chromsizes$size)
chromsizes$normnumber <- chromsizes$gene_number/min(chromsizes$gene_number)

se_per_chr <- as.data.frame(table(unique_mcl_annot_df$seqnames))
colnames(se_per_chr) <- c("chr", "freq")
se_per_chr  <- merge(se_per_chr,chromsizes, by.x = "chr", by.y = "chr")
se_per_chr$normfreq_size <- se_per_chr$freq/se_per_chr$normsizes #norm by size
se_per_chr$normfreq_number <- se_per_chr$freq/se_per_chr$normnumber #norm by gene number
write.csv(se_per_chr, "uniqueMCL_SEs_per_chromosome.csv")

se_per_chr <- se_per_chr[order(se_per_chr$freq,decreasing=TRUE),]
barplot(se_per_chr$freq, names.arg = se_per_chr$chr, ylab = "N MCL unique SEs") 
se_per_chr <- se_per_chr[order(se_per_chr$normfreq_size,decreasing=TRUE),]
barplot(se_per_chr$normfreq_size, names.arg = se_per_chr$chr, ylab = "N MCL unique SEs normalised to chromSize patients")
se_per_chr <- se_per_chr[order(se_per_chr$normfreq_number,decreasing=TRUE),]
barplot(se_per_chr$normfreq_number, names.arg = se_per_chr$chr, ylab = "N MCL unique SEs normalised to gene number patients")

###############################################################################
#> sessionInfo()
#R version 4.3.0 (2023-04-21)
#Platform: aarch64-apple-darwin20 (64-bit)
#Running under: macOS Ventura 13.3

#Matrix products: default
#BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
#LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

#locale:
#  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

#time zone: Europe/Paris
#tzcode source: internal

#attached base packages:
#  [1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#[1] rtracklayer_1.60.1          stringr_1.5.1               DESeq2_1.40.2              
#[4] SummarizedExperiment_1.30.2 MatrixGenerics_1.12.3       matrixStats_1.3.0          
#[7] clusterProfiler_4.8.2       enrichplot_1.20.0           DOSE_3.26.2                
#[10] org.Hs.eg.db_3.17.0         EnsDb.Hsapiens.v86_2.99.0   ensembldb_2.24.1           
#[13] AnnotationFilter_1.24.0     GenomicFeatures_1.52.2      AnnotationDbi_1.62.2       
#[16] Biobase_2.60.0              plyr_1.8.9                  ChIPpeakAnno_3.34.1        
#[19] ggrepel_0.9.5               ggvenn_0.1.10               ggplot2_3.5.1              
#[22] dplyr_1.1.4                 biomaRt_2.56.1              GenomicRanges_1.52.1       
#[25] GenomeInfoDb_1.36.4         IRanges_2.34.1              S4Vectors_0.38.2           
#[28] BiocGenerics_0.46.0          

#loaded via a namespace (and not attached):
#  [1] splines_4.3.0               BiocIO_1.10.0               bitops_1.0-7                ggplotify_0.1.2             filelock_1.0.3              tibble_3.2.1               
#[7] polyclip_1.10-6             graph_1.78.0                XML_3.99-0.16.1             lifecycle_1.0.4             lattice_0.22-6              vroom_1.6.5                
#[13] MASS_7.3-60                 magrittr_2.0.3              rmarkdown_2.27              yaml_2.3.8                  cowplot_1.1.3               DBI_1.2.2                  
#[19] RColorBrewer_1.1-3          abind_1.4-7                 zlibbioc_1.46.0             purrr_1.0.2                 ggraph_2.2.1                RCurl_1.98-1.14            
#[25] yulab.utils_0.1.4           tweenr_2.0.3                rappdirs_0.3.3              GenomeInfoDbData_1.2.10     ggrepel_0.9.5               tidytree_0.4.6             
#[31] codetools_0.2-20            DelayedArray_0.26.7         xml2_1.3.6                  ggforce_0.4.2               tidyselect_1.2.1            futile.logger_1.4.3        
#[37] aplot_0.2.2                 farver_2.1.2                viridis_0.6.5               matrixStats_1.3.0           BiocFileCache_2.8.0         GenomicAlignments_1.36.0   
#[43] jsonlite_1.8.8              multtest_2.56.0             tidygraph_1.3.1             survival_3.6-4              tools_4.3.0                 progress_1.2.3             
#[49] treeio_1.24.3               Rcpp_1.0.12                 glue_1.7.0                  gridExtra_2.3               xfun_0.44                   qvalue_2.32.0              
#[55] MatrixGenerics_1.12.3       withr_3.0.0                 formatR_1.14                fastmap_1.2.0               fansi_1.0.6                 digest_0.6.35              
#[61] R6_2.5.1                    gridGraphics_0.5-1          colorspace_2.1-1            GO.db_3.17.0                RSQLite_2.3.7               utf8_1.2.4                 
#[67] tidyr_1.3.1                 generics_0.1.3              data.table_1.15.4           rtracklayer_1.60.1          prettyunits_1.2.0           graphlayouts_1.1.1         
#[73] InteractionSet_1.28.1       httr_1.4.7                  S4Arrays_1.0.6              scatterpie_0.2.2            regioneR_1.32.0             pkgconfig_2.0.3            
#[79] gtable_0.3.5                blob_1.2.4                  XVector_0.40.0              shadowtext_0.1.3            htmltools_0.5.8.1           fgsea_1.26.0               
#[85] RBGL_1.76.0                 ProtGenerics_1.32.0         scales_1.3.0                png_0.1-8                   ggfun_0.1.4                 knitr_1.45                 
#[91] lambda.r_1.2.4              rstudioapi_0.16.0           tzdb_0.4.0                  reshape2_1.4.4              rjson_0.2.21                nlme_3.1-164               
#[97] curl_5.2.1                  cachem_1.1.0                stringr_1.5.1               parallel_4.3.0              HDO.db_0.99.1               restfulr_0.0.15            
#[103] pillar_1.9.0                vctrs_0.6.5                 dbplyr_2.5.0                evaluate_0.23               VennDiagram_1.7.3           cli_3.6.2                  
#[109] compiler_4.3.0              futile.options_1.0.1        Rsamtools_2.16.0            rlang_1.1.3                 crayon_1.5.2                labeling_0.4.3             
#[115] fs_1.6.4                    stringi_1.8.4               viridisLite_0.4.2           BiocParallel_1.34.2         munsell_0.5.1               Biostrings_2.68.1          
#[121] lazyeval_0.2.2              GOSemSim_2.26.1             Matrix_1.6-5                BSgenome_1.68.0             hms_1.1.3                   patchwork_1.2.0            
#[127] bit64_4.0.5                 KEGGREST_1.40.1             SummarizedExperiment_1.30.2 igraph_2.0.3                memoise_2.0.1               ggtree_3.8.2               
#[133] fastmatch_1.1-4             bit_4.0.5                   downloader_0.4              gson_0.1.0                  ape_5.8                    






