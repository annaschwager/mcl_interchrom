##################################
# Anna Schwager
# CNRS UMR9018, Institut Gustave Roussy 
# 2024 
##################################
library(DiffBind)
library(devtools)
library(patchwork)
library(karyoploteR)
library(EnsDb.Hsapiens.v86)
library(csaw)
library(profileplyr)
library(nucleR)
library(Cairo)
library(rtracklayer)
library(ggvenn)
library(VennDiagram)
library(ChIPseeker)
library(clusterProfiler)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

################### Differential binding analysis ##############################
samples <- read.csv2('meta.csv') #loading the samplesheet
dbObj <- dba(sampleSheet = samples)
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps = TRUE)

dba.plotPCA(dbObj, attributes = DBA_CONDITION, label = DBA_ID) #PCA plot using affinity (read count) data
plot(dbObj) #correlation plot using affinity (read count) data
dbObj_normalized <- dba.normalize(dbObj, background=TRUE) 
#normalisation using the default library-size based method (lib) + background normalisation is recommended for ATAC

#establishing contrast
dbObj_cont <- dba.contrast(dbObj_normalized, design="~Condition", minMembers = 2)
#filtering for blacklist
dbObj_cont <- dba.blacklist(dbObj_cont, greylist=FALSE)
#Removed 157 (of 59569) consensus peaks.

#finding differentially bound peaks
dbObj_cont <- dba.analyze(dbObj_cont, bGreylist=FALSE, design="~Condition")
dba.show(dbObj_cont, bContrasts = T) #show summary

#saving results
dba.report(dbObj_cont,  contrast=1)
write.csv2(dba.report(dbObj_cont,  contrast=1), "./results/diffbind_bnaive_mcl_dag.csv")

#creating a bed file with differentially enriched regions
export.bed(dba.report(dbObj_cont, contrast=1)[dba.report(dbObj_cont, contrast=1)$Fold > 0, ], con="results/up_naive_mcl.bed")
export.bed(dba.report(dbObj_cont, contrast=1)[dba.report(dbObj_cont, contrast=1)$Fold < 0, ], con="results/down_naive_mcl.bed")

down_mcl <- dba.report(dbObj_cont, contrast=1)[dba.report(dbObj_cont, contrast=1)$Fold > 0, ]
up_mcl <- dba.report(dbObj_cont, contrast=1)[dba.report(dbObj_cont, contrast=1)$Fold < 0, ]
# flipping the Fold, as the contrast happened to be naive vs mcl 

export.bed(up_mcl, con="results/h3k27ac_up_mcl_patiens.bed")
export.bed(up_mcl[up_mcl@seqnames@values == "chr11", ], con="results/h3k27ac_up_mcl_patiens_chr11.bed" )
export.bed(up_mcl[up_mcl@seqnames@values == "chr19", ], con="results/h3k27ac_up_mcl_patiens_chr19.bed" )
export.bed(up_mcl[up_mcl@seqnames@values == "chr17", ], con="results/h3k27ac_up_mcl_patiens_chr17.bed" )
export.bed(up_mcl[up_mcl@seqnames@values == "chr18", ], con="results/h3k27ac_up_mcl_patiens_chr18.bed" )
export.bed(up_mcl[up_mcl@seqnames@values == "chr22", ], con="results/h3k27ac_up_mcl_patiens_chr22.bed" )

export.bed(down_mcl, con="results/h3k27ac_down_mcl_patiens.bed")
export.bed(down_mcl[down_mcl@seqnames@values == "chr11", ], con="results/h3k27a_down_mcl_patiens_chr11.bed" )
export.bed(down_mcl[down_mcl@seqnames@values == "chr19", ], con="results/h3k27a_down_mcl_patiens_chr19.bed" )
export.bed(down_mcl[down_mcl@seqnames@values == "chr17", ], con="results/h3k27a_down_mcl_patiens_chr17.bed" )
export.bed(down_mcl[down_mcl@seqnames@values == "chr18", ], con="results/h3k27a_down_mcl_patiens_chr18.bed" )
export.bed(down_mcl[down_mcl@seqnames@values == "chr22", ], con="results/h3k27a_down_mcl_patiens_chr22.bed" )

#Affinity-based PCA with only the significant regions
dba.plotPCA(dbObj_cont, method=DBA_DESEQ2, attributes=DBA_CONDITION, label=DBA_ID)
dba.plotMA(dbObj_cont, contrast=1)

dba.plotBox(dbObj_cont, contrast=1)
dba.plotVolcano(dbObj_cont)

############################# Plotting profiles ################################
## the matrices were generated from the BED files with the differentially bound
## sites obtained at the previous step and the bigWig files from the pre-processing
## step using the computeMatrix function from deepTools

mcl_up.matrix <- import_deepToolsMat('matrix/mcl_up_both.matrix')
generateEnrichedHeatmap(mcl_up.matrix, matrices_color = c("white", "#8FBA31"), use_raster = T)

mcl_down.matrix <- import_deepToolsMat('matrix/mcl_down_both.matrix')
generateEnrichedHeatmap(mcl_down.matrix, matrices_color = c("white", "cornflowerblue"),  use_raster = T)

############################# Peak anotation ################################
covplot(up_mcl) #ChIP preaks over chromosomes
covplot(up_mcl, chrs = c("chr11", "chr14", "chr17", "chr18", "chr19", "chr22"))
covplot(down_mcl, chrs = c("chr11", "chr14", "chr17", "chr18", "chr19", "chr22"))
covplot(up_mcl, chrs = c("chr11", "chr14", "chr19"))

up_annot <- annotatePeak(up_mcl, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
down_annot <- annotatePeak(down_mcl, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
up_annot_chr19 <- annotatePeak(up_mcl[up_mcl@seqnames@values == "chr19", ], tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
up_annot_chr17 <- annotatePeak(up_mcl[up_mcl@seqnames@values == "chr17", ], tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
down_annot_chr19 <- annotatePeak(down_mcl[down_mcl@seqnames@values == "chr19", ], tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
down_annot_chr17 <- annotatePeak(down_mcl[down_mcl@seqnames@values == "chr17", ], tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
up_annot_chr1917 <- annotatePeak(up_mcl[up_mcl@seqnames@values == "chr19" | up_mcl@seqnames@values == "chr17" ], tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
down_annot_chr1917 <- annotatePeak(down_mcl[down_mcl@seqnames@values == "chr19" | down_mcl@seqnames@values == "chr17" ], tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")

plotAnnoPie(up_annot)
list <- list(up = up_annot, down = down_annot)
list19 <- list(up = up_annot_chr19, down = down_annot_chr19)
list17 <- list(up = up_annot_chr17, down = down_annot_chr17)
list1917 <- list(up = up_annot_chr1917, down = down_annot_chr1917)

plotAnnoBar(list)
plotAnnoBar(list17)
plotAnnoBar(list19)

plotDistToTSS(list)
plotDistToTSS(list19)

####################### Up peaks per chromosome ################################
chromsizes <- read.csv2("./chromSizes.csv", header = TRUE)
chromsizes$normsizes <- chromsizes$size/min(chromsizes$size)
chromsizes$normnumber <- chromsizes$gene_number/min(chromsizes$gene_number)

peaks_per_chr <- as.data.frame(table(up_mcl@seqnames))
peaks_per_chr <- peaks_per_chr[order(peaks_per_chr$Freq, decreasing = TRUE),]
peaks_per_chr <- peaks_per_chr[1:24, ]
colnames(peaks_per_chr) <- c("chr", "freq")
peaks_per_chr <- merge(peaks_per_chr,chromsizes, by.x = "chr", by.y = "chr")
peaks_per_chr$normfreq_size <- peaks_per_chr$freq/peaks_per_chr$normsizes #norm by size
peaks_per_chr$normfreq_number <- peaks_per_chr$freq/peaks_per_chr$normnumber #norm by gene number
write.csv(peaks_per_chr, "peaks_per_chr/up_peaks_per_chromosome_MCL_DAG.csv")

peaks_per_chr <- peaks_per_chr[order(peaks_per_chr$freq,decreasing=TRUE),]
barplot(peaks_per_chr$freq, names.arg = peaks_per_chr$chr, ylab = "N up peaks per chr")
peaks_per_chr <- peaks_per_chr[order(peaks_per_chr$normfreq_size,decreasing=TRUE),]
barplot(peaks_per_chr$normfreq_size, names.arg = peaks_per_chr$chr, ylab = "N up peaks normalised to chromSize")
peaks_per_chr <- peaks_per_chr[order(peaks_per_chr$normfreq_number,decreasing=TRUE),]
barplot(peaks_per_chr$normfreq_number, names.arg = peaks_per_chr$chr, ylab = "N up peaks normalised to gene number patients")

####################### Functional enrichment ################################
genes = lapply(list, function(i) as.data.frame(i)$geneId)
names(genes) = sub("_", "\n", names(genes))

genes19 = lapply(list19, function(i) as.data.frame(i)$geneId)
names(genes19) = sub("_", "\n", names(genes19))

genes17 = lapply(list17, function(i) as.data.frame(i)$geneId)
names(genes17) = sub("_", "\n", names(genes17))

genes1917 = lapply(list1917, function(i) as.data.frame(i)$geneId)
names(genes1917) = sub("_", "\n", names(genes1917))

universe <- read.csv("./MCL_cells_GOuniverse.csv")
compKEGG <- compareCluster(geneCluster   = genes,
                           fun           = "enrichKEGG",
                           pvalueCutoff  = 0.05,
                           universe = universe$ENTREZID,
                           pAdjustMethod = "BH")
write.csv2(compKEGG , "./results/kegg_enrichment_diffH3K27peaks_mcl_dag.csv")
dotplot(compKEGG, showCategory = 20, title = "KEGG Pathway Enrichment Analysis")

compGO <- compareCluster(geneCluster   = genes,
                           fun           = "enrichGO",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH",
                           universe = universe$ENTREZID,
                           OrgDb='org.Hs.eg.db')
write.csv2(compGO, "./results/go_enrichment_MF_diffH3K27peaks_mcl_dag.csv")
dotplot(compGO, showCategory = 20, title = "GO Enrichment Analysis MF")

groupGObp <- compareCluster(geneCluster   = genes,
                         fun           = "enrichGO",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         ont = "BP",
                         universe = universe$ENTREZID,
                         OrgDb='org.Hs.eg.db')
dotplot(groupGObp, showCategory = 10, title = "GO Enrichment Analysis BP")
write.csv2(groupGObp , "./results/go_enrichmentBP_diffH3K27peaks_mcl_dag.csv")

groupDO <- compareCluster(geneCluster   = genes,
                            fun           = "enrichDO",
                            pvalueCutoff  = 0.05,
                            universe = universe$ENTREZID,
                            pAdjustMethod = "BH")
dotplot(groupDO, showCategory = 20, title = "DO Enrichment Analysis")
write.csv2(groupDO, "./results/do_enrichment_diffH3K27peaks_mcl_dag.csv")

groupGObp19 <- compareCluster(geneCluster   = genes19,
                            fun           = "enrichGO",
                            pvalueCutoff  = 1,
                            universe = universe$ENTREZID,
                            pAdjustMethod = "BH",
                            ont = "BP",
                            OrgDb='org.Hs.eg.db')
dotplot(groupGObp19, showCategory = 9, title = "GO Enrichment Analysis BP chr19")
write.csv2(groupGObp19, "./results/go_enrichment_diffH3K27peaks_mcl_dag_chr19.csv")

groupDObp19 <- compareCluster(geneCluster   = genes19,
                              fun           = "enrichDO",
                              pvalueCutoff  = 1,
                              universe = universe$ENTREZID,
                              pAdjustMethod = "BH")
dotplot(groupDObp19, showCategory = 5, title = "DO Enrichment Analysis chr19")

groupGObp17 <- compareCluster(geneCluster   = genes17,
                              fun           = "enrichGO",
                              pvalueCutoff  = 1,
                              universe = universe$ENTREZID,
                              pAdjustMethod = "BH",
                              ont = "BP",
                              OrgDb='org.Hs.eg.db')
dotplot(groupGObp17, showCategory = 9, title = "GO Enrichment Analysis BP chr17")


groupGObp1917 <- compareCluster(geneCluster   = genes1917,
                              fun           = "enrichGO",
                              pvalueCutoff  = 0.05,
                              pAdjustMethod = "BH",
                              universe = universe$ENTREZID,
                              ont = "BP",
                              OrgDb='org.Hs.eg.db')
dotplot(groupGObp1917, showCategory = 10, title = "GO Enrichment Analysis BP chr19 + chr17")
write.csv2(groupGObp1917, "./results/go_enrichment_diffH3K27peaks_mcl_dag_chr1719.csv")


################################ Venn diagrams #################################
up_patients <- genes[[1]]
down_patients <- genes[[2]]
up_cells <- readRDS("~/Documents/Work/Projects/MCL/atac/diffbind_patients_dag/up_cells.rds")
down_cells <- readRDS("~/Documents/Work/Projects/MCL/atac/diffbind_patients_dag/down_cells.rds")

length(up_patients) #1248
length(down_patients) #3329
length(up_cells)
length(down_cells)

int <- list(A = up_cells, B = up_patients)
names(int) <- c("MCL cells vs c up","MCL patients vs c up")

v1 <- ggvenn(int, 
             fill_color = c("grey", 'grey'),
             stroke_size = 0.5, set_name_size = 4,
             show_percentage = F)

int2 <- list(A = down_cells, B = down_patients)
names(int2) <- c("MCL cells vs c down","MCL patients vs c down")

v2 <- ggvenn(int2, 
             fill_color = c("grey", 'grey'),
             stroke_size = 0.5, set_name_size = 4,
             show_percentage = F)

write.csv2(up_patients, "up_patients.csv")
write.csv2(up_cells, "up_cells.csv")

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
#  [1] TxDb.Hsapiens.UCSC.hg38.knownGene_3.17.0 clusterProfiler_4.8.2                   
#[3] ChIPseeker_1.36.0                        VennDiagram_1.7.3                       
#[5] futile.logger_1.4.3                      ggvenn_0.1.10                           
#[7] ggplot2_3.5.1                            dplyr_1.1.4                             
#[9] rtracklayer_1.60.1                       Cairo_1.6-2                             
#[11] nucleR_2.32.0                            csaw_1.33.0                             
#[13] EnsDb.Hsapiens.v86_2.99.0                ensembldb_2.24.1                        
#[15] AnnotationFilter_1.24.0                  GenomicFeatures_1.52.2                  
#[17] AnnotationDbi_1.62.2                     karyoploteR_1.26.0                      
#[19] regioneR_1.32.0                          patchwork_1.2.0                         
#[21] devtools_2.4.5                           usethis_2.2.3                           
#[23] DiffBind_3.10.1                          profileplyr_1.16.0                      
#[25] SummarizedExperiment_1.30.2              Biobase_2.60.0                          
#[27] GenomicRanges_1.52.1                     GenomeInfoDb_1.36.4                     
#[29] IRanges_2.34.1                           S4Vectors_0.38.2                        
#[31] MatrixGenerics_1.12.3                    matrixStats_1.3.0                       
#[33] BiocGenerics_0.46.0                     

#loaded via a namespace (and not attached):
#  [1] R.methodsS3_1.8.2                         dichromat_2.0-1                          
#[3] progress_1.2.3                            tiff_0.1-12                              
#[5] urlchecker_1.0.1                          nnet_7.3-19                              
#[7] DT_0.33                                   Biostrings_2.68.1                        
#[9] vctrs_0.6.5                               digest_0.6.35                            
#[11] png_0.1-8                                 shape_1.4.6.1                            
#[13] ggrepel_0.9.5                             mixsqp_0.3-54                            
#[15] org.Mm.eg.db_3.17.0                       deldir_2.0-4                             
#[17] magick_2.8.3                              MASS_7.3-60                              
#[19] reshape2_1.4.4                            SQUAREM_2021.1                           
#[21] httpuv_1.6.15                             foreach_1.5.2                            
#[23] qvalue_2.32.0                             withr_3.0.0                              
#[25] xfun_0.44                                 amap_0.8-19                              
#[27] ggfun_0.1.4                               ellipsis_0.3.2                           
#[29] memoise_2.0.1                             gson_0.1.0                               
#[31] profvis_0.3.8                             bamsignals_1.32.0                        
#[33] tidytree_0.4.6                            GlobalOptions_0.1.2                      
#[35] gtools_3.9.5                              R.oo_1.26.0                              
#[37] Formula_1.2-6                             prettyunits_1.2.0                        
#[39] KEGGREST_1.40.1                           promises_1.3.0                           
#[41] httr_1.4.7                                GreyListChIP_1.32.1                      
#[43] downloader_0.4                            restfulr_0.0.15                          
#[45] ashr_2.2-63                               rstudioapi_0.16.0                        
#[47] miniUI_0.1.1.1                            generics_0.1.3                           
#[49] DOSE_3.26.2                               base64enc_0.1-3                          
#[51] curl_5.2.1                                zlibbioc_1.46.0                          
#[53] EnrichedHeatmap_1.30.0                    ggraph_2.2.1                             
#[55] polyclip_1.10-6                           GenomeInfoDbData_1.2.10                  
#[57] xtable_1.8-4                              stringr_1.5.1                            
#[59] doParallel_1.0.17                         evaluate_0.23                            
#[61] S4Arrays_1.0.6                            systemPipeR_2.6.3                        
#[63] BiocFileCache_2.8.0                       preprocessCore_1.62.1                    
#[65] hms_1.1.3                                 irlba_2.3.5.1                            
#[67] colorspace_2.1-1                          filelock_1.0.3                           
#[69] magrittr_2.0.3                            later_1.3.2                              
#[71] viridis_0.6.5                             ggtree_3.8.2                             
#[73] lattice_0.22-6                            XML_3.99-0.16.1                          
#[75] shadowtext_0.1.3                          cowplot_1.1.3                            
#[77] Hmisc_5.1-2                               pillar_1.9.0                             
#[79] nlme_3.1-164                              iterators_1.0.14                         
#[81] caTools_1.18.2                            compiler_4.3.0                           
#[83] stringi_1.8.4                             GenomicAlignments_1.36.0                 
#[85] plyr_1.8.9                                crayon_1.5.2                             
#[87] abind_1.4-7                               BiocIO_1.10.0                            
#[89] truncnorm_1.0-9                           gridGraphics_0.5-1                       
#[91] emdbook_1.3.13                            locfit_1.5-9.9                           
#[93] graphlayouts_1.1.1                        org.Hs.eg.db_3.17.0                      
#[95] bit_4.0.5                                 fastmatch_1.1-4                          
#[97] codetools_0.2-20                          TxDb.Mmusculus.UCSC.mm10.knownGene_3.10.0
#[99] TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2   biovizBase_1.48.0                        
#[101] GetoptLong_1.0.5                          mime_0.12                                
#[103] splines_4.3.0                             circlize_0.4.16                          
#[105] Rcpp_1.0.12                               dbplyr_2.5.0                             
#[107] HDO.db_0.99.1                             interp_1.1-6                             
#[109] knitr_1.45                                blob_1.2.4                               
#[111] utf8_1.2.4                                clue_0.3-65                              
#[113] apeglm_1.22.1                             chipseq_1.50.0                           
#[115] fs_1.6.4                                  checkmate_2.3.1                          
#[117] pkgbuild_1.4.4                            ggplotify_0.1.2                          
#[119] rGREAT_2.2.0                              tibble_3.2.1                             
#[121] Matrix_1.6-5                              tweenr_2.0.3                             
#[123] pkgconfig_2.0.3                           pheatmap_1.0.12                          
#[125] tools_4.3.0                               cachem_1.1.0                             
#[127] RSQLite_2.3.7                             viridisLite_0.4.2                        
#[129] DBI_1.2.2                                 numDeriv_2022.9-1                        
#[131] fastmap_1.2.0                             rmarkdown_2.27                           
#[133] scales_1.3.0                              Rsamtools_2.16.0                         
#[135] coda_0.19-4.1                             VariantAnnotation_1.46.0                 
#[137] rpart_4.1.23                              farver_2.1.2                             
#[139] tidygraph_1.3.1                           scatterpie_0.2.2                         
#[141] yaml_2.3.8                                latticeExtra_0.6-30                      
#[143] foreign_0.8-86                            cli_3.6.2                                
#[145] purrr_1.0.2                               lifecycle_1.0.4                          
#[147] mvtnorm_1.2-5                             lambda.r_1.2.4                           
#[149] sessioninfo_1.2.2                         backports_1.5.0                          
#[151] BiocParallel_1.34.2                       gtable_0.3.5                             
#[153] rjson_0.2.21                              parallel_4.3.0                           
#[155] ape_5.8                                   limma_3.56.2                             
#[157] jsonlite_1.8.8                            edgeR_3.42.4                             
#[159] bitops_1.0-7                              bit64_4.0.5                              
#[161] yulab.utils_0.1.4                         soGGi_1.32.0                             
#[163] futile.options_1.0.1                      bezier_1.1.2                             
#[165] bdsmatrix_1.3-7                           metapod_1.7.0                            
#[167] GOSemSim_2.26.1                           R.utils_2.12.3                           
#[169] lazyeval_0.2.2                            shiny_1.8.1.1                            
#[171] htmltools_0.5.8.1                         enrichplot_1.20.0                        
#[173] GO.db_3.17.0                              formatR_1.14                             
#[175] rappdirs_0.3.3                            glue_1.7.0                               
#[177] TxDb.Mmusculus.UCSC.mm9.knownGene_3.2.2   XVector_0.40.0                           
#[179] RCurl_1.98-1.14                           treeio_1.24.3                            
#[181] BSgenome_1.68.0                           jpeg_0.1-10                              
#[183] gridExtra_2.3                             boot_1.3-30                              
#[185] igraph_2.0.3                              invgamma_1.1                             
#[187] R6_2.5.1                                  tidyr_1.3.1                              
#[189] DESeq2_1.40.2                             gplots_3.1.3.1                           
#[191] labeling_0.4.3                            cluster_2.1.6                            
#[193] bbmle_1.0.25.1                            pkgload_1.3.4                            
#[195] aplot_0.2.2                               DelayedArray_0.26.7                      
#[197] tidyselect_1.2.1                          plotrix_3.8-4                            
#[199] ProtGenerics_1.32.0                       htmlTable_2.4.2                          
#[201] ggforce_0.4.2                             xml2_1.3.6                               
#[203] munsell_0.5.1                             KernSmooth_2.23-24                       
#[205] data.table_1.15.4                         htmlwidgets_1.6.4                        
#[207] fgsea_1.26.0                              ComplexHeatmap_2.16.0                    
#[209] RColorBrewer_1.1-3                        hwriter_1.3.2.1                          
#[211] biomaRt_2.56.1                            rlang_1.1.3                              
#[213] remotes_2.5.0                             ShortRead_1.58.0                         
#[215] fansi_1.0.6          



