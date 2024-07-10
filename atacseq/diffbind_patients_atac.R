##################################
# Anna Schwager
# CNRS UMR9018, Institut Gustave Roussy 
# 2024 
##################################
library(DiffBind)
library(devtools)
library(patchwork)
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

################### Differential acessibility  analysis #########################
samples <- read.csv2('meta.csv')
dbObj <- dba(sampleSheet = samples)
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps = TRUE)

dba.plotPCA(dbObj, attributes = DBA_CONDITION, label = DBA_ID) #PCA plot using affinity (read count) data
plot(dbObj) #correlation plot using affinity (read count) data
dbObj_normalized <- dba.normalize(dbObj, background=TRUE) 
#normalisation using the default library-size based method (lib) + background normalisation as recommended for ATAC

### establishing contrast
dbObj_cont <- dba.contrast(dbObj_normalized, design="~Condition", minMembers = 2)
### filtering for blacklist
dbObj_cont <- dba.blacklist(dbObj_cont, greylist=FALSE)
#Removed 157 (of 59569) consensus peaks.

### differential accessibility analysis
dbObj_cont <- dba.analyze(dbObj_cont, bGreylist=FALSE, design="~Condition")
dba.show(dbObj_cont, bContrasts = T) #show summary

### affinity-based PCA with only the significant regions
dba.plotPCA(dbObj_cont, method=DBA_DESEQ2, attributes=DBA_CONDITION, label=DBA_ID)

### differential acessibility plots
dba.plotMA(dbObj_cont, contrast=1)
dba.plotBox(dbObj_cont, contrast=1)
dba.plotVolcano(dbObj_cont)

###saving results
dba.report(dbObj_cont,  contrast=1)
write.csv2(dba.report(dbObj_cont,  contrast=1), "./results/diffbind_bnaive_mcl_dag.csv")

######### Creating bed files with differentially accessible regions #############
export.bed(dba.report(dbObj_cont, contrast=1)[dba.report(dbObj_cont, contrast=1)$Fold > 0, ], con="results/up_c_mcl.bed")
export.bed(dba.report(dbObj_cont, contrast=1)[dba.report(dbObj_cont, contrast=1)$Fold < 0, ], con="results/down_c_mcl.bed")

down_mcl <- dba.report(dbObj_cont, contrast=1)[dba.report(dbObj_cont, contrast=1)$Fold > 0, ]
up_mcl <- dba.report(dbObj_cont, contrast=1)[dba.report(dbObj_cont, contrast=1)$Fold < 0, ]

export.bed(up_mcl[up_mcl@seqnames@values == "chr11", ], con="results/atac_up_mcl_patiens_chr11.bed" )
export.bed(up_mcl[up_mcl@seqnames@values == "chr14", ], con="results/atac_up_mcl_patiens_chr14.bed" )
export.bed(up_mcl[up_mcl@seqnames@values == "chr19", ], con="results/atac_up_mcl_patiens_chr19.bed" )

export.bed(down_mcl[down_mcl@seqnames@values == "chr11", ], con="results/atac_down_mcl_patiens_chr11.bed" )
export.bed(down_mcl[down_mcl@seqnames@values == "chr14", ], con="results/atac_down_mcl_patiens_chr11.bed" )
export.bed(down_mcl[down_mcl@seqnames@values == "chr19", ], con="results/atac_down_mcl_patiens_chr19.bed" )

########################### Plotting profiles ##################################
## the matrices were generated from the BED files with the differentially acessible
## sites obtained at the previous step and the bigWig files from the pre-processing
## step using the computeMatrix function from deepTools

mcl_down.matrix <- import_deepToolsMat('matrix/mcl_down.matrix')
generateEnrichedHeatmap(mcl_down.matrix, matrices_color = c("white", "cornflowerblue"))

mcl_up.matrix <- import_deepToolsMat('matrix/mcl_up.matrix')
generateEnrichedHeatmap(mcl_up.matrix, matrices_color = c("white", "indianred"), use_raster = T)

############################### Peak annotation ###############################

### peak locations on chromosomes
covplot(up_mcl, chrs = c("chr11", "chr14", "chr17", "chr18", "chr19", "chr22"))
covplot(down_mcl, chrs = c("chr11", "chr14", "chr17", "chr18", "chr19", "chr22"))
### peak locations relative to genomic features
up_annot <- annotatePeak(up_mcl, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
down_annot <- annotatePeak(down_mcl, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
plotAnnoPie(up_annot)
plotAnnoPie(down_annot)

list <- list(up = up_annot, down = down_annot)
plotAnnoBar(list)
plotDistToTSS(list)

## saving the genes associated with the up and down peaks for further
## intersection with the cell line data
genes = lapply(list, function(i) as.data.frame(i)$geneId)
names(genes) = sub("_", "\n", names(genes))

up_patients <- genes[[1]]
down_patients <- genes[[2]]
saveRDS(up_patients, "up_atac_assigned_genes_patients.rds")
saveRDS(down_patients, "down_atac_assigned_genes_patients.rds")

### functional annotation
compKEGG <- compareCluster(geneCluster   = genes,
                           fun           = "enrichKEGG",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH")
write.csv2(compKEGG , "./results/kegg_enrichment_diffatacpeaks_mcl_dag.csv")
dotplot(compKEGG, showCategory = 20, title = "KEGG Pathway Enrichment Analysis")

compGO <- compareCluster(geneCluster   = genes,
                         fun           = "enrichGO",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         OrgDb='org.Hs.eg.db')
write.csv2(compGO, "./results/go_enrichment_MF_diffatacpeaks_mcl_dag.csv")
dotplot(compGO, showCategory = 20, title = "GO Enrichment Analysis MF")

groupGObp <- compareCluster(geneCluster   = genes,
                            fun           = "enrichGO",
                            pvalueCutoff  = 0.05,
                            pAdjustMethod = "BH",
                            ont = "BP",
                            OrgDb='org.Hs.eg.db')
write.csv2(groupGObp , "./results/go_enrichmentBP_diffatacpeaks_mcl_dag.csv")
dotplot(groupGObp, showCategory = 10, title = "GO Enrichment Analysis BP")

groupDO <- compareCluster(geneCluster   = genes,
                          fun           = "enrichDO",
                          pvalueCutoff  = 0.05,
                          pAdjustMethod = "BH")
dotplot(groupDO, showCategory = 20, title = "DO Enrichment Analysis")
write.csv2(groupDO, "./results/do_enrichment_diffatacpeaks_mcl_dag.csv")

################################################################################
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
#[13] EnsDb.Hsapiens.v86_2.99.0                patchwork_1.2.0                         
#[15] devtools_2.4.5                           usethis_2.2.3                           
#[17] DiffBind_3.10.1                          profileplyr_1.16.0                      
#[19] SummarizedExperiment_1.30.2              MatrixGenerics_1.12.3                   
#[21] matrixStats_1.3.0                        ensembldb_2.24.1                        
#[23] AnnotationFilter_1.24.0                  GenomicFeatures_1.52.2                  
#[25] AnnotationDbi_1.62.2                     Biobase_2.60.0                          
#[27] GenomicRanges_1.52.1                     GenomeInfoDb_1.36.4                     
#[29] IRanges_2.34.1                           S4Vectors_0.38.2                        
#[31] BiocGenerics_0.46.0                     

#loaded via a namespace (and not attached):
#  [1] fs_1.6.4                                  ProtGenerics_1.32.0                      
#[3] bitops_1.0-7                              enrichplot_1.20.0                        
#[5] HDO.db_0.99.1                             httr_1.4.7                               
#[7] RColorBrewer_1.1-3                        doParallel_1.0.17                        
#[9] numDeriv_2022.9-1                         profvis_0.3.8                            
#[11] tools_4.3.0                               utf8_1.2.4                               
#[13] R6_2.5.1                                  DT_0.33                                  
#[15] lazyeval_0.2.2                            GetoptLong_1.0.5                         
#[17] apeglm_1.22.1                             urlchecker_1.0.1                         
#[19] withr_3.0.0                               prettyunits_1.2.0                        
#[21] gridExtra_2.3                             preprocessCore_1.62.1                    
#[23] cli_3.6.2                                 formatR_1.14                             
#[25] scatterpie_0.2.2                          labeling_0.4.3                           
#[27] SQUAREM_2021.1                            mvtnorm_1.2-5                            
#[29] mixsqp_0.3-54                             Rsamtools_2.16.0                         
#[31] yulab.utils_0.1.4                         gson_0.1.0                               
#[33] R.utils_2.12.3                            DOSE_3.26.2                              
#[35] sessioninfo_1.2.2                         plotrix_3.8-4                            
#[37] BSgenome_1.68.0                           invgamma_1.1                             
#[39] bbmle_1.0.25.1                            limma_3.56.2                             
#[41] rstudioapi_0.16.0                         RSQLite_2.3.7                            
#[43] generics_0.1.3                            gridGraphics_0.5-1                       
#[45] ggVennDiagram_1.5.2                       TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2  
#[47] shape_1.4.6.1                             BiocIO_1.10.0                            
#[49] hwriter_1.3.2.1                           gtools_3.9.5                             
#[51] GO.db_3.17.0                              Matrix_1.6-5                             
#[53] interp_1.1-6                              fansi_1.0.6                              
#[55] abind_1.4-7                               R.methodsS3_1.8.2                        
#[57] lifecycle_1.0.4                           edgeR_3.42.4                             
#[59] yaml_2.3.8                                gplots_3.1.3.1                           
#[61] qvalue_2.32.0                             BiocFileCache_2.8.0                      
#[63] blob_1.2.4                                promises_1.3.0                           
#[65] crayon_1.5.2                              bdsmatrix_1.3-7                          
#[67] miniUI_0.1.1.1                            lattice_0.22-6                           
#[69] cowplot_1.1.3                             KEGGREST_1.40.1                          
#[71] magick_2.8.3                              metapod_1.7.0                            
#[73] ComplexHeatmap_2.16.0                     pillar_1.9.0                             
#[75] soGGi_1.32.0                              fgsea_1.26.0                             
#[77] rjson_0.2.21                              boot_1.3-30                              
#[79] systemPipeR_2.6.3                         codetools_0.2-20                         
#[81] fastmatch_1.1-4                           glue_1.7.0                               
#[83] ShortRead_1.58.0                          downloader_0.4                           
#[85] ggfun_0.1.4                               GreyListChIP_1.32.1                      
#[87] remotes_2.5.0                             data.table_1.15.4                        
#[89] vctrs_0.6.5                               png_0.1-8                                
#[91] treeio_1.24.3                             org.Mm.eg.db_3.17.0                      
#[93] gtable_0.3.5                              amap_0.8-19                              
#[95] chipseq_1.50.0                            emdbook_1.3.13                           
#[97] cachem_1.1.0                              mime_0.12                                
#[99] TxDb.Mmusculus.UCSC.mm9.knownGene_3.2.2   S4Arrays_1.0.6                           
#[101] tidygraph_1.3.1                           coda_0.19-4.1                            
#[103] pheatmap_1.0.12                           iterators_1.0.14                         
#[105] ellipsis_0.3.2                            rGREAT_2.2.0                             
#[107] nlme_3.1-164                              ggtree_3.8.2                             
#[109] bit64_4.0.5                               progress_1.2.3                           
#[111] filelock_1.0.3                            irlba_2.3.5.1                            
#[113] KernSmooth_2.23-24                        colorspace_2.1-1                         
#[115] DBI_1.2.2                                 DESeq2_1.40.2                            
#[117] tidyselect_1.2.1                          bit_4.0.5                                
#[119] compiler_4.3.0                            curl_5.2.1                               
#[121] xml2_1.3.6                                TxDb.Mmusculus.UCSC.mm10.knownGene_3.10.0
#[123] DelayedArray_0.26.7                       shadowtext_0.1.3                         
#[125] scales_1.3.0                              caTools_1.18.2                           
#[127] rappdirs_0.3.3                            tiff_0.1-12                              
#[129] stringr_1.5.1                             digest_0.6.35                            
#[131] XVector_0.40.0                            htmltools_0.5.8.1                        
#[133] pkgconfig_2.0.3                           jpeg_0.1-10                              
#[135] dbplyr_2.5.0                              fastmap_1.2.0                            
#[137] rlang_1.1.3                               GlobalOptions_0.1.2                      
#[139] htmlwidgets_1.6.4                         shiny_1.8.1.1                            
#[141] EnrichedHeatmap_1.30.0                    farver_2.1.2                             
#[143] jsonlite_1.8.8                            BiocParallel_1.34.2                      
#[145] R.oo_1.26.0                               GOSemSim_2.26.1                          
#[147] RCurl_1.98-1.14                           magrittr_2.0.3                           
#[149] GenomeInfoDbData_1.2.10                   ggplotify_0.1.2                          
#[151] munsell_0.5.1                             Rcpp_1.0.12                              
#[153] ape_5.8                                   viridis_0.6.5                            
#[155] stringi_1.8.4                             ggraph_2.2.1                             
#[157] zlibbioc_1.46.0                           MASS_7.3-60                              
#[159] pkgbuild_1.4.4                            plyr_1.8.9                               
#[161] org.Hs.eg.db_3.17.0                       parallel_4.3.0                           
#[163] ggrepel_0.9.5                             deldir_2.0-4                             
#[165] Biostrings_2.68.1                         graphlayouts_1.1.1                       
#[167] splines_4.3.0                             hms_1.1.3                                
#[169] circlize_0.4.16                           locfit_1.5-9.9                           
#[171] igraph_2.0.3                              pkgload_1.3.4                            
#[173] reshape2_1.4.4                            biomaRt_2.56.1                           
#[175] futile.options_1.0.1                      XML_3.99-0.16.1                          
#[177] latticeExtra_0.6-30                       lambda.r_1.2.4                           
#[179] foreach_1.5.2                             tweenr_2.0.3                             
#[181] httpuv_1.6.15                             tidyr_1.3.1                              
#[183] purrr_1.0.2                               polyclip_1.10-6                          
#[185] clue_0.3-65                               ashr_2.2-63                              
#[187] ggforce_0.4.2                             xtable_1.8-4                             
#[189] restfulr_0.0.15                           tidytree_0.4.6                           
#[191] later_1.3.2                               viridisLite_0.4.2                        
#[193] truncnorm_1.0-9                           tibble_3.2.1                             
#[195] aplot_0.2.2                               memoise_2.0.1                            
#[197] GenomicAlignments_1.36.0                  cluster_2.1.6        




