##################################
# Anna Schwager
# CNRS UMR9018, Institut Gustave Roussy 
# 2024 
##################################
library(ggvenn)
library(ggrepel)
library(ChIPpeakAnno)
library(plyr)
library(EnsDb.Hsapiens.v86)
library(biomaRt)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(clusterProfiler)
library(DESeq2)
library(ggplot2)
library(stringr)

################ Reading the data ####################
degs <- read.csv("rnaseq_data/resMCL_naive_blood.csv")
up_h3k27AC_peaks <- toGRanges("diffbind_data/h3k27ac_up_mcl_patiens.bed")
abc <- read.delim("enhancer_predictions/MCL/EnhancerPredictions.txt")

## as one enhancer mayregulate several genes, there are duplicate enhancer rows
## need to remove them for the GRanges intersections
abc_unique <- abc %>% distinct(name, .keep_all = TRUE)
abc_unique_gr <- toGRanges(abc_unique)
remove(abc_unique)

################ Identifying the enhancers that are more active in MCL #########
## setting at least 50% overlap threshold ## 
up_enhancers <- findOverlapsOfPeaks(abc_unique_gr, up_h3k27AC_peaks, minoverlap = 0.5)
up_enhancers <- up_enhancers$mergedPeaks
head(up_enhancers)
length(up_enhancers)
#939

############### Retrieving the genes regulated by these enhancers ###########
up_enhancers_names <- up_enhancers$peakNames
up_enhancers_names <- up_enhancers_names@unlistData[grep('abc_', up_enhancers_names@unlistData)]
head(up_enhancers_names)
#[1] "abc_unique_gr__genic.chr21.14522713.14523882"
up_enhancers_names <- gsub('abc_unique_gr__', '', up_enhancers_names)
#[1] "genic.chr21.14522713.14523882" 

abc$name[1]
#[1] "intergenic|chr21:8989205-8989705"
abc$name <- gsub('[: | -]', '.', abc$name)
abc$name[1]
#[1] "intergenic.chr21.8989205.8989705"

df <- as.data.frame(up_enhancers_names)
colnames(df) <- "name"

abc_up <- inner_join(abc, df, by = "name")
#abc_up - table of gene-enhancer pairs for the enhancers which are more active in MCL

############### Retrieving gene names to match with the DEGs ###########
## adding gene symbols
geneIDs <- ensembldb::select(EnsDb.Hsapiens.v86,
                             keys = gsub('\\..*', '', abc_up$TargetGene),
                             keytype = "TXID",
                             columns = c("SYMBOL","TXID"))

geneIDs <- geneIDs %>% distinct(SYMBOL, .keep_all = TRUE)
enhancer_associated_genes <- degs[degs$Gene_Name %in% geneIDs$SYMBOL, ]$Gene_Name
row.names(degs) <- degs$X
degs <- degs[,-1]

########################## Volcano plot ##########################

volc <- ggplot(data=degs) +
        geom_point(aes(x=log2FoldChange, y=-log10(padj)), col="grey") +
        geom_point(data=degs[degs$log2FoldChange < -1 & degs$padj < 0.05, ],
                   aes(x=log2FoldChange, y=-log10(padj)), col="skyblue4") +
        geom_point(data=degs[degs$log2FoldChange > 1 & degs$padj < 0.05, ],
                   aes(x=log2FoldChange, y=-log10(padj)), col="indianred") +
        geom_point(data=degs[degs$Gene_Name %in% geneIDs$SYMBOL & degs$padj < 0.05 & degs$log2FoldChange > 1,],
                   aes(x=log2FoldChange, y=-log10(padj)), col="#8FBA31") +
        theme_minimal() +
        theme(plot.title = element_text(size = rel(1.5), hjust = 0.5)) +
        ggtitle("MCL patients vs control naive B cells")

degs <- degs[order(degs$log2FoldChange,decreasing=TRUE),]

volc + geom_text_repel(data = head(degs[degs$Gene_Name %in% geneIDs$SYMBOL & degs$padj < 0.001,], 10), 
                aes(label = Gene_Name, x = log2FoldChange, y = -log10(padj)), box.padding = unit(.7, "lines"), hjust= 0.30)


write.csv(degs[degs$Gene_Name %in% geneIDs$SYMBOL & degs$padj < 0.05 & degs$log2FoldChange > 1,], "upregulated_genes_associated_with_upregulated_enhancers.csv")

##############################################################################
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
#  [1] stats4    grid      stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] ggrepel_0.9.5               purrr_1.0.2                 stringr_1.5.1              
#[4] DESeq2_1.40.2               SummarizedExperiment_1.30.2 MatrixGenerics_1.12.3      
#[7] matrixStats_1.3.0           clusterProfiler_4.8.2       enrichplot_1.20.0          
#[10] DOSE_3.26.2                 org.Hs.eg.db_3.17.0         biomaRt_2.56.1             
#[13] EnsDb.Hsapiens.v86_2.99.0   ensembldb_2.24.1            AnnotationFilter_1.24.0    
#[16] GenomicFeatures_1.52.2      AnnotationDbi_1.62.2        Biobase_2.60.0             
#[19] plyr_1.8.9                  ChIPpeakAnno_3.34.1         GenomicRanges_1.52.1       
#[22] GenomeInfoDb_1.36.4         IRanges_2.34.1              S4Vectors_0.38.2           
#[25] BiocGenerics_0.46.0         ggvenn_0.1.10               ggplot2_3.5.1              
#[28] dplyr_1.1.4                

#loaded via a namespace (and not attached):
#  [1] RColorBrewer_1.1-3       jsonlite_1.8.8           rstudioapi_0.16.0        magrittr_2.0.3          
#[5] farver_2.1.2             fs_1.6.4                 BiocIO_1.10.0            zlibbioc_1.46.0         
#[9] vctrs_0.6.5              multtest_2.56.0          memoise_2.0.1            Rsamtools_2.16.0        
#[13] RCurl_1.98-1.14          ggtree_3.8.2             S4Arrays_1.0.6           progress_1.2.3          
#[17] lambda.r_1.2.4           curl_5.2.1               gridGraphics_0.5-1       futile.options_1.0.1    
#[21] cachem_1.1.0             GenomicAlignments_1.36.0 igraph_2.0.3             lifecycle_1.0.4         
#[25] pkgconfig_2.0.3          gson_0.1.0               Matrix_1.6-5             R6_2.5.1                
#[29] fastmap_1.2.0            GenomeInfoDbData_1.2.10  aplot_0.2.2              digest_0.6.35           
#[33] colorspace_2.1-1         patchwork_1.2.0          regioneR_1.32.0          RSQLite_2.3.7           
#[37] labeling_0.4.3           filelock_1.0.3           fansi_1.0.6              polyclip_1.10-6         
#[41] httr_1.4.7               abind_1.4-7              compiler_4.3.0           downloader_0.4          
#[45] bit64_4.0.5              withr_3.0.0              BiocParallel_1.34.2      viridis_0.6.5           
#[49] DBI_1.2.2                ggforce_0.4.2            MASS_7.3-60              rappdirs_0.3.3          
#[53] DelayedArray_0.26.7      rjson_0.2.21             HDO.db_0.99.1            tools_4.3.0             
#[57] scatterpie_0.2.2         ape_5.8                  glue_1.7.0               VennDiagram_1.7.3       
#[61] restfulr_0.0.15          InteractionSet_1.28.1    nlme_3.1-164             GOSemSim_2.26.1         
#[65] shadowtext_0.1.3         reshape2_1.4.4           fgsea_1.26.0             generics_0.1.3          
#[69] gtable_0.3.5             BSgenome_1.68.0          tidyr_1.3.1              data.table_1.15.4       
#[73] hms_1.1.3                tidygraph_1.3.1          xml2_1.3.6               utf8_1.2.4              
#[77] XVector_0.40.0           pillar_1.9.0             yulab.utils_0.1.4        splines_4.3.0           
#[81] tweenr_2.0.3             treeio_1.24.3            BiocFileCache_2.8.0      lattice_0.22-6          
#[85] survival_3.6-4           rtracklayer_1.60.1       bit_4.0.5                tidyselect_1.2.1        
#[89] RBGL_1.76.0              locfit_1.5-9.9           GO.db_3.17.0             Biostrings_2.68.1       
#[93] gridExtra_2.3            ProtGenerics_1.32.0      futile.logger_1.4.3      graphlayouts_1.1.1      
#[97] stringi_1.8.4            ggfun_0.1.4              lazyeval_0.2.2           yaml_2.3.8              
#[101] codetools_0.2-20         ggraph_2.2.1             tibble_3.2.1             qvalue_2.32.0           
#[105] graph_1.78.0             ggplotify_0.1.2          cli_3.6.2                munsell_0.5.1           
#[109] Rcpp_1.0.12              dbplyr_2.5.0             png_0.1-8                XML_3.99-0.16.1         
#[113] parallel_4.3.0           blob_1.2.4               prettyunits_1.2.0        bitops_1.0-7            
#[117] tidytree_0.4.6           viridisLite_0.4.2        scales_1.3.0             crayon_1.5.2            
#[121] rlang_1.1.3              cowplot_1.1.3            fastmatch_1.1-4          KEGGREST_1.40.1         
#[125] formatR_1.14      
