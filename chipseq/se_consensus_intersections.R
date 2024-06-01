##################################
# Anna Schwager
# CNRS UMR9018, Institut Gustave Roussy 
# 2024 
##################################
library(ggvenn)
library(ChIPpeakAnno)
library(plyr)

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
sizes <- data.frame(condition = c(rep("mcl", length(dfm$width)), rep("naive", length(dfn$width)), rep("germinal", length(dfg$width))), size = c(dfm$width, dfn$width, dfg$width))
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
                labs(title="Sizes of SEs + median")

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

##########################################
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
#  [1] grid      stats4    stats     graphics  grDevices utils     datasets 
#[8] methods   base     

#other attached packages:
#  [1] plyr_1.8.9           ChIPpeakAnno_3.34.1  ggvenn_0.1.10       
#[4] ggplot2_3.5.1        dplyr_1.1.4          GenomicRanges_1.52.1
#[7] GenomeInfoDb_1.36.4  IRanges_2.34.1       S4Vectors_0.38.2    
#[10] BiocGenerics_0.46.0 

#loaded via a namespace (and not attached):
# [1] DBI_1.2.2                   bitops_1.0-7               
#[3] RBGL_1.76.0                 formatR_1.14               
#[5] biomaRt_2.56.1              rlang_1.1.3                
#[7] magrittr_2.0.3              matrixStats_1.3.0          
#[9] compiler_4.3.0              RSQLite_2.3.7              
#[11] GenomicFeatures_1.52.2      png_0.1-8                  
#[13] vctrs_0.6.5                 stringr_1.5.1              
#[15] pkgconfig_2.0.3             crayon_1.5.2               
#[17] fastmap_1.2.0               dbplyr_2.5.0               
#[19] XVector_0.40.0              labeling_0.4.3             
#[21] utf8_1.2.4                  Rsamtools_2.16.0           
#[23] graph_1.78.0                bit_4.0.5                  
#[25] zlibbioc_1.46.0             cachem_1.1.0               
#[27] progress_1.2.3              blob_1.2.4                 
#[29] DelayedArray_0.26.7         BiocParallel_1.34.2        
#[31] parallel_4.3.0              prettyunits_1.2.0          
#[33] R6_2.5.1                    stringi_1.8.4              
#[35] rtracklayer_1.60.1          Rcpp_1.0.12                
#[37] SummarizedExperiment_1.30.2 VennDiagram_1.7.3          
#[39] Matrix_1.6-5                splines_4.3.0              
#[41] tidyselect_1.2.1            rstudioapi_0.16.0          
#[43] abind_1.4-7                 yaml_2.3.8                 
#[45] codetools_0.2-20            curl_5.2.1                 
#[49] regioneR_1.32.0             InteractionSet_1.28.1      
#[51] Biobase_2.60.0              withr_3.0.0                
#[53] KEGGREST_1.40.1             lambda.r_1.2.4             
#[55] survival_3.6-4              futile.logger_1.4.3        
#[57] BiocFileCache_2.8.0         xml2_1.3.6                 
#[59] Biostrings_2.68.1           pillar_1.9.0               
#[61] filelock_1.0.3              MatrixGenerics_1.12.3      
#[63] generics_0.1.3              RCurl_1.98-1.14            
#[65] hms_1.1.3                   munsell_0.5.1              
#[67] scales_1.3.0                glue_1.7.0                 
#[69] tools_4.3.0                 BiocIO_1.10.0              
#[71] BSgenome_1.68.0             GenomicAlignments_1.36.0   
#[73] XML_3.99-0.16.1             AnnotationDbi_1.62.2       
#[75] colorspace_2.1-1            GenomeInfoDbData_1.2.10    
#[77] restfulr_0.0.15             cli_3.6.2                  
#[79] rappdirs_0.3.3              futile.options_1.0.1       
#[81] fansi_1.0.6                 S4Arrays_1.0.6             
#[83] gtable_0.3.5                digest_0.6.35              
#[85] farver_2.1.2                rjson_0.2.21               
#[87] memoise_2.0.1               multtest_2.56.0            
#[89] lifecycle_1.0.4             httr_1.4.7                 
#[91] bit64_4.0.5                 MASS_7.3-60            