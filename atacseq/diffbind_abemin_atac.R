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
library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

################### Differential accessibility  analysis #########################
samples <- read.csv2('meta.csv')
dbObj <- dba(sampleSheet = samples)
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps = TRUE)
saveRDS(dbObj, "dbObj.RDS")

dba.plotPCA(dbObj, attributes = DBA_CONDITION, label = DBA_ID) #PCA plot using affinity (read count) data
dbObj_normalized <- dba.normalize(dbObj, background=TRUE) 
#normalisation using the default library-size based method (lib) + background no

### establishing contrast
dbObj_cont <- dba.contrast(dbObj_normalized, design="~Condition", minMembers = 2)
dbObj_cont <- dba.contrast(dbObj_cont,
                           group1 = dbObj_cont$masks$BLAS,
                           group2 = dbObj_cont$masks$GRANTA,
                           name1="BLAS", name2="GRANTA")

dbObj_cont <- dba.contrast(dbObj_cont,
                           group2 = dbObj_cont$masks$GRANTA,
                           group1 = dbObj_cont$masks$GRANTA_Min,
                           name2="GRANTA", name1="GRANTA_Min")

dbObj_cont <- dba.contrast(dbObj_cont,
                           group2 = dbObj_cont$masks$GRANTA,
                           group1 = dbObj_cont$masks$GRANTA_Abe,
                           name2="GRANTA", name1="GRANTA_Abe")

dbObj_cont <- dba.contrast(dbObj_cont,
                           group2 = dbObj_cont$masks$BLAS,
                           group1 = dbObj_cont$masks$BLAS_Abe,
                           name2="BLAS", name1="BLAS_Abe")

dbObj_cont <- dba.contrast(dbObj_cont,
                           group2 = dbObj_cont$masks$BLAS,
                           group1 = dbObj_cont$masks$BLAS_Min,
                           name2="BLAS", name1="BLAS_Min")

### filtering for blacklist
dbObj_cont <- dba.blacklist(dbObj_cont, greylist=FALSE)
#Removed 25 (of 128549) consensus peaks.

### differential binding analysis
dbObj_cont <- dba.analyze(dbObj_cont, bGreylist=FALSE)
dba.show(dbObj_cont, bContrasts = TRUE)
saveRDS(dbObj_cont, "dbObj_cont.RDS")
#dba.report(dbObj_cont,  contrast=1)

#saving results
write.csv2(dba.report(dbObj_cont,  contrast=1), "./results/diffbind_blas_granta.csv")
write.csv2(dba.report(dbObj_cont,  contrast=2), "./results/diffbind_grantamin_granta.csv")
write.csv2(dba.report(dbObj_cont,  contrast=3), "./results/diffbind_grantaabe_granta.csv")
write.csv2(dba.report(dbObj_cont,  contrast=4), "./results/diffbind_blasabe_blas.csv")
write.csv2(dba.report(dbObj_cont,  contrast=5), "./results/diffbind_blasmin_blas.csv")

######### Creating bed files with differentially accessible regions #############
up_min_granta <- dba.report(dbObj_cont, contrast=2)[dba.report(dbObj_cont, contrast=2)$Fold > 0, ]
down_min_granta <- dba.report(dbObj_cont, contrast=2)[dba.report(dbObj_cont, contrast=2)$Fold < 0, ]

up_abe_granta <- dba.report(dbObj_cont, contrast=3)[dba.report(dbObj_cont, contrast=3)$Fold > 0, ]
down_abe_granta <- dba.report(dbObj_cont, contrast=3)[dba.report(dbObj_cont, contrast=3)$Fold < 0, ]

up_min_blas <- dba.report(dbObj_cont, contrast=5)[dba.report(dbObj_cont, contrast=5)$Fold > 0, ]
down_min_blas <- dba.report(dbObj_cont, contrast=5)[dba.report(dbObj_cont, contrast=5)$Fold < 0, ]

up_abe_blas <- dba.report(dbObj_cont, contrast=4)[dba.report(dbObj_cont, contrast=4)$Fold > 0, ]
down_abe_blas <- dba.report(dbObj_cont, contrast=4)[dba.report(dbObj_cont, contrast=4)$Fold < 0, ]

up_granta <- dba.report(dbObj_cont, contrast=1)[dba.report(dbObj_cont, contrast=1)$Fold < 0, ]
down_granta <- dba.report(dbObj_cont, contrast=1)[dba.report(dbObj_cont, contrast=1)$Fold > 0, ]

export.bed(up_min_granta, con="results/up_min_granta.bed")
export.bed(down_min_granta, con="results/down_min_granta.bed")
export.bed(up_abe_granta, con="results/up_abe_granta.bed")
export.bed(down_abe_granta, con="results/down_abe_granta.bed")
export.bed(up_min_blas, con="results/up_min_blas.bed")
export.bed(down_min_blas, con="results/down_min_blas.bed")
export.bed(up_abe_blas, con="results/up_abe_blas.bed")
export.bed(down_abe_blas, con="results/down_abe_blas.bed")
export.bed(up_granta, con="results/up_granta.bed")
export.bed(down_granta, con="results/down_granta.bed")

########################### Plotting profiles ##################################
bigWigs <- list.files("./bigwig/mergedReplicate/", full.names = TRUE)

granta_up_profile <- BamBigwig_to_chipProfile(bigWigs,
                                              testRanges = "./results/up_granta.bed",
                                              style = "point",
                                              format = "bigwig",
                                              distanceAround = 2000)
granta_up_profile <- as_profileplyr(granta_up_profile)
generateEnrichedHeatmap(granta_up_profile, matrices_color = c("white", "indianred"))

granta_down_profile <- BamBigwig_to_chipProfile(bigWigs,
                                              testRanges = "./results/down_granta.bed",
                                              style = "point",
                                              format = "bigwig",
                                              distanceAround = 2000)
granta_down_profile <- as_profileplyr(granta_down_profile)
generateEnrichedHeatmap(granta_down_profile, matrices_color = c("white", "cornflowerblue"))

min_up_profile <- BamBigwig_to_chipProfile(bigWigs,
                                             testRanges = "./results/up_min_granta.bed",
                                             style = "point",
                                             format = "bigwig",
                                             distanceAround = 2000)
min_up_profile <- as_profileplyr(min_up_profile)
generateEnrichedHeatmap(min_up_profile, matrices_color = c("white", "indianred"))

min_down_profile <- BamBigwig_to_chipProfile(bigWigs,
                                          testRanges = "./results/down_min_granta.bed",
                                          style = "point",
                                          format = "bigwig",
                                          distanceAround = 2000)
min_down_profile <- as_profileplyr(min_down_profile)
generateEnrichedHeatmap(min_down_profile, matrices_color = c("white", "cornflowerblue"))

mcl_up_profile <- BamBigwig_to_chipProfile(bigWigs,
                                           testRanges = "./results/down_c_mcl.bed",
                                           style = "point",
                                           format = "bigwig",
                                           distanceAround = 2000)
mcl_up_profile <- as_profileplyr(mcl_up_profile)
generateEnrichedHeatmap(mcl_up_profile, matrices_color = c("white", "indianred"))

mcl_down_profile <- BamBigwig_to_chipProfile(bigWigs,
                                           testRanges = "./results/up_c_mcl.bed",
                                           style = "point",
                                           format = "bigwig",
                                           distanceAround = 2000)
mcl_down_profile <- as_profileplyr(mcl_down_profile)
generateEnrichedHeatmap(mcl_down_profile, matrices_color = c("white", "cornflowerblue"))

############################### Venn diagrams #################################
dba.show(dbObj_cont, bContrasts = TRUE)
par(mfrow = c(1, 2))

dba.plotVenn(dbObj_cont, contrast = c(1,2), bLoss=TRUE, bAll=FALSE, main = "up in MCL down by Min") 
dba.plotVenn(dbObj_cont, contrast = c(1,2), bGain=TRUE, bAll=FALSE, main = "down in MCL up by Min") 

dba.plotVenn(dbObj_cont, contrast = c(1,3), bLoss=TRUE, bAll=FALSE, main = "up in MCL down by Abe") 
dba.plotVenn(dbObj_cont, contrast = c(1,3), bGain=TRUE, bAll=FALSE, main = "down in MCL up by Abe") 

dba.plotVenn(dbObj_cont, contrast = c(2,5), bLoss=TRUE, bAll=FALSE, main = "down by Min")
dba.plotVenn(dbObj_cont, contrast = c(2,5), bGain=TRUE, bAll=FALSE, main = "up by Min") 

dba.plotVenn(dbObj_cont, contrast = c(3,4), bLoss=TRUE, bAll=FALSE, main = "down by Abe")
dba.plotVenn(dbObj_cont, contrast = c(3,4), bGain=TRUE, bAll=FALSE, main = "up by Abe") 

### loading GRanges from the MCL patients vs B naive comparisons
up_mcl_granges <- readRDS("./up_mcl_granges.RDS")
down_mcl_granges <- readRDS("./down_mcl_granges.RDS")

overlap <- findOverlapsOfPeaks(up_mcl_granges, up_granta, down_min_granta)
venn <- makeVennDiagram(overlap)

overlap2 <- findOverlapsOfPeaks(up_mcl_granges, up_granta, down_abe_granta)
venn2 <- makeVennDiagram(overlap2)

overlap3 <- findOverlapsOfPeaks(down_mcl_granges, down_granta, up_min_granta)
venn3 <- makeVennDiagram(overlap3)

overlap4 <- findOverlapsOfPeaks(down_mcl_granges, down_granta, up_abe_granta)
venn4 <- makeVennDiagram(overlap4)

############################### Peak annotation ###############################
### peak locations on chromosomes
covplot(up_granta, chrs = c("chr11", "chr14", "chr17", "chr18", "chr19", "chr22"))
covplot(down_granta, chrs = c("chr11", "chr14", "chr17", "chr18", "chr19", "chr22"))

### peak locations relative to genomic features
down_min_granta_annot <- annotatePeak(down_min_granta, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
down_min_blas_annot <- annotatePeak(down_min_blas, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
down_abe_granta_annot <- annotatePeak(down_abe_granta, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
down_abe_blas_annot <- annotatePeak(down_abe_blas, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")

up_min_granta_annot <- annotatePeak(up_min_granta, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
up_min_blas_annot <- annotatePeak(up_min_blas, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
up_abe_granta_annot <- annotatePeak(up_abe_granta, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
up_abe_blas_annot <- annotatePeak(up_abe_blas, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")

up_granta_annot <- annotatePeak(up_granta, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
down_granta_annot <- annotatePeak(down_granta, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")

### peak distributions relative to genomic features
list1 <- list(up = up_min_granta_annot, down = down_min_granta_annot)
plotAnnoBar(list1)
plotDistToTSS(list1)

list2 <- list(up = up_abe_granta_annot, down = down_abe_granta_annot)
plotAnnoBar(list2)

list3 <- list(up = up_min_blas_annot, down = down_min_blas_annot)
plotAnnoBar(list3)

list4 <- list(up = up_abe_blas_annot, down = down_abe_blas_annot)
plotAnnoBar(list4)

list5 <- list(up = up_granta_annot, down = down_granta_annot)
plotAnnoBar(list5)

### functional enrichment analysis 
universe <- read.csv("./MCL_cells_GOuniverse.csv")
#all the genes detected in a parallel RNA-seq experiment

genes1 = lapply(list1, function(i) as.data.frame(i)$geneId)
genes2 = lapply(list2, function(i) as.data.frame(i)$geneId)
genes5 = lapply(list5, function(i) as.data.frame(i)$geneId)

names(genes1) = sub("_", "\n", names(genes1))
names(genes2) = sub("_", "\n", names(genes2))
names(genes5) = sub("_", "\n", names(genes5))

### minnelide
compKEGG1 <- compareCluster(geneCluster   = genes1,
                           fun           = "enrichKEGG",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH",
                           universe = universe)
write.csv2(compKEGG , "./results/kegg_enrichment_diffatacpeaks_minnelide_in_granta.csv")
dotplot(compKEGG1, showCategory = 10, title = "KEGG Pathway Enrichment Analysis GRANTA+Min vs GRANTA")

compGO1 <- compareCluster(geneCluster   = genes1,
                         fun           = "enrichGO",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         OrgDb='org.Hs.eg.db',
                         universe = universe)
write.csv2(compGO, "./results/go_enrichment_MF_diffatacpeaks_minnelide_in_granta.csv")
dotplot(compGO1, showCategory = 10, title = "GO Enrichment Analysis MF GRANTA+Min vs GRANTA")

groupGObp1 <- compareCluster(geneCluster   = genes1,
                            fun           = "enrichGO",
                            pvalueCutoff  = 0.05,
                            pAdjustMethod = "BH",
                            ont = "BP",
                            OrgDb='org.Hs.eg.db',
                            universe = universe)
write.csv2(groupGObp1, "./results/go_enrichment_BP_diffatacpeaks_minnelide_in_granta.csv")
dotplot(groupGObp1, showCategory = 10, title = "GO Enrichment Analysis BP GRANTA+Min vs GRANTA")

groupDO1 <- compareCluster(geneCluster   = genes1,
                          fun           = "enrichDO",
                          pvalueCutoff  = 0.05,
                          pAdjustMethod = "BH",
                          universe = universe)
write.csv2(groupDO, "./results/do_enrichment_diffatacpeaks_minnelide_in_granta.csv")
dotplot(groupDO1, showCategory = 15, title = "DO Enrichment Analysis GRANTA+Min vs GRANTA")

### abelaciclib
compKEGG2 <- compareCluster(geneCluster   = genes2,
                            fun           = "enrichKEGG",
                            pvalueCutoff  = 0.05,
                            pAdjustMethod = "BH")
# no enrichment
compGO2 <- compareCluster(geneCluster   = genes2,
                          fun           = "enrichGO",
                          pvalueCutoff  = 0.05,
                          pAdjustMethod = "BH",
                          OrgDb='org.Hs.eg.db')
# no enrichment 
groupGObp2 <- compareCluster(geneCluster   = genes2,
                             fun           = "enrichGO",
                             pvalueCutoff  = 0.05,
                             pAdjustMethod = "BH",
                             ont = "BP",
                             OrgDb='org.Hs.eg.db')
write.csv2(groupGObp1, "./results/go_enrichment_BP_diffatacpeaks_abema_in_granta.csv")
dotplot(groupGObp2, showCategory = 10, title = "GO Enrichment Analysis BP GRANTA+Abe vs GRANTA")

groupDO2 <- compareCluster(geneCluster   = genes2,
                           fun           = "enrichDO",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH")
# no enrichment 

####################### Up peaks per chromosome ################################
chromsizes <- read.csv2("./chromSizes.csv", header = TRUE)
chromsizes$normsizes <- chromsizes$size/min(chromsizes$size)
chromsizes$normnumber <- chromsizes$gene_number/min(chromsizes$gene_number)

### granta up
peaks_per_chr_granta_up <- as.data.frame(table(up_granta@seqnames))
peaks_per_chr_granta_up <- peaks_per_chr_granta_up[order(peaks_per_chr_granta_up$Freq, decreasing = TRUE),]
peaks_per_chr_granta_up <- peaks_per_chr_granta_up[1:23, ]
colnames(peaks_per_chr_granta_up) <- c("chr", "freq")
peaks_per_chr_granta_up <- merge(peaks_per_chr_granta_up,chromsizes, by.x = "chr", by.y = "chr")
peaks_per_chr_granta_up$normfreq_size <- peaks_per_chr_granta_up$freq/peaks_per_chr_granta_up$normsizes #norm by size
peaks_per_chr_granta_up$normfreq_number <- peaks_per_chr_granta_up$freq/peaks_per_chr_granta_up$normnumber #norm by gene number
write.csv(peaks_per_chr_granta_up, "results/peaks_per_chr_granta_up_atac.csv")

par(mfrow = c(1, 3))
peaks_per_chr_granta_up <- peaks_per_chr_granta_up[order(peaks_per_chr_granta_up$freq,decreasing=TRUE),]
barplot(peaks_per_chr_granta_up$freq, names.arg = peaks_per_chr_granta_up$chr, ylab = "N up peaks per chr")
peaks_per_chr_granta_up <- peaks_per_chr_granta_up[order(peaks_per_chr_granta_up$normfreq_size,decreasing=TRUE),]
barplot(peaks_per_chr_granta_up$normfreq_size, names.arg = peaks_per_chr_granta_up$chr, ylab = "N up peaks normalised to chromSize")
peaks_per_chr_granta_up <- peaks_per_chr_granta_up[order(peaks_per_chr_granta_up$normfreq_number,decreasing=TRUE),]
barplot(peaks_per_chr_granta_up$normfreq_number, names.arg = peaks_per_chr_granta_up$chr, ylab = "N up peaks normalised to gene number patients")

### minnelide down
peaks_per_chr_min_down <- as.data.frame(table(down_min_granta@seqnames))
peaks_per_chr_min_down <- peaks_per_chr_min_down[order(peaks_per_chr_min_down$Freq, decreasing = TRUE),]
peaks_per_chr_min_down <- peaks_per_chr_min_down[1:23, ]
colnames(peaks_per_chr_min_down) <- c("chr", "freq")
peaks_per_chr_min_down <- merge(peaks_per_chr_min_down,chromsizes, by.x = "chr", by.y = "chr")
peaks_per_chr_min_down$normfreq_size <- peaks_per_chr_min_down$freq/peaks_per_chr_min_down$normsizes #norm by size
peaks_per_chr_min_down$normfreq_number <- peaks_per_chr_min_down$freq/peaks_per_chr_min_down$normnumber #norm by gene number
write.csv(peaks_per_chr_min_down, "results/peaks_per_chr_min_down_atac.csv")

par(mfrow = c(1, 3))
peaks_per_chr_min_down <- peaks_per_chr_min_down[order(peaks_per_chr_min_down$freq,decreasing=TRUE),]
barplot(peaks_per_chr_min_down$freq, names.arg = peaks_per_chr_min_down$chr, ylab = "N up peaks per chr")
peaks_per_chr_min_down <- peaks_per_chr_min_down[order(peaks_per_chr_min_down$normfreq_size,decreasing=TRUE),]
barplot(peaks_per_chr_min_down$normfreq_size, names.arg = peaks_per_chr_min_down$chr, ylab = "N up peaks normalised to chromSize")
peaks_per_chr_min_down <- peaks_per_chr_min_down[order(peaks_per_chr_min_down$normfreq_number,decreasing=TRUE),]
barplot(peaks_per_chr_min_down$normfreq_number, names.arg = peaks_per_chr_min_down$chr, ylab = "N up peaks normalised to gene number patients")

### minnelide up
peaks_per_chr_min_up <- as.data.frame(table(up_min_granta@seqnames))
peaks_per_chr_min_up <- peaks_per_chr_min_up[order(peaks_per_chr_min_up$Freq, decreasing = TRUE),]
peaks_per_chr_min_up <- peaks_per_chr_min_up[1:23, ]
colnames(peaks_per_chr_min_up) <- c("chr", "freq")
peaks_per_chr_min_up <- merge(peaks_per_chr_min_up,chromsizes, by.x = "chr", by.y = "chr")
peaks_per_chr_min_up$normfreq_size <- peaks_per_chr_min_up$freq/peaks_per_chr_min_up$normsizes #norm by size
peaks_per_chr_min_up$normfreq_number <- peaks_per_chr_min_up$freq/peaks_per_chr_min_up$normnumber #norm by gene number
write.csv(peaks_per_chr_min_up, "results/peaks_per_chr_min_up_atac.csv")

peaks_per_chr_min_up <- peaks_per_chr_min_up[order(peaks_per_chr_min_up$freq,decreasing=TRUE),]
barplot(peaks_per_chr_min_up$freq, names.arg = peaks_per_chr_min_up$chr, ylab = "N up peaks per chr")
peaks_per_chr_min_up <- peaks_per_chr_min_up[order(peaks_per_chr_min_up$normfreq_size,decreasing=TRUE),]
barplot(peaks_per_chr_min_up$normfreq_size, names.arg = peaks_per_chr_min_up$chr, ylab = "N up peaks normalised to chromSize")
peaks_per_chr_min_up <- peaks_per_chr_min_up[order(peaks_per_chr_min_up$normfreq_number,decreasing=TRUE),]
barplot(peaks_per_chr_min_up$normfreq_number, names.arg = peaks_per_chr_min_up$chr, ylab = "N up peaks normalised to gene number patients")

### abemaciclib down
peaks_per_chr_abe_down <- as.data.frame(table(down_abe_granta@seqnames))
peaks_per_chr_abe_down <- peaks_per_chr_abe_down[order(peaks_per_chr_abe_down$Freq, decreasing = TRUE),]
# too little peaks (1-8 peaks per chr), no peaks on chr19 at all, will not plot

### abemaciclib up
peaks_per_chr_abe_up <- as.data.frame(table(up_abe_granta@seqnames))
peaks_per_chr_abe_up <- peaks_per_chr_abe_up[order(peaks_per_chr_abe_up$Freq, decreasing = TRUE),]
peaks_per_chr_abe_up <- peaks_per_chr_abe_up[1:23, ]
colnames(peaks_per_chr_abe_up) <- c("chr", "freq")
peaks_per_chr_abe_up <- merge(peaks_per_chr_abe_up,chromsizes, by.x = "chr", by.y = "chr")
peaks_per_chr_abe_up$normfreq_size <- peaks_per_chr_abe_up$freq/peaks_per_chr_abe_up$normsizes #norm by size
peaks_per_chr_abe_up$normfreq_number <- peaks_per_chr_abe_up$freq/peaks_per_chr_abe_up$normnumber #norm by gene number
write.csv(peaks_per_chr_abe_up, "results/peaks_per_chr_abe_up_atac.csv")

peaks_per_chr_abe_up <- peaks_per_chr_abe_up[order(peaks_per_chr_abe_up$freq,decreasing=TRUE),]
barplot(peaks_per_chr_abe_up$freq, names.arg = peaks_per_chr_abe_up$chr, ylab = "N up peaks per chr")
peaks_per_chr_abe_up <- peaks_per_chr_abe_up[order(peaks_per_chr_abe_up$normfreq_size,decreasing=TRUE),]
barplot(peaks_per_chr_abe_up$normfreq_size, names.arg = peaks_per_chr_abe_up$chr, ylab = "N up peaks normalised to chromSize")
peaks_per_chr_abe_up <- peaks_per_chr_abe_up[order(peaks_per_chr_abe_up$normfreq_number,decreasing=TRUE),]
barplot(peaks_per_chr_abe_up$normfreq_number, names.arg = peaks_per_chr_abe_up$chr, ylab = "N up peaks normalised to gene number patients")

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
#  [1] ChIPpeakAnno_3.34.1                      org.Hs.eg.db_3.17.0                      TxDb.Hsapiens.UCSC.hg38.knownGene_3.17.0
#[4] clusterProfiler_4.8.2                    ChIPseeker_1.36.0                        VennDiagram_1.7.3                       
#[7] futile.logger_1.4.3                      ggvenn_0.1.10                            ggplot2_3.5.1                           
#[10] dplyr_1.1.4                              rtracklayer_1.60.1                       Cairo_1.6-2                             
#[13] nucleR_2.32.0                            csaw_1.33.0                              EnsDb.Hsapiens.v86_2.99.0               
#[16] ensembldb_2.24.1                         AnnotationFilter_1.24.0                  GenomicFeatures_1.52.2                  
#[19] AnnotationDbi_1.62.2                     patchwork_1.2.0                          devtools_2.4.5                          
#[22] usethis_2.2.3                            DiffBind_3.10.1                          profileplyr_1.16.0                      
#[25] SummarizedExperiment_1.30.2              Biobase_2.60.0                           GenomicRanges_1.52.1                    
#[28] GenomeInfoDb_1.36.4                      IRanges_2.34.1                           S4Vectors_0.38.2                        
#[31] MatrixGenerics_1.12.3                    matrixStats_1.3.0                        BiocGenerics_0.46.0                     

#loaded via a namespace (and not attached):
#  [1] ProtGenerics_1.32.0                       fs_1.6.4                                  bitops_1.0-7                             
#[4] enrichplot_1.20.0                         HDO.db_0.99.1                             httr_1.4.7                               
#[7] RColorBrewer_1.1-3                        doParallel_1.0.17                         InteractionSet_1.28.1                    
#[10] numDeriv_2022.9-1                         profvis_0.3.8                             tools_4.3.0                              
#[13] utf8_1.2.4                                R6_2.5.1                                  DT_0.33                                  
#[16] lazyeval_0.2.2                            apeglm_1.22.1                             GetoptLong_1.0.5                         
#[19] urlchecker_1.0.1                          withr_3.0.0                               prettyunits_1.2.0                        
#[22] gridExtra_2.3                             preprocessCore_1.62.1                     cli_3.6.2                                
#[25] formatR_1.14                              scatterpie_0.2.2                          labeling_0.4.3                           
#[28] SQUAREM_2021.1                            mvtnorm_1.2-5                             mixsqp_0.3-54                            
#[31] Rsamtools_2.16.0                          yulab.utils_0.1.4                         gson_0.1.0                               
#[34] R.utils_2.12.3                            DOSE_3.26.2                               sessioninfo_1.2.2                        
#[37] plotrix_3.8-4                             BSgenome_1.68.0                           invgamma_1.1                             
#[40] bbmle_1.0.25.1                            limma_3.56.2                              rstudioapi_0.16.0                        
#[43] RSQLite_2.3.7                             generics_0.1.3                            gridGraphics_0.5-1                       
#[46] TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2   shape_1.4.6.1                             BiocIO_1.10.0                            
#[49] hwriter_1.3.2.1                           gtools_3.9.5                              GO.db_3.17.0                             
#[52] Matrix_1.6-5                              interp_1.1-6                              fansi_1.0.6                              
#[55] abind_1.4-7                               R.methodsS3_1.8.2                         lifecycle_1.0.4                          
#[58] edgeR_3.42.4                              yaml_2.3.8                                gplots_3.1.3.1                           
#[61] qvalue_2.32.0                             BiocFileCache_2.8.0                       blob_1.2.4                               
#[64] promises_1.3.0                            crayon_1.5.2                              bdsmatrix_1.3-7                          
#[67] miniUI_0.1.1.1                            lattice_0.22-6                            cowplot_1.1.3                            
#[70] KEGGREST_1.40.1                           metapod_1.7.0                             ComplexHeatmap_2.16.0                    
#[73] pillar_1.9.0                              soGGi_1.32.0                              fgsea_1.26.0                             
#[76] rjson_0.2.21                              boot_1.3-30                               systemPipeR_2.6.3                        
#[79] codetools_0.2-20                          fastmatch_1.1-4                           glue_1.7.0                               
#[82] ShortRead_1.58.0                          downloader_0.4                            ggfun_0.1.4                              
#[85] GreyListChIP_1.32.1                       remotes_2.5.0                             data.table_1.15.4                        
#[88] vctrs_0.6.5                               png_0.1-8                                 treeio_1.24.3                            
#[91] org.Mm.eg.db_3.17.0                       gtable_0.3.5                              amap_0.8-19                              
#[94] chipseq_1.50.0                            emdbook_1.3.13                            cachem_1.1.0                             
#[97] TxDb.Mmusculus.UCSC.mm9.knownGene_3.2.2   S4Arrays_1.0.6                            mime_0.12                                
#[100] tidygraph_1.3.1                           coda_0.19-4.1                             survival_3.6-4                           
#[103] pheatmap_1.0.12                           iterators_1.0.14                          ellipsis_0.3.2                           
#[106] rGREAT_2.2.0                              nlme_3.1-164                              ggtree_3.8.2                             
#[109] bit64_4.0.5                               progress_1.2.3                            filelock_1.0.3                           
#[112] irlba_2.3.5.1                             KernSmooth_2.23-24                        colorspace_2.1-1                         
#[115] DBI_1.2.2                                 DESeq2_1.40.2                             tidyselect_1.2.1                         
#[118] bit_4.0.5                                 compiler_4.3.0                            curl_5.2.1                               
#[121] graph_1.78.0                              xml2_1.3.6                                TxDb.Mmusculus.UCSC.mm10.knownGene_3.10.0
#[124] DelayedArray_0.26.7                       shadowtext_0.1.3                          scales_1.3.0                             
#[127] caTools_1.18.2                            RBGL_1.76.0                               rappdirs_0.3.3                           
#[130] tiff_0.1-12                               stringr_1.5.1                             digest_0.6.35                            
#[133] XVector_0.40.0                            htmltools_0.5.8.1                         pkgconfig_2.0.3                          
#[136] jpeg_0.1-10                               regioneR_1.32.0                           dbplyr_2.5.0                             
#[139] fastmap_1.2.0                             rlang_1.1.3                               GlobalOptions_0.1.2                      
#[142] htmlwidgets_1.6.4                         shiny_1.8.1.1                             EnrichedHeatmap_1.30.0                   
#[145] farver_2.1.2                              jsonlite_1.8.8                            BiocParallel_1.34.2                      
#[148] R.oo_1.26.0                               GOSemSim_2.26.1                           RCurl_1.98-1.14                          
#[151] magrittr_2.0.3                            GenomeInfoDbData_1.2.10                   ggplotify_0.1.2                          
#[154] munsell_0.5.1                             Rcpp_1.0.12                               ape_5.8                                  
#[157] viridis_0.6.5                             stringi_1.8.4                             ggraph_2.2.1                             
#[160] zlibbioc_1.46.0                           MASS_7.3-60                               pkgbuild_1.4.4                           
#[163] plyr_1.8.9                                parallel_4.3.0                            ggrepel_0.9.5                            
#[166] deldir_2.0-4                              Biostrings_2.68.1                         graphlayouts_1.1.1                       
#[169] splines_4.3.0                             multtest_2.56.0                           hms_1.1.3                                
#[172] circlize_0.4.16                           locfit_1.5-9.9                            igraph_2.0.3                             
#[175] pkgload_1.3.4                             reshape2_1.4.4                            biomaRt_2.56.1                           
#[178] futile.options_1.0.1                      XML_3.99-0.16.1                           latticeExtra_0.6-30                      
#[181] lambda.r_1.2.4                            foreach_1.5.2                             tweenr_2.0.3                             
#[184] httpuv_1.6.15                             tidyr_1.3.1                               purrr_1.0.2                              
#[187] polyclip_1.10-6                           clue_0.3-65                               ashr_2.2-63                              
#[190] ggforce_0.4.2                             xtable_1.8-4                              restfulr_0.0.15                          
#[193] tidytree_0.4.6                            later_1.3.2                               viridisLite_0.4.2                        
#[196] truncnorm_1.0-9                           tibble_3.2.1                              aplot_0.2.2                              
#[199] memoise_2.0.1                             GenomicAlignments_1.36.0                  cluster_2.1.6     



