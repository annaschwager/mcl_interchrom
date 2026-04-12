##################################
# Anna Schwager
# CNRS UMR9018, Institut Gustave Roussy 
# 2024 
##################################
library(DESeq2)
library(karyoploteR)
library(ggplot2)
library(ggrepel)
library(ggvenn)
library(VennDiagram)
library(org.Hs.eg.db)
library(ggplotify)
library(ggupset)
library(ggimage)
library(gridExtra)
library(gtable)
library(latticeExtra)
library(karyoploteR)
library(EnsDb.Hsapiens.v86)
library(csaw)
library(profileplyr)
library(nucleR)
library(Cairo)
library(rtracklayer)
library(clusterProfiler)
library(tidyverse)
library(fgsea)
library(msigdbr)
library(dplyr)
library(ggpubr)
library(GenomicRanges)
library(AnnotationDbi)
library(scales)
library(zoo)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

################### Loading the count table ##############################
setwd("~/Documents/Work/Projects/MCL/rnaseq/DAG_MCL")
counts = read.table("counts.txt", header = T)
colnames(counts)
colnames(counts) <- gsub('X.home.anna.data1.Projects.MCL.rnaseq_dag.star_salmon.', '', colnames(counts))
colnames(counts) <- gsub('.markdup.sorted.bam', '', colnames(counts))
rownames(counts) <- counts$Geneid
counts <- counts[,7:30]
head(counts)
str(counts)

######################### Quality checks ################################
par(mfrow=c(3,3))
#Histograms of log transformed counts to see reads distribution
hist(log(as.numeric(counts[,2])+1), main = "germinal_1", xlab = "log(counts+1)")
hist(log(as.numeric(counts[,3])+1), main = "germinal_2", xlab = "log(counts+1)")
hist(log(as.numeric(counts[,4])+1), main = "germinal_3", xlab = "log(counts+1)")

hist(log(as.numeric(counts[,5])+1), main = "GRANTA_1", xlab = "log(counts+1)")
hist(log(as.numeric(counts[,6])+1), main = "GRANTA_2", xlab = "log(counts+1)")
hist(log(as.numeric(counts[,7])+1), main = "GRANTA_3", xlab = "log(counts+1)")

hist(log(as.numeric(counts[,8])+1), main = "MCL_DAG_1", xlab = "log(counts+1)")
hist(log(as.numeric(counts[,9])+1), main = "MCL_DAG_2", xlab = "log(counts+1)")
hist(log(as.numeric(counts[,10])+1), main = "MCL_DAG_3", xlab = "log(counts+1)")

hist(log(as.numeric(counts[,8])+1), main = "MCL_our_1", xlab = "log(counts+1)")
hist(log(as.numeric(counts[,9])+1), main = "MCL_our_2", xlab = "log(counts+1)")
hist(log(as.numeric(counts[,10])+1), main = "MCL_our_3", xlab = "log(counts+1)")

hist(log(as.numeric(counts[,8])+1), main = "naive_blood_1", xlab = "log(counts+1)")
hist(log(as.numeric(counts[,9])+1), main = "naive_blood_2", xlab = "log(counts+1)")
hist(log(as.numeric(counts[,10])+1), main = "naive_blood_3", xlab = "log(counts+1)")

hist(log(as.numeric(counts[,8])+1), main = "RPMI_1", xlab = "log(counts+1)")
hist(log(as.numeric(counts[,9])+1), main = "RPMI_2", xlab = "log(counts+1)")
hist(log(as.numeric(counts[,10])+1), main = "RPMI_3", xlab = "log(counts+1)")

#comparing library sizes
boxplot(log(counts + 1))

#filtering low counts
colSums(counts)
counts_per_gene <- rowSums(counts)
table(counts_per_gene > 0)
counts <- counts[which(counts_per_gene > 0), ]
dim(counts)
boxplot(log(counts + 1))

#exporting the counts table for downstream applications
counts_export <- write.csv2(counts, "counts.csv")

################### Differential expression analysis with DESeq2 ##############
## MCL patients from EGA and our patients are treated separately to assess batch 
cond <- factor(gsub("_[0-9]", "", colnames(counts)))

dds <- DESeqDataSetFromMatrix(counts, colData = DataFrame(cond), design = ~ cond)
colData(dds)
dds <- DESeq(dds)

rld <- rlog(dds)
plotPCA(rld, intgroup="cond") #coloring the samples by biological replicate
plotPCA(rld, intgroup="cond", returnData = T)
# our patients and patients from EGA cluster together; ok to join

resMCL_our_MCL_DAG <- results(dds, contrast = c("cond", "MCL_our", "MCL_DAG"))
resMCL_DAG_naive_blood <- results(dds, contrast = c("cond", "MCL_DAG", "naive_blood"))
resMCL_our_naive_blood <- results(dds, contrast = c("cond", "MCL_our", "naive_blood"))
resGRANTA_naive_blood <- results(dds, contrast = c("cond", "GRANTA", "naive_blood"))
resGRANTA_MCL_DAG <- results(dds, contrast = c("cond", "GRANTA", "MCL_DAG"))
resGRANTA_MCL_our <- results(dds, contrast = c("cond", "GRANTA", "MCL_our"))
resMCL_DAG_germinal <- results(dds, contrast = c("cond", "MCL_DAG", "germinal"))

summary(resMCL_our_MCL_DAG)
summary(resMCL_DAG_naive_blood)
summary(resMCL_our_naive_blood)
summary(resGRANTA_naive_blood)
summary(resGRANTA_MCL_DAG)
summary(resGRANTA_MCL_our)
summary(resMCL_DAG_germinal)

resMCL_DAG_germinal <- resMCL_DAG_germinal[complete.cases(resMCL_DAG_germinal),]  #remove any rows with NA (row counts and outliers)
resMCL_DAG_germinal <- resMCL_DAG_germinal[order(resMCL_DAG_germinal$padj),] #order by adjusted p-value
write.csv(resMCL_DAG_germinal, "results/resMCL_DAG_germinal.csv")

resMCL_our_MCL_DAG <- resMCL_our_MCL_DAG[complete.cases(resMCL_our_MCL_DAG),]  #remove any rows with NA (row counts and outliers)
resMCL_our_MCL_DAG <- resMCL_our_MCL_DAG[order(resMCL_our_MCL_DAG$padj),] #order by adjusted p-value
write.csv(resMCL_our_MCL_DAG, "results/resMCL_our_MCL_DAG.csv")

resMCL_DAG_naive_blood <- resMCL_DAG_naive_blood[complete.cases(resMCL_DAG_naive_blood),]  #remove any rows with NA (row counts and outliers)
resMCL_DAG_naive_blood <- resMCL_DAG_naive_blood[order(resMCL_DAG_naive_blood$padj),] #order by adjusted p-value
write.csv(resMCL_DAG_naive_blood, "results/resMCL_DAG_naive_blood.csv")

resMCL_our_naive_blood <- resMCL_our_naive_blood[complete.cases(resMCL_our_naive_blood),]  #remove any rows with NA (row counts and outliers)
resMCL_our_naive_blood <- resMCL_our_naive_blood[order(resMCL_our_naive_blood$padj),] #order by adjusted p-value
write.csv(resMCL_our_naive_blood, "results/resMCL_our_naive_blood.csv")

resGRANTA_naive_blood <- resGRANTA_naive_blood[complete.cases(resGRANTA_naive_blood),]  #remove any rows with NA (row counts and outliers)
resGRANTA_naive_blood <- resGRANTA_naive_blood[order(resGRANTA_naive_blood$padj),] #order by adjusted p-value
write.csv(resGRANTA_naive_blood, "results/resGRANTA_naive_blood.csv")

resGRANTA_MCL_DAG <- resGRANTA_MCL_DAG[complete.cases(resGRANTA_MCL_DAG),]  #remove any rows with NA (row counts and outliers)
resGRANTA_MCL_DAG <- resGRANTA_MCL_DAG[order(resGRANTA_MCL_DAG$padj),] #order by adjusted p-value
write.csv(resGRANTA_MCL_DAG, "results/resGRANTA_MCL_DAG.csv")

resGRANTA_MCL_our <- resGRANTA_MCL_our[complete.cases(resGRANTA_MCL_our),]  #remove any rows with NA (row counts and outliers)
resGRANTA_MCL_our <- resGRANTA_MCL_our[order(resGRANTA_MCL_our$padj),] #order by adjusted p-value
write.csv(resGRANTA_MCL_our, "results/resGRANTA_MCL_our.csv")

#filtering the results
resMCL_DAG_naive_blood_sign <- resMCL_DAG_naive_blood[resMCL_DAG_naive_blood[,'padj'] < 0.05, ]
resMCL_DAG_naive_blood_up <- resMCL_DAG_naive_blood_sign[resMCL_DAG_naive_blood_sign[,'log2FoldChange'] > 1, ]
resMCL_DAG_naive_blood_down <- resMCL_DAG_naive_blood_sign[resMCL_DAG_naive_blood_sign[,'log2FoldChange'] < -1, ]
summary(resMCL_DAG_naive_blood_up)
summary(resMCL_DAG_naive_blood_down)
write.csv(resMCL_DAG_naive_blood_up, "results/resMCL_DAG_naive_blood_up.csv")
write.csv(resMCL_DAG_naive_blood_down, "results/resMCL_DAG_naive_blood_down.csv")

resMCL_our_naive_blood_sign <- resMCL_our_naive_blood[resMCL_our_naive_blood[,'padj'] < 0.05, ]
resMCL_our_naive_blood_up <- resMCL_our_naive_blood_sign[resMCL_our_naive_blood_sign[,'log2FoldChange'] > 1, ]
resMCL_our_naive_blood_down <- resMCL_our_naive_blood_sign[resMCL_our_naive_blood_sign[,'log2FoldChange'] < -1, ]
summary(resMCL_our_naive_blood_up)
summary(resMCL_our_naive_blood_down)
write.csv(resMCL_our_naive_blood_up, "results/resMCL_our_naive_blood_up.csv")
write.csv(resMCL_our_naive_blood_down, "results/resMCL_our_naive_blood_down.csv")

resGRANTA_naive_blood_sign <- resGRANTA_naive_blood[resGRANTA_naive_blood[,'padj'] < 0.05, ]
resGRANTA_naive_blood_up <- resGRANTA_naive_blood_sign[resGRANTA_naive_blood_sign[,'log2FoldChange'] > 1, ]
resGRANTA_naive_blood_down <- resGRANTA_naive_blood_sign[resGRANTA_naive_blood_sign[,'log2FoldChange'] < -1, ]
summary(resGRANTA_naive_blood_up)
summary(resGRANTA_naive_blood_down)
write.csv(resGRANTA_naive_blood_up, "results/resGRANTA_naive_blood_up.csv")
write.csv(resGRANTA_naive_blood_down, "results/resGRANTA_naive_blood_down.csv")


################### Differential expression analysis with DESeq2 ##############
## MCL patients from EGA and our patients are treated together
cond <- factor(gsub("_[0-9]", "", colnames(counts)))
cond <- factor(gsub("MCL_DAG", "MCL", cond))
cond <- factor(gsub("MCL_our", "MCL", cond))

dds2 <- DESeqDataSetFromMatrix(counts, colData = DataFrame(cond), design = ~ cond)
colData(dds2)
dds2 <- DESeq(dds2)
resMCL_naive_blood <- results(dds2, contrast = c("cond", "MCL", "naive_blood"))
summary(resMCL_naive_blood)
resMCL_naive_blood <- resMCL_naive_blood[complete.cases(resMCL_naive_blood),]  #remove any rows with NA (row counts and outliers)
resMCL_naive_blood <- resMCL_naive_blood[order(resMCL_naive_blood$padj),] #order by adjusted p-value
write.csv(resMCL_naive_blood, "results/resMCL_naive_blood.csv")

resMCL_germinal <- results(dds2, contrast = c("cond", "MCL", "germinal"))
summary(resMCL_germinal)
resMCL_germinal <- resMCL_germinal[complete.cases(resMCL_germinal),]  #remove any rows with NA (row counts and outliers)
resMCL_germinal <- resMCL_germinal[order(resMCL_germinal$padj),] #order by adjusted p-value
write.csv(resMCL_germinal, "results/resMCL_germinal.csv")

#filtering the results
resMCL_naive_blood_sign <- resMCL_naive_blood[resMCL_naive_blood[,'padj'] < 0.05, ]
resMCL_naive_blood_up <- resMCL_naive_blood_sign[resMCL_naive_blood_sign[,'log2FoldChange'] > 1, ]
resMCL_naive_blood_down <- resMCL_naive_blood_sign[resMCL_naive_blood_sign[,'log2FoldChange'] < -1, ]
summary(resMCL_naive_blood_up)
summary(resMCL_naive_blood_down)
write.csv(resMCL_naive_blood_up, "results/resMCL_naive_blood_up.csv")
write.csv(resMCL_naive_blood_down, "results/resMCL_naive_blood_down.csv")

resMCL_germinal_sign <- resMCL_germinal[resMCL_germinal[,'padj'] < 0.05, ]
resMCL_germinal_up <- resMCL_germinal_sign[resMCL_germinal_sign[,'log2FoldChange'] > 1, ]
resMCL_germinal_down <- resMCL_germinal_sign[resMCL_germinal_sign[,'log2FoldChange'] < -1, ]
summary(resMCL_germinal_up)
summary(resMCL_germinal_down)
write.csv(resMCL_germinal_up, "results/resMCL_germinal_up.csv")
write.csv(resMCL_germinal_down, "results/resMCL_germinal_down.csv")

############################### Venn diagrams ################################
intersection_dag_ours_granta <- list(A=rownames(resMCL_DAG_naive_blood_up),
                                  B=rownames(resMCL_our_naive_blood_up),
                                  C=rownames(resGRANTA_naive_blood_up))

names(intersection_dag_ours_granta) <- c("MCL_DAG_vs_naive_up",
                                         "MCL_this_study_vs_naive_up",
                                         "GRANTA_vs_naive_up")

venn <- ggvenn(intersection_dag_ours_granta,
               stroke_size = 0.5, set_name_size = 4,
               show_percentage = FALSE)

intersection_dag_ours_granta2 <- list(A=rownames(resMCL_DAG_naive_blood_down),
                                     B=rownames(resMCL_our_naive_blood_down),
                                     C=rownames(resGRANTA_naive_blood_down))

names(intersection_dag_ours_granta2) <- c("MCL_DAG_vs_naive_down",
                                         "MCL_this_study_vs_naive_down",
                                         "GRANTA_vs_naive_down")

venn2 <- ggvenn(intersection_dag_ours_granta2,
               stroke_size = 0.5, set_name_size = 4,
               show_percentage = FALSE)

intersection_MCL_granta_up <- list(A=rownames(resMCL_naive_blood_up),
                                B=rownames(resGRANTA_naive_blood_up))

names(intersection_MCL_granta_up) <- c("MCL_vs_naive_up",
                                    "GRANTA_vs_naive_up")

venn3 <- ggvenn(intersection_MCL_granta_up,
               stroke_size = 0.5, set_name_size = 4,
               show_percentage = FALSE)

intersection_MCL_granta_down <- list(A=rownames(resMCL_naive_blood_down),
                                   B=rownames(resGRANTA_naive_blood_down))

names(intersection_MCL_granta_down) <- c("MCL_vs_naive_down",
                                       "GRANTA_vs_naive_down")

venn4 <- ggvenn(intersection_MCL_granta_down,
                stroke_size = 0.5, set_name_size = 4,
                show_percentage = FALSE)

#Intersection files
degs_intersect_up = intersect(rownames(resMCL_DAG_naive_blood_up), rownames(resMCL_our_naive_blood_up))
degs_intersect_up = intersect(degs_intersect_up, rownames(resGRANTA_naive_blood_up))
length(degs_intersect_up)
write.csv(resMCL_DAG_naive_blood_up[degs_intersect_up, ], "up_degs_intersect_DAG_our_GRANTA_padjofDAG.csv")
write.csv(resMCL_our_naive_blood_up[degs_intersect_up, ], "up_degs_intersect_DAG_our_GRANTA_padjofour.csv")
write.csv(resGRANTA_naive_blood_up[degs_intersect_up, ], "up_degs_intersect_DAG_our_GRANTA_padjofGRANTA.csv")

degs_intersect_down = intersect(rownames(resMCL_DAG_naive_blood_down), rownames(resMCL_our_naive_blood_down))
degs_intersect_down = intersect(degs_intersect_down, rownames(resGRANTA_naive_blood_down))
length(degs_intersect_down)
write.csv(resMCL_DAG_naive_blood_down[degs_intersect_down, ], "down_degs_intersect_DAG_our_GRANTA_padjofDAG.csv")
write.csv(resMCL_our_naive_blood_down[degs_intersect_down, ], "down_degs_intersect_DAG_our_GRANTA_padjofour.csv")
write.csv(resGRANTA_naive_blood_down[degs_intersect_down, ], "down_degs_intersect_DAG_our_GRANTA_padjofGRANTA.csv")

degs_intersect_up_broad = intersect(rownames(resMCL_naive_blood_up), rownames(resGRANTA_naive_blood_up))
length(degs_intersect_up_broad)
write.csv(resMCL_naive_blood_up[degs_intersect_up_broad, ], "up_degs_intersect_MCL_GRANTA_padjofMCL.csv")
write.csv(resGRANTA_naive_blood_up[degs_intersect_up_broad, ], "up_degs_intersect_MCL_GRANTA_padjofGRANTA.csv")

degs_intersect_down_broad = intersect(rownames(resMCL_naive_blood_down), rownames(resGRANTA_naive_blood_down))
length(degs_intersect_down_broad)
write.csv(resMCL_naive_blood_down[degs_intersect_down_broad, ], "down_degs_intersect_MCL_GRANTA_padjofMCL.csv")
write.csv(resGRANTA_naive_blood_down[degs_intersect_down_broad, ], "down_degs_intersect_MCL_GRANTA_padjofGRANTA.csv")


############## Plotting the genes mapped to chromosomes ########################
ensdb <- EnsDb.Hsapiens.v86

##retrieving the ranges of DEGs for resMCL_naive_blood
all.genes.ens <- genes(ensdb)
head(all.genes.ens)
seqlevelsStyle(all.genes.ens) <- "UCSC"  #we need to do that, because kp genome is from UCSC

head(resMCL_naive_blood)
resMCL_naive_blood$Gene_Name <- rownames(resMCL_naive_blood)
mcols(all.genes.ens) <- resMCL_naive_blood[all.genes.ens$gene_name, c("log2FoldChange", "stat", "padj", "Gene_Name")]
head(all.genes.ens)
keep <- !is.na(mcols(all.genes.ens)$log2FoldChange)
all.genes.ens.MCL <- all.genes.ens[keep,] #remove NAs and keep only DEGs
head(all.genes.ens.MCL)
all.genes.ens.MCL_all <- all.genes.ens.MCL
all.genes.ens.MCL <- all.genes.ens.MCL[all.genes.ens.MCL$padj < 0.05]
head(all.genes.ens.MCL)
write.csv(all.genes.ens.MCL, "MCL_vs_c_DAG.csv")
saveRDS(all.genes.ens.MCL, "all.genes.ens.MCL.RDS")


col.over <- "indianred" 
col.under <- "cornflowerblue" 
col <- rep(col.under, length(all.genes.ens.MCL))
col[all.genes.ens.MCL$log2FoldChange>1] <- col.over

all.genes.ens <- genes(ensdb)
head(all.genes.ens)
seqlevelsStyle(all.genes.ens) <- "UCSC"  #we need to do that, because kp genome is from UCSC
head(resGRANTA_naive_blood)
resGRANTA_naive_blood$Gene_Name <- rownames(resGRANTA_naive_blood)
mcols(all.genes.ens) <- resGRANTA_naive_blood[all.genes.ens$gene_name, c("log2FoldChange", "stat", "padj", "Gene_Name")]
head(all.genes.ens)
keep <- !is.na(mcols(all.genes.ens)$log2FoldChange)
all.genes.ens.GRANTA <- all.genes.ens[keep,] #remove NAs and keep only DEGs
head(all.genes.ens.GRANTA)
all.genes.ens.GRANTA_all <- all.genes.ens.GRANTA
all.genes.ens.GRANTA <- all.genes.ens.GRANTA[all.genes.ens.GRANTA$padj < 0.05]
all.genes.ens.GRANTA001 <- all.genes.ens.GRANTA[all.genes.ens.GRANTA$padj < 0.01]
head(all.genes.ens.GRANTA)
saveRDS(all.genes.ens.GRANTA, "all.genes.ens.GRANTA.RDS")

col.over <- "indianred" 
col.under <- "cornflowerblue" 
col <- rep(col.under, length(all.genes.ens.GRANTA))
col[all.genes.ens.GRANTA$log2FoldChange>1] <- col.over
col001 <- rep(col.under, length(all.genes.ens.GRANTA001))
col001[all.genes.ens.GRANTA001$log2FoldChange>1] <- col.over

### MCL_broad vs naive
kp <- plotKaryotype(genome="hg38", chromosomes = c("chr11", "chr14"), plot.type = 2)
max(all.genes.ens.MCL$log2FoldChange)
min(all.genes.ens.MCL$log2FoldChange)
kpAxis(kp, ymax = 10, ymin = -10, )
kpAddLabels(kp, labels = "log2FC", srt=90, pos = 1, label.margin = 0.07, cex=0.9)
kpAddLabels(kp, labels="MCL DAG and ours", r1=0.1, r0=0, data.panel = 1, side = "right", cex=0.9)
kp <- kpPlotDensity(kp, all.genes.ens, data.panel = 2)
kpPoints(kp, data = all.genes.ens.MCL, y=all.genes.ens.MCL$log2FoldChange,
         ymax = 10, ymin = -10, col=col)

kp <- plotKaryotype(genome="hg38", chromosomes = c("chr19"), plot.type = 2)
max(all.genes.ens.MCL$log2FoldChange)
min(all.genes.ens.MCL$log2FoldChange)
kpAxis(kp, ymax = 10, ymin = -10, )
kpAddLabels(kp, labels = "log2FC", srt=90, pos = 1, label.margin = 0.07, cex=0.9)
kpAddLabels(kp, labels="MCL DAG and ours vs B naive", r1=0.1, r0=0, data.panel = 1, side = "right", cex=0.9)
kp <- kpPlotDensity(kp, all.genes.ens, data.panel = 2)
kpPoints(kp, data = all.genes.ens.MCL, y=all.genes.ens.MCL$log2FoldChange,
         ymax = 10, ymin = -10, col=col)

### GRANTA vs naive
kp <- plotKaryotype(genome="hg38", chromosomes = c("chr11", "chr14"), plot.type = 2)
max(all.genes.ens.GRANTA$log2FoldChange)
min(all.genes.ens.GRANTA$log2FoldChange)
kpAxis(kp, ymax = 15, ymin = -15)
kpAddLabels(kp, labels = "log2FC", srt=90, pos = 1, label.margin = 0.07, cex=0.9)
kpAddLabels(kp, labels="GRANTA", r1=0.1, r0=0, data.panel = 1, side = "right", cex=0.9)
kp <- kpPlotDensity(kp, all.genes.ens, data.panel = 2)

kpPoints(kp, data = all.genes.ens.GRANTA, y=all.genes.ens.GRANTA$log2FoldChange,
         ymax = 15, ymin = -15, col=col)
kpPoints(kp, data = all.genes.ens.GRANTA001, y=all.genes.ens.GRANTA001$log2FoldChange,
         ymax = 15, ymin = -15, col=col001)

kp <- plotKaryotype(genome="hg38", chromosomes = c("chr19"), plot.type = 2)
max(all.genes.ens.GRANTA001$log2FoldChange)
min(all.genes.ens.GRANTA001$log2FoldChange)
kpAxis(kp, ymax = 10, ymin = -10, )
kpAddLabels(kp, labels = "log2FC", srt=90, pos = 1, label.margin = 0.07, cex=0.9)
kpAddLabels(kp, labels="GRANTA vs B naive", r1=0.1, r0=0, data.panel = 1, side = "right", cex=0.9)
kp <- kpPlotDensity(kp, all.genes.ens, data.panel = 2)
kpPoints(kp, data = all.genes.ens.GRANTA001, y=all.genes.ens.GRANTA001$log2FoldChange,
         ymax = 10, ymin = -10, col=col001)

############## Barplots DEGs per chromosome ########################
### DEGs per chromosome, MCL
chromsizes <- read.csv2("./chromSizes.csv", header = TRUE)
chromsizes$normsizes <- chromsizes$size/min(chromsizes$size)
chromsizes$normnumber <- chromsizes$gene_number/min(chromsizes$gene_number)

##  all DEGS
df_MCL_001 <- as.data.frame(table(as.vector(all.genes.ens.MCL@seqnames)))
df_MCL_001 <- df_MCL_001[153:176,]
colnames(df_MCL_001) <- c("chr", "freq")
df_MCL_001 <- merge(df_MCL_001,chromsizes, by.x = "chr", by.y = "chr")
df_MCL_001$normfreq_size <- df_MCL_001$freq/df_MCL_001$normsizes #norm by size
df_MCL_001$normfreq_number <- df_MCL_001$freq/df_MCL_001$normnumber #norm by gene number
write.csv(df_MCL_001, "results/degs_per_chromosome_MCL_padj001.csv")


## up DEGs
df_MCL_001_up <- as.data.frame(table(as.vector(all.genes.ens.MCL[all.genes.ens.MCL@elementMetadata@listData[["log2FoldChange"]] > 1, ]@seqnames)))
df_MCL_001_up <- df_MCL_001_up[117:139,]
colnames(df_MCL_001_up) <- c("chr", "freq")
df_MCL_001_up <- merge(df_MCL_001_up,chromsizes, by.x = "chr", by.y = "chr")
df_MCL_001_up$normfreq_size <- df_MCL_001_up$freq/df_MCL_001_up$normsizes #norm by size
df_MCL_001_up$normfreq_number <- df_MCL_001_up$freq/df_MCL_001_up$normnumber #norm by gene number
write.csv(df_MCL_001_up, "results/up_degs_per_chromosome_MCL_padj001.csv")

df_MCL_001_up <- df_MCL_001_up[order(df_MCL_001_up$freq,decreasing=TRUE),]
barplot(df_MCL_001_up$freq, names.arg = df_MCL_001_up$chr, ylab = "N up DEGs patients") #up DEGs raw
df_MCL_001_up <- df_MCL_001_up[order(df_MCL_001_up$normfreq_size,decreasing=TRUE),]
barplot(df_MCL_001_up$normfreq_size, names.arg = df_MCL_001_up$chr, ylab = "N up DEGs normalised to chromSize patients")
df_MCL_001_up <- df_MCL_001_up[order(df_MCL_001_up$normfreq_number,decreasing=TRUE),]
barplot(df_MCL_001_up$normfreq_number, names.arg = df_MCL_001_up$chr, ylab = "N up DEGs normalised to gene number patients")


### DEGs per chromosome, GRANTA
df_GRANTA_001 <- as.data.frame(table(as.vector(all.genes.ens.GRANTA001@seqnames)))
df_GRANTA_001 <- df_GRANTA_001[165:188,]
colnames(df_GRANTA_001) <- c("chr", "freq")
df_GRANTA_001 <- merge(df_GRANTA_001,chromsizes, by.x = "chr", by.y = "chr")
df_GRANTA_001$normfreq_size <- df_GRANTA_001$freq/df_GRANTA_001$normsizes #norm by size
df_GRANTA_001$normfreq_number <- df_GRANTA_001$freq/df_GRANTA_001$normnumber #norm by gene number
write.csv(df_MCL_001_up, "results/degs_per_chromosome_GRANTA_padj001.csv")

## up DEGs
df_GRANTA_001_up <- as.data.frame(table(as.vector(all.genes.ens.GRANTA001[all.genes.ens.GRANTA001@elementMetadata@listData[["log2FoldChange"]] > 1, ]@seqnames)))
df_GRANTA_001_up <- df_GRANTA_001_up[110:132,]
colnames(df_GRANTA_001_up) <- c("chr", "freq")
df_GRANTA_001_up <- merge(df_GRANTA_001_up,chromsizes, by.x = "chr", by.y = "chr")
df_GRANTA_001_up$normfreq_size <- df_GRANTA_001_up$freq/df_GRANTA_001_up$normsizes #norm by size
df_GRANTA_001_up$normfreq_number <- df_GRANTA_001_up$freq/df_GRANTA_001_up$normnumber #norm by gene number
write.csv(df_GRANTA_001_up, "results/up_degs_per_chromosome_GRANTA_padj001.csv")

df_GRANTA_001_up <- df_GRANTA_001_up[order(df_GRANTA_001_up$freq,decreasing=TRUE),]
barplot(df_GRANTA_001_up$freq, names.arg = df_GRANTA_001_up$chr, ylab = "N up DEGs granta") #up DEGs raw
df_GRANTA_001_up <- df_GRANTA_001_up[order(df_GRANTA_001_up$normfreq_size,decreasing=TRUE),]
barplot(df_GRANTA_001_up$normfreq_size, names.arg = df_GRANTA_001_up$chr, ylab = "N up DEGs normalised to chromSize granta")
df_GRANTA_001_up <- df_GRANTA_001_up[order(df_GRANTA_001_up$normfreq_number,decreasing=TRUE),]
barplot(df_GRANTA_001_up$normfreq_number, names.arg = df_GRANTA_001_up$chr, ylab = "N up DEGs normalised to gene number granta")

############## GO enrichment on chromosomes ########################
all.genes.ens.MCL.19 <- all.genes.ens.MCL[all.genes.ens.MCL@seqnames == "chr19", ]
all.genes.ens.MCL.19$Gene_Name

eg.19 = bitr(all.genes.ens.MCL.19$Gene_Name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(eg.19)
eg.universe = bitr(row.names(counts), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
write.csv(eg.universe, "results/MCL_cells_GOuniverse.csv")


ego <- enrichGO(gene          = eg.19$ENTREZID,
                universe      = eg.universe$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.1,
                qvalueCutoff  = 0.1,
                readable      = TRUE)
dotplot(ego, showCategory = 10, title = "GO Enrichment Analysis BP DEGs on chr19")

chr19_up = intersect(rownames(resMCL_naive_blood_up), eg.19$SYMBOL)
chr19_down = intersect(rownames(resMCL_naive_blood_down), eg.19$SYMBOL)
eg.19.up = bitr(chr19_up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

ego_up <- enrichGO(gene          = eg.19.up$ENTREZID,
                   universe      = eg.universe$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.1,
                   qvalueCutoff  = 0.1,
                   readable      = TRUE)
dotplot(ego_up, showCategory = 10, title = "GO Enrichment Analysis BP upregulated DEGS on chr19")

################### Exporting expression data for ABC #######################
MCL_DAG_counts <- counts[,7:11]
a <- c("r", "r", "R", "R", "R")

ddsm <- DESeqDataSetFromMatrix(MCL_DAG_counts, colData = DataFrame(a), design = ~ a)
colData(ddsm)
ddsm <- DESeq(ddsm)
normalized_counts_MCL = counts(ddsm, normalized=TRUE)
write(normalized_counts_MCL, "results/MCL_norm_counts.txt")

a <- as.data.frame(normalized_counts_MCL)
a$average <- sum(a$MCL_DAG_1, a$MCL_DAG_2, a$MCL_DAG_3, a$MCL_DAG_4, a$MCL_DAG_5)/5

normalized_counts_MCL$average <- sum(normalized_counts_MCL$MCL_)

naive_counts <- counts[,16:21]
condition <-factor(c("n", "n", "n", "N", "N", "N"))
condition <- DataFrame(condition)
row.names(condition) <- colnames(naive_counts)

ddsm_n <- DESeqDataSetFromMatrix(naive_counts, colData = condition, ~ condition)
ddsm_n <- DESeq(ddsm_n)
normalized_counts_naive = counts(ddsm_n, normalized=TRUE)
write.csv2(normalized_counts_naive, "results/naive_norm_counts.csv")

n <- as.data.frame(normalized_counts_naive)
n$average <- rowMeans(n)
naive_average <- data.frame(row.names(n), n$average)
write.csv2(naive_average, "results/naive_norm_counts_averagevalues.csv")


################# Observed vs expected N up genes #############################
### MCL
## Prepare per-chromosome table 
genes_df <- as.data.frame(all.genes.ens.MCL_all) %>%
  transmute(
    chr = as.character(seqnames),
    Gene_Name = Gene_Name,
    padj,
    log2FoldChange
  ) %>%
  filter(grepl("^chr([0-9]+|X|Y)$", chr)) %>%
  distinct(chr, Gene_Name, .keep_all = TRUE) 

df_panelA <- genes_df %>%
  group_by(chr) %>%
  summarise(
    n_genes_chr = n(),
    n_up_obs    = sum(padj < 0.05 & log2FoldChange > 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    total_genes = sum(n_genes_chr),
    total_up    = sum(n_up_obs),
    p0          = total_up / total_genes,
    n_up_expected = total_up * (n_genes_chr / total_genes),
    enrichment    = n_up_obs / n_up_expected,
    excess        = n_up_obs - n_up_expected
  ) %>%
  rowwise() %>%
  mutate(
    p_enrich = binom.test(
      x = as.integer(n_up_obs),
      n = as.integer(n_genes_chr),
      p = p0,
      alternative = "greater"
    )$p.value
  ) %>%
  ungroup() %>%
  mutate(
    p_adj = p.adjust(p_enrich, method = "BH"),
    p_lab = ifelse(p_adj < 1e-3, formatC(p_adj, format = "e", digits = 1),
                   signif(p_adj, 2))
  )

## Plot absoolute counts
df_left <- df_panelA %>%
  arrange(n_up_obs) %>%
  mutate(chr = factor(chr, levels = chr))

df_long_left <- df_left %>%
  select(chr, n_up_obs, n_up_expected) %>%
  pivot_longer(
    cols = c(n_up_obs, n_up_expected),
    names_to = "type",
    values_to = "count"
  ) %>%
  mutate(type = recode(type,
                       n_up_obs = "Observed",
                       n_up_expected = "Expected"))

p_left <- ggplot(df_long_left, aes(x = chr, y = count, fill = type)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8) +
  coord_flip() +
  scale_fill_manual(values = c(Expected = "grey70", Observed = "#F04E4E")) +
  theme_classic() +
  labs(x = NULL, y = "Number of upregulated genes", fill = NULL) +
  geom_text(
    data = df_left %>% mutate(y_lab = pmax(n_up_obs, n_up_expected) * 1.03),
    aes(x = chr, y = y_lab, label = paste0("FDR=", p_lab)),
    inherit.aes = FALSE,
    size = 3,
    hjust = 0
  ) +
  expand_limits(y = max(pmax(df_left$n_up_obs, df_left$n_up_expected)) * 1.15)

## Plot difference
df_right <- df_panelA %>%
  arrange(excess) %>%
  mutate(chr = factor(chr, levels = chr),
         is_chr19 = as.character(chr) == "chr19")

p_right <- ggplot(df_right, aes(x = chr, y = excess)) +
  geom_col(aes(fill = is_chr19), width = 0.8) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c(`TRUE` = "#7B4EA3",
                               `FALSE` = "grey75"),
                    guide = "none") +
  theme_classic() +
  labs(x = NULL, y = "Excess up genes (Observed − Expected)") +
  geom_text(
    aes(
      y = excess + ifelse(excess >= 0, 1, -1),
      label = paste0("FDR=", p_lab),
      hjust = ifelse(excess >= 0, -0.05, 1.05)
    ),
    size = 3
  ) +
  expand_limits(y = max(df_right$excess) * 1.15)

p_left
p_right


### GRANTA
## Prepare per-chromosome table 
genes_df <- as.data.frame(all.genes.ens.GRANTA_all) %>%
  transmute(
    chr = as.character(seqnames),
    Gene_Name = Gene_Name,
    padj,
    log2FoldChange
  ) %>%
  filter(grepl("^chr([0-9]+|X|Y)$", chr)) %>%
  distinct(chr, Gene_Name, .keep_all = TRUE) 

df_panelA <- genes_df %>%
  group_by(chr) %>%
  summarise(
    n_genes_chr = n(),
    n_up_obs    = sum(padj < 0.05 & log2FoldChange > 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    total_genes = sum(n_genes_chr),
    total_up    = sum(n_up_obs),
    p0          = total_up / total_genes,
    n_up_expected = total_up * (n_genes_chr / total_genes),
    enrichment    = n_up_obs / n_up_expected,
    excess        = n_up_obs - n_up_expected
  ) %>%
  rowwise() %>%
  mutate(
    p_enrich = binom.test(
      x = as.integer(n_up_obs),
      n = as.integer(n_genes_chr),
      p = p0,
      alternative = "greater"
    )$p.value
  ) %>%
  ungroup() %>%
  mutate(
    p_adj = p.adjust(p_enrich, method = "BH"),
    p_lab = ifelse(p_adj < 1e-3, formatC(p_adj, format = "e", digits = 1),
                   signif(p_adj, 2))
  )

## Plot absoolute counts
df_left <- df_panelA %>%
  arrange(n_up_obs) %>%
  mutate(chr = factor(chr, levels = chr))

df_long_left <- df_left %>%
  select(chr, n_up_obs, n_up_expected) %>%
  pivot_longer(
    cols = c(n_up_obs, n_up_expected),
    names_to = "type",
    values_to = "count"
  ) %>%
  mutate(type = recode(type,
                       n_up_obs = "Observed",
                       n_up_expected = "Expected"))

p_left_GRANTA <- ggplot(df_long_left, aes(x = chr, y = count, fill = type)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8) +
  coord_flip() +
  scale_fill_manual(values = c(Expected = "grey70", Observed = "#F04E4E")) +
  theme_classic() +
  labs(x = NULL, y = "Number of upregulated genes", fill = NULL) +
  geom_text(
    data = df_left %>% mutate(y_lab = pmax(n_up_obs, n_up_expected) * 1.03),
    aes(x = chr, y = y_lab, label = paste0("FDR=", p_lab)),
    inherit.aes = FALSE,
    size = 3,
    hjust = 0
  ) +
  expand_limits(y = max(pmax(df_left$n_up_obs, df_left$n_up_expected)) * 1.15)

## Plot difference
df_right <- df_panelA %>%
  arrange(excess) %>%
  mutate(chr = factor(chr, levels = chr),
         is_chr19 = as.character(chr) == "chr19")

p_right_GRANTA <- ggplot(df_right, aes(x = chr, y = excess)) +
  geom_col(aes(fill = is_chr19), width = 0.8) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c(`TRUE` = "#7B4EA3",
                               `FALSE` = "grey75"),
                    guide = "none") +
  theme_classic() +
  labs(x = NULL, y = "Excess up genes (Observed − Expected)") +
  geom_text(
    aes(
      y = excess + ifelse(excess >= 0, 1, -1),
      label = paste0("FDR=", p_lab),
      hjust = ifelse(excess >= 0, -0.05, 1.05)
    ),
    size = 3
  ) +
  expand_limits(y = max(df_right$excess) * 1.15)

p_left_GRANTA
p_right_GRANTA

################ Positional enrichment GSEA ####################################
### Helper functions
add_symbol_to_gr <- function(gr, ensg_from = c("names", "mcols"), ensg_col = NULL) {
  ensg_from <- match.arg(ensg_from)
  
  if (ensg_from == "names") {
    ensg <- names(gr)
  } else {
    if (is.null(ensg_col)) stop("Provide ensg_col when ensg_from='mcols'")
    ensg <- as.character(mcols(gr)[[ensg_col]])
  }
  
  if (is.null(ensg) || all(is.na(ensg)) || all(ensg == "")) {
    stop("ENSG IDs not found in GRanges names() or specified mcols column.")
  }
  
  ensg_clean <- sub("\\..*$", "", ensg)
  
  sym <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = ensg_clean,
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )
  
  mcols(gr)$SYMBOL <- unname(sym)
  gr
}

make_chr_pathways_from_gr_symbols <- function(gr,
                                              symbol_col = "SYMBOL",
                                              keep_chr = c(as.character(1:22), "X", "Y"),
                                              prefix = "CHR") {
  
  symbols <- as.character(mcols(gr)[[symbol_col]])
  chrs <- as.character(seqnames(gr))
  chrs <- gsub("^chr", "", chrs, ignore.case = TRUE)
  
  ok <- !is.na(symbols) & symbols != "" & chrs %in% keep_chr
  symbols <- symbols[ok]
  chrs <- chrs[ok]
  
  pathways <- split(unique(symbols), paste0(prefix, chrs))
  
  ## order CHR1..CHR22,X,Y
  ord <- order(match(names(pathways), paste0(prefix, keep_chr)))
  pathways[ord]
}


make_ranks_from_deseq <- function(res, stat_col = "stat") {
  df <- as.data.frame(res)
  genes <- rownames(df)
  
  ok <- !is.na(genes) & genes != "" & !is.na(df[[stat_col]]) & is.finite(df[[stat_col]])
  df <- df[ok, , drop = FALSE]
  genes <- genes[ok]
  
  ## de-duplicate symbols if needed
  df$gene <- genes
  df <- df[order(abs(df[[stat_col]]), decreasing = TRUE), ]
  df <- df[!duplicated(df$gene), ]
  
  ranks <- df[[stat_col]]
  names(ranks) <- df$gene
  sort(ranks, decreasing = TRUE)
}

run_fgsea_minimal <- function(pathways, ranks, out_csv = NULL,
                              minSize = 10, maxSize = 50000,
                              collapse_leadingEdge = TRUE) {
  
  fg <- fgsea(pathways = pathways, stats = ranks,
              minSize = minSize, maxSize = maxSize)
  
  fg <- as.data.frame(fg)
  fg <- fg[order(fg$padj, -abs(fg$NES)), ]
  
  if (!is.null(out_csv)) {
    fg_out <- fg
    if ("leadingEdge" %in% colnames(fg_out)) {
      if (collapse_leadingEdge) {
        fg_out$leadingEdge <- vapply(fg_out$leadingEdge,
                                     function(x) paste(x, collapse = ";"),
                                     FUN.VALUE = character(1))
      } else {
        fg_out$leadingEdge <- NULL
      }
    }
    write.csv(fg_out, out_csv, row.names = FALSE)
  }
  
  fg
}

plot_chr_bar <- function(fg_df,
                         main = "Chromosome enrichment (GSEA)",
                         ylab = "NES") {
  
  df <- fg_df
  df <- df[!is.na(df$NES) & !is.na(df$padj), ]
  df$pathway_clean <- tolower(df$pathway)
  df$log10padj <- -log10(df$padj)
  
  
  df$pathway_clean <- factor(
    df$pathway_clean,
    levels = df$pathway_clean[order(df$NES)]
  )
  
  ggplot(df, aes(x = pathway_clean, y = NES, fill = log10padj)) +
    geom_col() +
    coord_flip() +
    scale_fill_gradient(
      low = "grey90",
      high = "red3",
      name = "-log10(padj)"
    ) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = main, x = "Chromosome", y = ylab) +
    theme_classic() +
    theme(
      axis.text.y = element_text(color = "black"),
      axis.title.y = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5)
    )
}

### Rank genes and add annotation
## MCL all
ranks_MCL <- make_ranks_from_deseq(resMCL_naive_blood, stat_col = "stat")
gr_sym_MCL <- add_symbol_to_gr(all.genes.ens.MCL_all, ensg_from = "names")

## GRANTA
ranks_GRANTA <- make_ranks_from_deseq(resGRA_naive_blood, stat_col = "stat")
gr_sym_GRANTA <- add_symbol_to_gr(all.genes.ens.GRANTA_all, ensg_from = "names")

### Prepare chromosome sets 
chr_sets_sym_MCL <- make_chr_pathways_from_gr_symbols(gr_sym_MCL, symbol_col = "SYMBOL")
chr_sets_sym_MCL <- lapply(chr_sets_sym_MCL, function(gs) intersect(gs, names(ranks_MCL)))

chr_sets_sym_GRANTA <- make_chr_pathways_from_gr_symbols(gr_sym_GRANTA, symbol_col = "SYMBOL")
chr_sets_sym_GRANTA <- lapply(chr_sets_sym_GRANTA, function(gs) intersect(gs, names(ranks_GRANTA)))


### Run GSEA
fg_chr_genomewide_MCL <- run_fgsea_minimal(
  pathways = chr_sets_sym_MCL,
  ranks    = ranks_MCL,
  out_csv  = "results/chrGSEA_genomewide_MCL_vs_naive_blood.csv",
  minSize  = 10,
  maxSize  = 50000
)

fg_chr_genomewide_GRANTA <- run_fgsea_minimal(
  pathways = chr_sets_sym_GRANTA,
  ranks    = ranks_GRANTA,
  out_csv  = "results/chrGSEA_genomewide_GRANTA_vs_naive_blood.csv",
  minSize  = 10,
  maxSize  = 50000
)

plot_chr_bar(
  fg_chr_genomewide_MCL,
  main = "Chromosome enrichment (GSEA) - MCL vs naive_blood"
)

plot_chr_bar(
  fg_chr_genomewide_GRANTA,
  main = "Chromosome enrichment (GSEA) - GRANTA vs naive_blood"
)

plotEnrichment(chr_sets_sym[["CHR19"]], ranks_MCL)
plotEnrichment(chr_sets_sym[["CHR19"]], ranks_GRANTA)

########################## Tables with chr19 up genes ##########################
save_chr19_up_table_from_gr <- function(gr,
                                        out_csv,
                                        lfc_min = 1,
                                        padj_max = 0.05,
                                        chr_focus = "chr19",
                                        gene_col = "Gene_Name") {
  df <- as.data.frame(gr) %>%
    mutate(chr = as.character(seqnames)) %>%
    filter(chr == chr_focus) %>%
    filter(!is.na(.data[[gene_col]])) %>%
    filter(!is.na(log2FoldChange), !is.na(padj)) %>%
    filter(log2FoldChange > lfc_min, padj < padj_max) %>%
    transmute(
      gene = as.character(.data[[gene_col]]),
      gene_id = rownames(.),
      log2FoldChange,
      padj,
      chr,
      start = start,
      end   = end,
      strand = as.character(strand)
    ) %>%
    # one gene = one row, matching Venn "unique()"
    group_by(gene) %>%
    arrange(padj, desc(log2FoldChange), .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    arrange(padj, desc(log2FoldChange))
  
  write.csv2(df, out_csv, row.names = FALSE)
  df
}

# MCL table 
chr19_up_tbl_MCL <- save_chr19_up_table_from_gr(
  gr       = all.genes.ens.MCL_all,
  out_csv  = "results/chr19_up_genes_MCL_vs_naive_blood_coordinates.csv",
  lfc_min  = 1,
  padj_max = 0.05,
  chr_focus = "chr19",
  gene_col = "Gene_Name"
)

# GRANTA table 
chr19_up_tbl_GRANTA <- save_chr19_up_table_from_gr(
  gr       = all.genes.ens.GRANTA_all,
  out_csv  = "results/chr19_up_genes_GRANTA_vs_naive_blood_coordinates.csv",
  lfc_min  = 1,
  padj_max = 0.05,
  chr_focus = "chr19",
  gene_col = "Gene_Name"
)

get_chr19_up_set <- function(gr, lfc_min = 1, padj_max = 0.05, gene_col = "Gene_Name") {
  df <- as.data.frame(gr) %>%
    mutate(chr = as.character(seqnames)) %>%
    filter(chr == "chr19") %>%
    filter(!is.na(.data[[gene_col]])) %>%
    filter(!is.na(log2FoldChange), !is.na(padj)) %>%
    filter(log2FoldChange > lfc_min, padj < padj_max)
  
  unique(df[[gene_col]])
}

set_MCL    <- get_chr19_up_set(all.genes.ens.MCL_all,    lfc_min = 1, padj_max = 0.05, gene_col = "Gene_Name")
set_GRANTA <- get_chr19_up_set(all.genes.ens.GRANTA_all, lfc_min = 1, padj_max = 0.05, gene_col = "Gene_Name")

p_venn_up_chr19 <- ggvenn(
  list(MCL = set_MCL, GRANTA = set_GRANTA),
  fill_color = c("#993333", "#4DBBD5FF"),
  stroke_size = 0.8,
  set_name_size = 5,
  text_size = 5
)


################### Expression along chr19 #####################################
cts <- as.matrix(counts)
mcl_cols    <- grep("^MCL_", colnames(cts), value = TRUE)
ctrl_cols   <- grep("^naive_blood", colnames(cts), value = TRUE)
granta_cols <- grep("^GRANTA", colnames(cts), value = TRUE)

### size-factor normalization using Control + MCL + GRANTA
use_cols <- c(mcl_cols, ctrl_cols, granta_cols)

coldata <- data.frame(
  sample = use_cols,
  group  = dplyr::case_when(
    use_cols %in% mcl_cols    ~ "MCL",
    use_cols %in% ctrl_cols   ~ "Control",
    use_cols %in% granta_cols ~ "GRANTA"
  ),
  row.names = use_cols
)

dds <- DESeqDataSetFromMatrix(
  countData = cts[, use_cols, drop = FALSE],
  colData   = coldata,
  design    = ~ 1
)

dds <- estimateSizeFactors(dds)
norm <- counts(dds, normalized = TRUE)  

mcl_mean_expr    <- rowMeans(norm[, mcl_cols,    drop = FALSE], na.rm = TRUE)
ctrl_mean_expr   <- rowMeans(norm[, ctrl_cols,   drop = FALSE], na.rm = TRUE)
granta_mean_expr <- rowMeans(norm[, granta_cols, drop = FALSE], na.rm = TRUE)

### keep only chr19 genes
gr <- all.genes.ens.MCL_all  
gene_id_col <- if ("gene_id" %in% names(mcols(gr))) "gene_id" else
  if ("Gene_Name" %in% names(mcols(gr))) "Gene_Name" 

gene_id <- mcols(gr)[[gene_id_col]]

keep <- as.character(seqnames(gr)) == "chr19" & gene_id %in% rownames(cts)
gr19 <- gr[keep]
gene_id19 <- gene_id[keep]

df19 <- data.frame(
  gene_id      = gene_id19,
  pos          = start(gr19) + width(gr19) %/% 2,
  expr_MCL     = mcl_mean_expr[gene_id19],
  expr_Control = ctrl_mean_expr[gene_id19],
  expr_GRANTA  = granta_mean_expr[gene_id19],
  log2FC       = mcols(gr19)$log2FoldChange,
  stringsAsFactors = FALSE
) %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  arrange(pos)

### smoothing
k <- 101  

df19 <- df19 %>%
  mutate(
    expr_MCL_s     = zoo::rollmean(expr_MCL,     k = k, fill = NA, align = "center"),
    expr_Control_s = zoo::rollmean(expr_Control, k = k, fill = NA, align = "center"),
    expr_GRANTA_s  = zoo::rollmean(expr_GRANTA,  k = k, fill = NA, align = "center"),
    log2FC_s       = zoo::rollmean(log2FC,       k = k, fill = NA, align = "center")
  )

### probe coordinates
probe_start <- 478637
probe_end   <- 702132

probe_chr19 <- GRanges(
  seqnames = "chr19",
  ranges = IRanges(start = probe_start, end = probe_end)
)
probe_width <- width(probe_chr19)

### Plots
y_min <- min(df19$expr_Control_s, df19$expr_MCL_s, df19$expr_GRANTA_s, na.rm = TRUE)
y_max <- max(df19$expr_Control_s, df19$expr_MCL_s, df19$expr_GRANTA_s, na.rm = TRUE)
y_height <- (y_max - y_min) * 0.05

p_along19_patients <- ggplot(df19, aes(x = pos)) +
  geom_rect(
    xmin = probe_start,
    xmax = probe_end,
    ymin = y_min,
    ymax = y_min + y_height,
    inherit.aes = FALSE,
    fill = "#A275B3",
    alpha = 0.8
  ) +
  geom_line(aes(y = expr_Control_s), linewidth = 1, color = "#7CD3F7") +
  geom_line(aes(y = expr_MCL_s),     linewidth = 1, color = "#F26767") +
  scale_y_continuous(name = "Mean normalized expression") +
  scale_x_continuous(
    breaks = seq(0, 6e7, by = 5e6),
    labels = scales::label_number(scale = 1e-6, suffix = " Mb")
  ) +
  theme_classic() +
  labs(x = "chr19 coordinate", linetype = NULL) + 
  ggtitle("MCL patients")



p_along19_GRANTA <- ggplot(df19, aes(x = pos)) +
  geom_rect(
    xmin = probe_start,
    xmax = probe_end,
    ymin = y_min,
    ymax = y_min + y_height,
    inherit.aes = FALSE,
    fill = "#A275B3",
    alpha = 0.8
  ) +
  geom_line(aes(y = expr_Control_s), linewidth = 1, color = "#7CD3F7") +
  geom_line(aes(y = expr_GRANTA_s),  linewidth = 1, color = "#F26767") +
  scale_y_continuous(name = "Mean normalized expression") +
  scale_x_continuous(
    breaks = seq(0, 6e7, by = 5e6),
    labels = scales::label_number(scale = 1e-6, suffix = " Mb")
  ) +
  theme_classic() +
  labs(x = "chr19 coordinate", linetype = NULL) + 
  ggtitle("GRANTA")

################### Expression vs distance to FISH probe #######################
### Save tables for whole chromosome chr19
## MLC
saveRDS(all.genes.ens.MCL_all, "results/all.genes.ens.MCL_all.RDS")

chr19_genes_MCL <- all.genes.ens.MCL_all[seqnames(all.genes.ens.MCL_all) == "chr19"]
dtn <- distanceToNearest(chr19_genes_MCL, probe_chr19)
chr19_genes_MCL$dist_probe <- mcols(dtn)$distance

df_chr19_MCL <- as.data.frame(chr19_genes_MCL)
df_chr19_MCL$group <- ifelse(df_chr19_MCL$dist_probe <= 1e6, "Near (±1Mb)", "Far")
write.csv2(df_chr19_MCL, "results/chr19_degs_dist_to_probe_patients.csv")

## GRANTA
chr19_genes_GRANTA <- all.genes.ens.GRANTA_all[seqnames(all.genes.ens.GRANTA_all) == "chr19"]
dtn <- distanceToNearest(chr19_genes_GRANTA, probe_chr19)
chr19_genes_GRANTA$dist_probe <- mcols(dtn)$distance

df_chr19_GRANTA <- as.data.frame(chr19_genes_GRANTA)
df_chr19_GRANTA$group <- ifelse(df_chr19_GRANTA$dist_probe <= 1e6, "Near (±1Mb)", "Far")
write.csv2(df_chr19_GRANTA, "results/chr19_degs_dist_to_probe_GRANTA.csv")


### Calculate stats for the genes on the same chr arm (chr19p)
#"https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz"
cyto <- read.table(
  "cytoBand.txt",
  sep = "\t",
  header = FALSE,
  stringsAsFactors = FALSE
)
colnames(cyto) <- c("chr", "start", "end", "band", "gieStain")
cyto19 <- subset(cyto, chr == "chr19")
p_bands <- cyto19[grepl("^p", cyto19$band), ]
p_end <- max(p_bands$end)

### MCL patients
chr19_genes_MCL <- all.genes.ens.MCL_all[seqnames(all.genes.ens.MCL_all) == "chr19"]

# Restrict genes to chr19p using midpoint
mid <- start(chr19_genes_MCL) + (width(chr19_genes_MCL) %/% 2)
chr19p_genes_MCL <- chr19_genes_MCL[mid <= p_end]

dtn <- distanceToNearest(chr19p_genes_MCL, probe_chr19)
chr19p_genes_MCL$dist_probe <- mcols(dtn)$distance
chr19p_genes_MCL$log_dist <- log10(chr19p_genes_MCL$dist_probe + 1)

df_chr19p_MCL <- as.data.frame(chr19p_genes_MCL)
df_chr19p_MCL$group <- ifelse(df_chr19p_MCL$dist_probe <= 1e6, "Near (±1Mb)", "Far")
write.csv2(df_chr19p_MCL, "results/chr19p_degs_dist_to_probe_patients.csv")

## Linear model
lm_chr19p_MCL <- lm(log2FoldChange ~ log_dist, data=df_chr19p_MCL)
summary(lm_chr19p_MCL)

## Spearman correlation
sp_chr19p_MCL <- cor.test(df_chr19p_MCL$log_dist, df_chr19p_MCL$log2FoldChange, method="spearman")
sp_chr19p_MCL

rho_MCL  <- -0.2119405 
pval_MCL <- 1.986e-07
n_MCL    <- nrow(df_chr19p_MCL)

label_txt_MCL <- paste0(
  "Spearman \u03C1 = ", sprintf("%.3f", rho_MCL),
  "\n", "p = ", formatC(pval_MCL, format = "e", digits = 2),
  "\n", "n = ", n_MCL
)

## Random permutation test
# Observed Spearman rho
obs_rho <- rho_MCL
chr19p_len <- p_end 

set.seed(1)
B <- 10000  
rho_null <- numeric(B)

for (i in seq_len(B)) {
  # random start, keeping probe fully inside chr19p
  rs <- sample.int(chr19p_len - probe_width + 1, 1)
  rand_probe <- GRanges("chr19", IRanges(rs, rs + probe_w - 1))
  
  dtn_i <- distanceToNearest(chr19p_genes_MCL, rand_probe)
  log_dist_i <- log10(mcols(dtn_i)$distance + 1)
  
  rho_null[i] <- suppressWarnings(cor(log_dist_i, df_chr19p_MCL$log2FoldChange, method="spearman"))
}

# Empirical two-sided p-value
p_emp <- (sum(abs(rho_null) >= abs(obs_rho), na.rm=TRUE) + 1) / (sum(!is.na(rho_null)) + 1)
p_emp

quantile(rho_null, c(0.025, 0.5, 0.975), na.rm=TRUE)

## Plots
# Loess
p_linear_patients <- ggplot(df_chr19p_MCL, aes(x = log_dist, y = log2FoldChange)) +
  geom_point(alpha = 0.15, size = 1, shape = 16) +
  geom_smooth(method = "loess", color = "#B393C4", linewidth = 1.2) +
  theme_classic() +
  labs(x = "Distance to chr19 probe (log10 bp)",
       y = "log2FC") +
  annotate("text",
           x = min(df_chr19$log_dist, na.rm = TRUE) + 0.2,
           y = max(df_chr19$log2FoldChange, na.rm = TRUE),
           hjust = 0, vjust = 1,
           label = label_txt, size = 4)

# Empirical p
null_df <- data.frame(rho_null = rho_null)

p_emp_p_patients <- ggplot(null_df, aes(x = rho_null)) +
  stat_ecdf(geom = "step", linewidth = 0.9) +
  geom_vline(xintercept = obs_rho, color = "red", linewidth = 1) +
  theme_classic() +
  labs(x = "Spearman ρ (random probe positions)",
       y = "Empirical CDF") +
  annotate("text",
           x = min(rho_null, na.rm = TRUE),
           y = 1, hjust = 0, vjust = 1.2,
           label = paste0("Observed ρ = ", sprintf("%.3f", obs_rho),
                          "\nEmpirical p = ", sprintf("%.3f", p_emp),
                          "\nB = ", length(rho_null)),
           size = 4) +
  coord_cartesian(clip = "off")


# Binary ±1 Mb
summary_stats <- df_chr19p_MCL %>%
  group_by(group) %>%
  summarise(
    mean   = mean(log2FoldChange, na.rm = TRUE),
    median = median(log2FoldChange, na.rm = TRUE)
  )

p_binary_patients <- ggplot(df_chr19p_MCL, aes(x = group, y = log2FoldChange)) +
  geom_boxplot() +
  theme_classic() +
  stat_compare_means(method = "wilcox.test",
                     label = "p.format") +
  geom_text(data = summary_stats,
            aes(x = group,
                y = max(df_chr19p_MCL$log2FoldChange, na.rm = TRUE) * 0.9,
                label = paste0("Mean = ", round(mean, 2),
                               "\nMedian = ", round(median, 2))),
            inherit.aes = FALSE,
            size = 4) +
  ggtitle("MCL patients")

### GRANTA
chr19_genes_GRANTA <- all.genes.ens.GRANTA_all[seqnames(all.genes.ens.GRANTA_all) == "chr19"]
# Restrict genes to chr19p using midpoint
mid <- start(chr19_genes_GRANTA) + (width(chr19_genes_GRANTA) %/% 2)
chr19p_genes_GRANTA <- chr19_genes_GRANTA[mid <= p_end]

# Distance to probe
dtn <- distanceToNearest(chr19p_genes_GRANTA, probe_chr19)
chr19p_genes_GRANTA$dist_probe <- mcols(dtn)$distance
chr19p_genes_GRANTA$log_dist <- log10(chr19p_genes_GRANTA$dist_probe + 1)

df_chr19p_GRANTA <- as.data.frame(chr19p_genes_GRANTA)
df_chr19p_GRANTA$group <- ifelse(df_chr19p_GRANTA$dist_probe <= 1e6, "Near (±1Mb)", "Far")
write.csv2(df_chr19p_GRANTA, "results/chr19p_degs_dist_to_probe_GRANTA.csv.csv")

## Linear model
lm_chr19p_GRANTA <- lm(log2FoldChange ~ log_dist, data=df_chr19p_GRANTA)
summary(lm_chr19p_GRANTA)

## Spearman correlation
sp_chr19p_GRANTA <- cor.test(df_chr19p_GRANTA$log_dist, df_chr19p_GRANTA$log2FoldChange, method="spearman")
sp_chr19p_GRANTA

rho_GRANTA  <- -0.2307418 
pval_GRANTA <- 6.348e-06
n_GRANTA   <- nrow(df_chr19p_GRANTA)

label_txt_GRANTA <- paste0(
  "Spearman \u03C1 = ", sprintf("%.3f", rho_GRANTA),
  "\n", "p = ", formatC(pval_GRANTA, format = "e", digits = 2),
  "\n", "n = ", n_GRANTA
)

## Random permutation test
# Observed Spearman rho
obs_rho_GRANTA <- rho_GRANTA
chr19p_len <- p_end 

set.seed(1)
B <- 10000  
rho_null_GRANTA <- numeric(B)

for (i in seq_len(B)) {
  rs <- sample.int(chr19p_len - probe_width + 1, 1)
  rand_probe <- GRanges("chr19", IRanges(rs, rs + probe_w - 1))
  
  dtn_i <- distanceToNearest(chr19p_genes_GRANTA, rand_probe)
  log_dist_i <- log10(mcols(dtn_i)$distance + 1)
  
  rho_null_GRANTA[i] <- suppressWarnings(cor(log_dist_i, df_chr19p_GRANTA$log2FoldChange, method="spearman"))
}

# Empirical two-sided p-value
p_emp_GRANTA <- (sum(abs(rho_null_GRANTA) >= abs(obs_rho_GRANTA), na.rm=TRUE) + 1) / (sum(!is.na(rho_null_GRANTA)) + 1)
p_emp_GRANTA

quantile(rho_null_GRANTA, c(0.025, 0.5, 0.975), na.rm=TRUE)

## Plots
#Loess
p_linear_GRANTA <- ggplot(df_chr19p_GRANTA, aes(x = log_dist, y = log2FoldChange)) +
  geom_point(alpha = 0.15, size = 1, shape = 16) +
  geom_smooth(method = "loess", color = "#B393C4", linewidth = 1.2) +
  theme_classic() +
  labs(x = "Distance to chr19 probe (log10 bp)",
       y = "log2FC") +
  annotate("text",
           x = min(df_chr19$log_dist, na.rm = TRUE) + 0.2,
           y = max(df_chr19$log2FoldChange, na.rm = TRUE),
           hjust = 0, vjust = 1,
           label = label_txt, size = 4)

# Empirical p
null_df_GRANTA <- data.frame(rho_null = rho_null_GRANTA)

p_emp_p_GRANTA <- ggplot(null_df_GRANTA, aes(x = rho_null)) +
  stat_ecdf(geom = "step", linewidth = 0.9) +
  geom_vline(xintercept = obs_rho, color = "red", linewidth = 1) +
  theme_classic() +
  labs(x = "Spearman ρ (random probe positions)",
       y = "Empirical CDF") +
  annotate("text",
           x = min(rho_null_GRANTA, na.rm = TRUE),
           y = 1, hjust = 0, vjust = 1.2,
           label = paste0("Observed ρ = ", sprintf("%.3f", obs_rho_GRANTA),
                          "\nEmpirical p = ", sprintf("%.3f", p_emp),
                          "\nB = ", length(rho_null_GRANTA)),
           size = 4) +
  coord_cartesian(clip = "off")

# Binary ±1 Mb
summary_stats <- df_chr19p_GRANTA %>%
  group_by(group) %>%
  summarise(
    mean   = mean(log2FoldChange, na.rm = TRUE),
    median = median(log2FoldChange, na.rm = TRUE)
  )

p_binary_GRANTA <- ggplot(df_chr19p_GRANTA, aes(x = group, y = log2FoldChange)) +
  geom_boxplot() +
  theme_classic() +
  stat_compare_means(method = "wilcox.test",
                     label = "p.format") +
  geom_text(data = summary_stats,
            aes(x = group,
                y = max(df_chr19p_GRANTA$log2FoldChange, na.rm = TRUE) * 0.9,
                label = paste0("Mean = ", round(mean, 2),
                               "\nMedian = ", round(median, 2))),
            inherit.aes = FALSE,
            size = 4) +
  ggtitle("GRANTA")


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
#
#attached base packages:
#  [1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] lubridate_1.9.3                          forcats_1.0.0                            stringr_1.5.1                            purrr_1.0.2                             
#[5] readr_2.1.5                              tidyr_1.3.1                              tibble_3.2.1                             tidyverse_2.0.0                         
#[9] clusterProfiler_4.8.2                    TxDb.Hsapiens.UCSC.hg38.knownGene_3.17.0 rtracklayer_1.60.1                       Cairo_1.6-2                             
#[13] nucleR_2.32.0                            profileplyr_1.16.0                       csaw_1.33.0                              EnsDb.Hsapiens.v86_2.99.0               
#[17] ensembldb_2.24.1                         AnnotationFilter_1.24.0                  GenomicFeatures_1.52.2                   latticeExtra_0.6-30                     
#[21] lattice_0.22-6                           gtable_0.3.5                             gridExtra_2.3                            ggimage_0.3.3                           
#[25] ggupset_0.3.0                            ggplotify_0.1.2                          org.Hs.eg.db_3.17.0                      AnnotationDbi_1.62.2                    
#[29] VennDiagram_1.7.3                        futile.logger_1.4.3                      ggvenn_0.1.10                            dplyr_1.1.4                             
#[33] ggrepel_0.9.5                            ggplot2_3.5.1                            karyoploteR_1.26.0                       regioneR_1.32.0                         
#[37] DESeq2_1.40.2                            SummarizedExperiment_1.30.2              Biobase_2.60.0                           MatrixGenerics_1.12.3                   
#[41] matrixStats_1.3.0                        GenomicRanges_1.52.1                     GenomeInfoDb_1.36.4                      IRanges_2.34.1                          
#[45] S4Vectors_0.38.2                         BiocGenerics_0.46.0                     

#loaded via a namespace (and not attached):
#  [1] fs_1.6.4                                  ProtGenerics_1.32.0                       bitops_1.0-7                              enrichplot_1.20.0                        
#[5] doParallel_1.0.17                         HDO.db_0.99.1                             httr_1.4.7                                RColorBrewer_1.1-3                       
#[9] tools_4.3.0                               backports_1.5.0                           DT_0.33                                   utf8_1.2.4                               
#[13] R6_2.5.1                                  lazyeval_0.2.2                            GetoptLong_1.0.5                          withr_3.0.0                              
#[17] prettyunits_1.2.0                         preprocessCore_1.62.1                     cli_3.6.2                                 formatR_1.14                             
#[21] scatterpie_0.2.2                          labeling_0.4.3                            Rsamtools_2.16.0                          yulab.utils_0.1.4                        
#[25] gson_0.1.0                                foreign_0.8-86                            R.utils_2.12.3                            DOSE_3.26.2                              
#[29] dichromat_2.0-1                           plotrix_3.8-4                             BSgenome_1.68.0                           limma_3.56.2                             
#[33] rstudioapi_0.16.0                         RSQLite_2.3.7                             generics_0.1.3                            gridGraphics_0.5-1                       
#[37] TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2   shape_1.4.6.1                             BiocIO_1.10.0                             hwriter_1.3.2.1                          
#[41] gtools_3.9.5                              GO.db_3.17.0                              Matrix_1.6-5                              interp_1.1-6                             
#[45] fansi_1.0.6                               abind_1.4-7                               R.methodsS3_1.8.2                         lifecycle_1.0.4                          
#[49] yaml_2.3.8                                edgeR_3.42.4                              gplots_3.1.3.1                            qvalue_2.32.0                            
#[53] BiocFileCache_2.8.0                       blob_1.2.4                                promises_1.3.0                            crayon_1.5.2                             
#[57] cowplot_1.1.3                             KEGGREST_1.40.1                           magick_2.8.3                              ComplexHeatmap_2.16.0                    
#[61] pillar_1.9.0                              knitr_1.45                                metapod_1.7.0                             soGGi_1.32.0                             
#[65] fgsea_1.26.0                              rjson_0.2.21                              boot_1.3-30                               codetools_0.2-20                         
#[69] fastmatch_1.1-4                           glue_1.7.0                                ShortRead_1.58.0                          downloader_0.4                           
#[73] ggfun_0.1.4                               data.table_1.15.4                         vctrs_0.6.5                               png_0.1-8                                
#[77] treeio_1.24.3                             org.Mm.eg.db_3.17.0                       chipseq_1.50.0                            cachem_1.1.0                             
#[81] xfun_0.44                                 mime_0.12                                 TxDb.Mmusculus.UCSC.mm9.knownGene_3.2.2   S4Arrays_1.0.6                           
#[85] tidygraph_1.3.1                           pheatmap_1.0.12                           iterators_1.0.14                          rGREAT_2.2.0                             
#[89] nlme_3.1-164                              ggtree_3.8.2                              bit64_4.0.5                               progress_1.2.3                           
#[93] filelock_1.0.3                            KernSmooth_2.23-24                        rpart_4.1.23                              colorspace_2.1-1                         
#[97] DBI_1.2.2                                 Hmisc_5.1-2                               nnet_7.3-19                               tidyselect_1.2.1                         
#[101] bit_4.0.5                                 compiler_4.3.0                            curl_5.2.1                                htmlTable_2.4.2                          
#[105] bezier_1.1.2                              xml2_1.3.6                                TxDb.Mmusculus.UCSC.mm10.knownGene_3.10.0 DelayedArray_0.26.7                      
#[109] shadowtext_0.1.3                          checkmate_2.3.1                           scales_1.3.0                              caTools_1.18.2                           
#[113] ChIPseeker_1.36.0                         rappdirs_0.3.3                            tiff_0.1-12                               digest_0.6.35                            
#[117] rmarkdown_2.27                            XVector_0.40.0                            htmltools_0.5.8.1                         pkgconfig_2.0.3                          
#[121] jpeg_0.1-10                               base64enc_0.1-3                           dbplyr_2.5.0                              fastmap_1.2.0                            
#[125] rlang_1.1.3                               GlobalOptions_0.1.2                       htmlwidgets_1.6.4                         shiny_1.8.1.1                            
#[129] EnrichedHeatmap_1.30.0                    farver_2.1.2                              jsonlite_1.8.8                            BiocParallel_1.34.2                      
#[133] R.oo_1.26.0                               GOSemSim_2.26.1                           VariantAnnotation_1.46.0                  RCurl_1.98-1.14                          
#[137] magrittr_2.0.3                            Formula_1.2-6                             GenomeInfoDbData_1.2.10                   patchwork_1.2.0                          
#[141] munsell_0.5.1                             Rcpp_1.0.12                               ape_5.8                                   bamsignals_1.32.0                        
#[145] viridis_0.6.5                             stringi_1.8.4                             ggraph_2.2.1                              zlibbioc_1.46.0                          
#[149] MASS_7.3-60                               plyr_1.8.9                                parallel_4.3.0                            deldir_2.0-4                             
#[153] Biostrings_2.68.1                         graphlayouts_1.1.1                        splines_4.3.0                             hms_1.1.3                                
#[157] circlize_0.4.16                           locfit_1.5-9.9                            igraph_2.0.3                              reshape2_1.4.4                           
#[161] biomaRt_2.56.1                            futile.options_1.0.1                      XML_3.99-0.16.1                           evaluate_0.23                            
#[165] biovizBase_1.48.0                         lambda.r_1.2.4                            tzdb_0.4.0                                httpuv_1.6.15                            
#[169] foreach_1.5.2                             tweenr_2.0.3                              polyclip_1.10-6                           clue_0.3-65                              
#[173] ggforce_0.4.2                             xtable_1.8-4                              restfulr_0.0.15                           tidytree_0.4.6                           
#[177] later_1.3.2                               viridisLite_0.4.2                         aplot_0.2.2                               memoise_2.0.1                            
#[181] GenomicAlignments_1.36.0                  cluster_2.1.6                             timechange_0.3.0      


