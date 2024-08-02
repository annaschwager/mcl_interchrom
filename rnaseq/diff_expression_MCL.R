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


