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
library(gridExtra)
library(gtable)
library(latticeExtra)
library(biomaRt)
library(FactoMineR)
library(factoextra)
library(clusterProfiler)
library(fgsea)
library(msigdbr)
library(enrichplot)
library(DOSE)
library(dplyr)
library(ggpubr)
library(GenomicRanges)

################### Loading the count table ##############################
counts = read.table("salmon.merged.gene_counts.tsv", header = T)
rownames(counts) <- counts$gene_id
counts <- counts[,-1]
counts <- counts[,-1]

######################### Quality checks ################################
par(mfrow=c(3,3))
#Histograms of log transformed counts to see reads distribution
hist(log(as.numeric(counts[,2])+1), main = "BLAS_1", xlab = "log(counts+1)")
hist(log(as.numeric(counts[,3])+1), main = "BLAS_1", xlab = "log(counts+1)")
hist(log(as.numeric(counts[,4])+1), main = "BLAS_1", xlab = "log(counts+1)")

hist(log(as.numeric(counts[,5])+1), main = "BLAS_abe_1", xlab = "log(counts+1)")
hist(log(as.numeric(counts[,6])+1), main = "BLAS_abe_2", xlab = "log(counts+1)")
hist(log(as.numeric(counts[,7])+1), main = "BLAS_abe_3", xlab = "log(counts+1)")

hist(log(as.numeric(counts[,8])+1), main = "BLAS_min_1", xlab = "log(counts+1)")
hist(log(as.numeric(counts[,9])+1), main = "BLAS_min_2", xlab = "log(counts+1)")
hist(log(as.numeric(counts[,10])+1), main = "BLAS_min_3", xlab = "log(counts+1)")

hist(log(as.numeric(counts[,8])+1), main = "GRANTA_1", xlab = "log(counts+1)")
hist(log(as.numeric(counts[,9])+1), main = "GRANTA_2", xlab = "log(counts+1)")
hist(log(as.numeric(counts[,10])+1), main = "GRANTA_3", xlab = "log(counts+1)")

hist(log(as.numeric(counts[,8])+1), main = "GRANTA_abe_1", xlab = "log(counts+1)")
hist(log(as.numeric(counts[,9])+1), main = "GRANTA_abe_2", xlab = "log(counts+1)")
hist(log(as.numeric(counts[,10])+1), main = "GRANTA_abe_3", xlab = "log(counts+1)")

hist(log(as.numeric(counts[,8])+1), main = "GRANTA_min_1", xlab = "log(counts+1)")
hist(log(as.numeric(counts[,9])+1), main = "GRANTA_min_2", xlab = "log(counts+1)")
hist(log(as.numeric(counts[,10])+1), main = "GRANTA_min_3", xlab = "log(counts+1)")
par(mfrow=c(1,1))

#comparing library sizes
boxplot(log(counts + 1))

#filtering low counts
colSums(counts)
counts_per_gene <- rowSums(counts)
table(counts_per_gene > 0)
counts <- counts[which(counts_per_gene > 0), ]
dim(counts)
boxplot(log(counts + 1))

################### Differential expression analysis with DESeq2 ##############
cond <- factor(gsub("_[0-9]", "", colnames(counts)))

for (i in 1:18) {
  counts[, i] <- as.numeric(gsub("\\..*", "", counts[, i]))
}

dds <- DESeqDataSetFromMatrix(counts, colData = DataFrame(cond), design = ~ cond)
colData(dds)
dds <- DESeq(dds)

### annotating and saving results
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") 

resGmin_G <- results(dds, contrast = c("cond", "GRANTA_min", "GRANTA"))
summary(resGmin_G)
resGmin_G  <- resGmin_G[complete.cases(resGmin_G ),]  
resGmin_G  <- resGmin_G[order(resGmin_G $padj),]
Gmin_G <- getBM(attributes = c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol', 'genecards', 'ensembl_gene_id', 'description', 'wikigene_description'),
                filters = 'hgnc_symbol', values = row.names(resGmin_G),
                mart = ensembl)
Gmin_G_df <- as.data.frame(resGmin_G)
Gmin_G_df$id <- row.names(resGmin_G)
Gmin_G_df <- left_join(Gmin_G_df, Gmin_G, by = join_by("id" == "hgnc_symbol"))
write.csv2(Gmin_G_df, "results/resGmin_G.csv")

resBmin_B <- results(dds, contrast = c("cond", "BLAS_min", "BLAS"))
summary(resBmin_B)
resBmin_B <- resBmin_B[complete.cases(resBmin_B),]  
resBmin_B  <- resBmin_B[order(resBmin_B$padj),]
Bmin_B <- getBM(attributes = c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol', 'genecards', 'ensembl_gene_id', 'description', 'wikigene_description'),
                filters = 'hgnc_symbol', values = row.names(resBmin_B),
                mart = ensembl)
Bmin_B_df <- as.data.frame(resBmin_B)
Bmin_B_df$id <- row.names(resBmin_B)
Bmin_B_df <- left_join(Bmin_B_df, Bmin_B, by = join_by("id" == "hgnc_symbol"))
write.csv2(Bmin_B_df, "results/resBmin_B.csv")

resG_B <- results(dds, contrast = c("cond", "GRANTA", "BLAS"))
summary(resG_B)
resG_B <- resG_B[complete.cases(resG_B),]  
resG_B  <- resG_B[order(resG_B$padj),]
G_B <- getBM(attributes = c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol', 'genecards', 'ensembl_gene_id', 'description', 'wikigene_description'),
             filters = 'hgnc_symbol', values = row.names(resG_B),
             mart = ensembl)
G_B_df <- as.data.frame(resG_B)
G_B_df$id <- row.names(resG_B)
G_B_df <- left_join(G_B_df, G_B, by = join_by("id" == "hgnc_symbol"))
write.csv2(G_B_df, "results/resG_B.csv")

resGabe_G <- results(dds, contrast = c("cond", "GRANTA_abe", "GRANTA"))
summary(resGabe_G)
resGabe_G <- resGabe_G[complete.cases(resGabe_G),]  
resGabe_G  <- resGabe_G[order(resGabe_G$padj),]
Gabe_G <- getBM(attributes = c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol', 'genecards', 'ensembl_gene_id', 'description', 'wikigene_description'),
                filters = 'hgnc_symbol', values = row.names(resGabe_G),
                mart = ensembl)
Gabe_G_df <- as.data.frame(resGabe_G)
Gabe_G_df$id <- row.names(resGabe_G)
Gabe_G_df <- left_join(Gabe_G_df, Gabe_G, by = join_by("id" == "hgnc_symbol"))
write.csv2(Gabe_G_df, "results/resGabe_G.csv")

resBabe_B <- results(dds, contrast = c("cond", "BLAS_abe", "BLAS"))
summary(resBabe_B)
resBabe_B <- resBabe_B[complete.cases(resBabe_B),]  
resBabe_B  <- resBabe_B[order(resBabe_B$padj),]
Babe_B <- getBM(attributes = c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol', 'genecards', 'ensembl_gene_id', 'description', 'wikigene_description'),
                filters = 'hgnc_symbol', values = row.names(resBabe_B),
                mart = ensembl)
Babe_B_df <- as.data.frame(resBabe_B)
Babe_B_df$id <- row.names(resBabe_B)
Babe_B_df <- left_join(Babe_B_df, Babe_B, by = join_by("id" == "hgnc_symbol"))
write.csv2(Babe_B_df, "results/resBabe_B.csv")

##filtering 
resG_B_sign <- resG_B[resG_B[,'padj'] < 0.05, ]
resG_B_up <- resG_B_sign[resG_B_sign[,'log2FoldChange'] > 1, ]
resG_B_down <- resG_B_sign[resG_B_sign[,'log2FoldChange'] < -1, ]

resGmin_G_sign <- resGmin_G[resGmin_G[,'padj'] < 0.05, ]
resGmin_G_up <- resGmin_G_sign[resGmin_G_sign[,'log2FoldChange'] > 1, ]
resGmin_G_down <- resGmin_G_sign[resGmin_G_sign[,'log2FoldChange'] < -1, ]
summary(resGmin_G_up)
summary(resGmin_G_down)

resGabe_G_sign <- resGabe_G[resGabe_G[,'padj'] < 0.05, ]
resGabe_G_up <- resGabe_G_sign[resGabe_G_sign[,'log2FoldChange'] > 1, ]
resGabe_G_down <- resGabe_G_sign[resGabe_G_sign[,'log2FoldChange'] < -1, ]
summary(resGabe_G_up)
summary(resGabe_G_down)

resBmin_B_sign <- resBmin_B[resBmin_B[,'padj'] < 0.05, ]
resBmin_B_up <- resBmin_B_sign[resBmin_B_sign[,'log2FoldChange'] > 1, ]
resBmin_B_down <- resBmin_B_sign[resBmin_B_sign[,'log2FoldChange'] < -1, ]

resBabe_B_sign <- resBabe_B[resBabe_B[,'padj'] < 0.05, ]
resBabe_B_up <- resBabe_B_sign[resBabe_B_sign[,'log2FoldChange'] > 1, ]
resBabe_B_down <- resBabe_B_sign[resBabe_B_sign[,'log2FoldChange'] < -1, ]

################################ Venn diagrams ###############################
granta_treatments <- list(A = row.names(resG_B_up),
                          B = row.names(resGmin_G_down),
                          C = row.names(resGabe_G_down))
names(granta_treatments) <- c("GRANTA vs BLAS up",
                              "GRANTA + Min vs GRANTA down",
                              "GRANTA + Abe vs GRANTA down")

p1 <- ggvenn(granta_treatments, 
                fill_color = c("indianred", 'lavender', 'lavender'),
                stroke_size = 0.5, set_name_size = 4,
                show_percentage = FALSE)

blas_granta_treatments_down <- list(A = row.names(resGmin_G_down),
                               B = row.names(resGabe_G_down),
                               C = row.names(resBmin_B_down),
                               D = row.names(resBabe_B_down))
names(blas_granta_treatments_down) <- c("GRANTA + Min vs GRANTA down",
                                        "GRANTA + Abe vs GRANTA down",
                                        "BLAS + Min vs BLAS down",
                                        "BLAS + Abe vs BLAS down")
p2 <- ggvenn(blas_granta_treatments_down, 
             fill_color = c('lavender', 'lavender', 'lavender', 'lavender'),
             stroke_size = 0.5, set_name_size = 4,
             show_percentage = FALSE)


blas_granta_treatments_up <- list(A = row.names(resGmin_G_up),
                                    B = row.names(resGabe_G_up),
                                    C = row.names(resBmin_B_up),
                                    D = row.names(resBabe_B_up))
names(blas_granta_treatments_up) <- c("GRANTA + Min vs GRANTA up",
                                        "GRANTA + Abe vs GRANTA up",
                                        "BLAS + Min vs BLAS up",
                                        "BLAS + Abe vs BLAS up")
p3 <- ggvenn(blas_granta_treatments_up, 
             fill_color = rep("indianred", 4),
             stroke_size = 0.5, set_name_size = 4,
             show_percentage = FALSE)

### Intersections with the patients results 

# loading the results from the MCL patients vs B naive comparison
MCL_vs_naive <- read.csv("./resMCL_naive_blood.csv")
rownames(MCL_vs_naive) <- MCL_vs_naive$X
mcl_up <- MCL_vs_naive[MCL_vs_naive$log2FoldChange > 1 & MCL_vs_naive$padj < 0.05,]

# loading enhancer predictions in MCL
abc <- read.delim("./EnhancerPredictions.txt")
abc$TargetGene <- gsub("\\.[0-9]", "", abc$TargetGene)
# retriving gene names
abc_annot <- getBM(attributes = c('hgnc_symbol', 'ensembl_transcript_id'),
                   filters = 'ensembl_transcript_id', values = abc$TargetGene,
                   mart = ensembl)
abc <- left_join(abc, abc_annot, by = join_by("TargetGene" == "ensembl_transcript_id"))
#as many enhancers can regulate many genes, will keep only unique gene rows
genes_regulated_by_ehnancers <- abc %>% distinct(hgnc_symbol, .keep_all = TRUE)

#MCL patients and GRANTA cells
mcl_granta <- list(A = mcl_up$Gene_Name,
                   B = row.names(resG_B_up),
                   C = genes_regulated_by_ehnancers$hgnc_symbol)

names(mcl_granta) <- c("MCL_patients_vs_naive_up",
                       "GRANTA_vs_BLAS_up",
                       "Genes regulated by enhancers in MCL patients")

venn_mcl_granta <- ggvenn(mcl_granta,
                         stroke_size = 0.5, set_name_size = 4,
                         show_percentage = FALSE)

# enhancers and Min
enhancers_min  <- list(A = row.names(resG_B_up),
                       B = genes_regulated_by_ehnancers$hgnc_symbol,
                       C = row.names(resGmin_G_down))

names(enhancers_min) <- c("GRANTA_vs_BLAS_up",
                          "Genes regulated by enhancers in MCL patients",
                          "GRANTA + Min vs GRANTA down")

venn_enhancers_min  <- ggvenn(enhancers_min,
                              stroke_size = 0.5, set_name_size = 4,
                              show_percentage = FALSE)

# enhancers and Abe
enhancers_abe  <- list(A = row.names(resG_B_up),
                       B = row.names(resGabe_G_down),
                       C = genes_regulated_by_ehnancers$hgnc_symbol)

names(enhancers_abe) <- c("GRANTA_vs_BLAS_up",
                          "GRANTA + Abe vs GRANTA down",
                          "Genes regulated by enhancers in MCL patients")

venn_enhancers_abe  <- ggvenn(enhancers_abe,
                              stroke_size = 0.5, set_name_size = 4,
                              show_percentage = FALSE)


######################## Hypergeometric test on intersections ##############
granta_min <- list(A = row.names(resG_B_up),
                   B = row.names(resGmin_G_down))

names(granta_min) <- c("GRANTA vs BLAS up",
                       "GRANTA + Min vs GRANTA down")

p1_min <- ggvenn(granta_min, 
                 fill_color = c("indianred", 'lavender'),
                 stroke_size = 0.5, set_name_size = 4,
                 show_percentage = FALSE)

N = length(rownames(counts))
A = 2382+479
B = 1084+479
x = 479

phyper(q = x - 1,
       m = A,
       n = N-A,
       k = B,
       lower.tail = FALSE)

########################## PCA ##########################
rld <- rlog(dds)
plotPCA(rld, intgroup="cond")
plotPCA(rld, intgroup="cond", returnData = T)

#more detailed PCA using FactoMineR
resPCA <- PCA(t(assay(rld)),graph=FALSE) #using the rlog transformed counts obtained before
barplot(resPCA$eig[,2],main="Variance explained by all PCs", ylab ="explained variance, %")
# most of the variance is explained by the first 4 components

#PC1 and PC2
g1_cond <- fviz_pca_ind(resPCA, habillage = colData(dds)$cond, geom="point", repel = "TRUE", title ="PC1 and PC2", pointsize = 3)
g1_cond <- g1_cond + theme(axis.text=element_text(size=12), axis.title=element_text(size=14),title=element_text(size=14))
#PC3 and PC4
g2_cond <- fviz_pca_ind(resPCA, axes=c(3,4), habillage = colData(dds)$cond, geom="point", repel = "TRUE", title ="PC3 and PC4", pointsize = 3)
g2_cond <- g2_cond + theme(axis.text=element_text(size=12), axis.title=element_text(size=14),title=element_text(size=14))

########################## Intersections  ##########################
up_in_granta_down_by_min <- intersect(row.names(resG_B_up), row.names(resGmin_G_down))
up_in_granta_down_by_min <-resG_B[up_in_granta_down_by_min, ]
write.csv2(up_in_granta_down_by_min, "results/up_in_granta_down_by_min.csv") 

up_in_granta_down_by_abe <- intersect(row.names(resG_B_up), row.names(resGabe_G_down))
up_in_granta_down_by_abe <-resG_B[up_in_granta_down_by_abe, ]
write.csv2(up_in_granta_down_by_abe, "results/up_in_granta_down_by_abe.csv") 

################################## Volcano plots  ##############################
volc_granta_c_min <- ggplot(data=as.data.frame(resG_B)) + geom_point(aes(x=log2FoldChange, y=-log10(padj)), col="grey") +
  geom_point(data=as.data.frame(resG_B_down), aes(x=log2FoldChange, y=-log10(padj)), col="skyblue4") +
  geom_point(data=as.data.frame(resG_B_up), aes(x=log2FoldChange, y=-log10(padj)), col="indianred") +
  geom_point(data=as.data.frame(up_in_granta_down_by_min), aes(x=log2FoldChange, y=-log10(padj)), col="rosybrown2") +
  theme_minimal() +
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5)) +
  ggtitle("GRANTA vs BLAS downs by Min")

volc_granta_c_abe <- ggplot(data=as.data.frame(resG_B)) + geom_point(aes(x=log2FoldChange, y=-log10(padj)), col="grey") +
  geom_point(data=as.data.frame(resG_B_down), aes(x=log2FoldChange, y=-log10(padj)), col="skyblue4") +
  geom_point(data=as.data.frame(resG_B_up), aes(x=log2FoldChange, y=-log10(padj)), col="indianred") +
  geom_point(data=as.data.frame(up_in_granta_down_by_abe), aes(x=log2FoldChange, y=-log10(padj)), col="rosybrown2") +
  theme_minimal() +
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5)) +
  ggtitle("GRANTA vs BLAS downs by Abe")

volc_granta_granta_min <- ggplot(data=as.data.frame(resGmin_G)) + geom_point(aes(x=log2FoldChange, y=-log10(padj)), col="grey") +
  geom_point(data=as.data.frame(resGmin_G_down), aes(x=log2FoldChange, y=-log10(padj)), col="skyblue4") +
  geom_point(data=as.data.frame(resGmin_G_up), aes(x=log2FoldChange, y=-log10(padj)), col="indianred") +
  theme_minimal() +
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5)) +
  ggtitle("GRANTA + Min vs GRANTA")

volc_granta_granta_abe <- ggplot(data=as.data.frame(resGabe_G)) + geom_point(aes(x=log2FoldChange, y=-log10(padj)), col="grey") +
  geom_point(data=as.data.frame(resGabe_G_down), aes(x=log2FoldChange, y=-log10(padj)), col="skyblue4") +
  geom_point(data=as.data.frame(resGabe_G_up), aes(x=log2FoldChange, y=-log10(padj)), col="indianred") +
  theme_minimal() +
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5)) +
  ggtitle("GRANTA + Abe vs GRANTA")
 
############################## Enrichment analysis ################################
### GO ###
listA <- list(up = rownames(resGabe_G_up), down = rownames(resGabe_G_down))
listM <- list(up = rownames(resGmin_G_up), down = rownames(resGmin_G_down))

compGO_A <- compareCluster(geneCluster   = listA,
                           fun           = "enrichGO",
                           pvalueCutoff  = 0.05,
                           keyType       = "SYMBOL",
                           pAdjustMethod = "BH",
                           OrgDb='org.Hs.eg.db',
                           ont           = "BP",
                           universe = rownames(counts))
write.csv2(compGO_A, "./results/go_enrichment_BP_granta_Abema.csv")
dotplot(compGO_A, showCategory = 10, title = "go_enrichment_BP_Abema")

compGO_M <- compareCluster(geneCluster   = listM,
                            fun           = "enrichGO",
                            pvalueCutoff  = 0.05,
                            keyType       = "SYMBOL",
                            pAdjustMethod = "BH",
                            OrgDb='org.Hs.eg.db',
                            ont           = "BP",
                            universe = rownames(counts))
write.csv2(compGO_M, "./results/go_enrichment_BP__granta_Min.csv")
dotplot(compGO_M, showCategory = 10, title = "go_enrichment_BP_Min")

compGO_Amf <- compareCluster(geneCluster   = listA,
                           fun           = "enrichGO",
                           pvalueCutoff  = 0.05,
                           keyType       = "SYMBOL",
                           pAdjustMethod = "BH",
                           OrgDb='org.Hs.eg.db',
                           ont           = "MF",
                           universe = rownames(counts))
write.csv2(compGO_Amf, "./results/go_enrichment_MF_granta_Abema.csv")
dotplot(compGO_Amf, showCategory = 10, title = "go_enrichment_MF_Abema")

compGO_Mmf <- compareCluster(geneCluster   = listM,
                           fun           = "enrichGO",
                           pvalueCutoff  = 0.05,
                           keyType       = "SYMBOL",
                           pAdjustMethod = "BH",
                           OrgDb='org.Hs.eg.db',
                           ont           = "MF",
                           universe = rownames(counts))
write.csv2(compGO_Mmf, "./results/go_enrichment_MF_granta_Min.csv")
dotplot(compGO_Mmf, showCategory = 10, title = "go_enrichment_MF_Min")

### GSEA ###
abe_gene_list <- resGabe_G_001$log2FoldChange
names(abe_gene_list) <- row.names(resGabe_G_001)
abe_gene_list = sort(abe_gene_list, decreasing = TRUE)

gse_abe_bp <- gseGO(geneList=abe_gene_list, 
                    ont ="BP", 
                    keyType = "SYMBOL", 
                    minGSSize = 3, 
                    maxGSSize = 800, 
                    pvalueCutoff = "0.05", 
                    verbose = TRUE, 
                    OrgDb = org.Hs.eg.db, 
                    pAdjustMethod = "BH")
write.csv2(gse_abe_bp, "results/gse_abe_granta_bp.csv")

gse_abe_bp <- pairwise_termsim(gse_abe_bp)
treeplot(gse_abe_bp, showCategory = 40, color = 'NES')
treeplot(gse_abe_bp, showCategory = 60, color = 'NES')
#dotplot(gse_abe_bp, showCategory = 20)


min_gene_list <- resGmin_G_001$log2FoldChange
names(min_gene_list) <- row.names(resGmin_G_001)
min_gene_list= sort(min_gene_list, decreasing = TRUE)
gse_min_bp <- gseGO(geneList=min_gene_list, 
                      ont ="BP", 
                      keyType = "SYMBOL", 
                      minGSSize = 3, 
                      maxGSSize = 800, 
                      pvalueCutoff = "none", 
                      verbose = TRUE, 
                      OrgDb = org.Hs.eg.db, 
                      pAdjustMethod = "BH")
write.csv2(gse_min_bp, "results/gse_min_granta_bp.csv")

gse_min_bp <- pairwise_termsim(gse_min_bp)
treeplot(gse_min_bp, showCategory = 40, color = 'NES')
treeplot(gse_min_bp, showCategory = 60, color = 'NES')


### Selected plots
df_go <- as.data.frame(compGO_M)

df_go_down <- df_go %>%
  filter(grepl("down", Cluster, ignore.case = TRUE)) %>%   # adjust if needed
  mutate(
    GeneRatio_num = sapply(GeneRatio, function(x) eval(parse(text = x)))
  ) %>%
  arrange(p.adjust) %>%
  slice_head(n =10) %>%
  mutate(
    Description = factor(Description, levels = rev(Description))
  )

ggplot(df_go_down,
       aes(x = reorder(Description, GeneRatio_num),
           y = GeneRatio_num,
           fill = p.adjust)) +
  geom_col(width = 0.8) +
  coord_flip() +
  theme_classic() +
  scale_fill_gradient(
    low = "#3E667C",
    high = "#FFCCCC",
    trans = "reverse",
    name = "adj. p-value"
  ) +
  labs(
    title = "GO enrichment of down genes",
    x = "GO Biological Process",
    y = "Gene ratio"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 10)
  )

############## Barplots DEGs per chromosome ########################
chromsizes <- read.csv2("./chromSizes.csv", header = TRUE)
chromsizes$normsizes <- chromsizes$size/min(chromsizes$size)
chromsizes$normnumber <- chromsizes$gene_number/min(chromsizes$gene_number)
chromsizes$chr <- gsub("chr", "", chromsizes$chr)

## down DEGs Abe 
df_a_down <- as.data.frame(table(Gabe_G_df[Gabe_G_df$log2FoldChange < -1 & Gabe_G_df$padj < 0.05,]$chromosome_name))
df_a_down <- df_a_down[1:22,]
colnames(df_a_down) <- c("chr", "freq")
df_a_down <- merge(df_a_down,chromsizes, by.x = "chr", by.y = "chr")
df_a_down$normfreq_size <- df_a_down$freq/df_a_down$normsizes #norm by size
df_a_down$normfreq_number <- df_a_down$freq/df_a_down$normnumber #norm by gene number
write.csv(df_a_down, "results/down_degs_per_chromosome_abe.csv")

df_a_down <- df_a_down[order(df_a_down$freq,decreasing=TRUE),]
barplot(df_a_down$freq, names.arg = df_a_down$chr, ylab = "N down DEGs Abe") 
df_a_down <- df_a_down[order(df_a_down$normfreq_size,decreasing=TRUE),]
barplot(df_a_down$normfreq_size, names.arg = df_a_down$chr, ylab = "N down DEGs Abe normalised to chromSize patients")
df_a_down <- df_a_down[order(df_a_down$normfreq_number,decreasing=TRUE),]
barplot(df_a_down$normfreq_number, names.arg = df_a_down$chr, ylab = "N down DEGs Abe normalised to gene number patients")

## down DEGs Min
df_m_down <- as.data.frame(table(Gmin_G_df[Gmin_G_df$log2FoldChange < -1 & Gmin_G_df$padj < 0.05,]$chromosome_name))
df_m_down <- df_m_down[1:22,]
colnames(df_m_down) <- c("chr", "freq")
df_m_down <- merge(df_m_down,chromsizes, by.x = "chr", by.y = "chr")
df_m_down$normfreq_size <- df_m_down$freq/df_m_down$normsizes #norm by size
df_m_down$normfreq_number <- df_m_down$freq/df_m_down$normnumber #norm by gene number
write.csv(df_m_down, "results/down_degs_per_chromosome_min.csv")

df_m_down <- df_m_down[order(df_m_down$freq,decreasing=TRUE),]
barplot(df_m_down$freq, names.arg = df_m_down$chr, ylab = "N down DEGs Min50") 
df_m_down <- df_m_down[order(df_m_down$normfreq_size,decreasing=TRUE),]
barplot(df_m_down$normfreq_size, names.arg = df_m_down$chr, ylab = "N down DEGs Min50 normalised to chromSize patients")
df_m_down <- df_m_down[order(df_m_down$normfreq_number,decreasing=TRUE),]
barplot(df_m_down$normfreq_number, names.arg = df_m_down$chr, ylab = "N down DEGs Min50 normalised to gene number patients")

########################### Chr19 vs others ###################################
df_plot <- Gmin_G_df %>%
  mutate(
    chromosome_name = as.character(chromosome_name),
    group = ifelse(chromosome_name == "19", "chr19", "other")
  ) %>%
  filter(!grepl("PATCH|CHR_", chromosome_name)) %>%
  filter(!is.na(group))

#Remove points above 99.9th percentile
threshold <- quantile(df_plot$log2FoldChange, 0.999)
df_plot_trim <- df_plot %>%
  filter(log2FoldChange < threshold)

p_binary_granta <- ggplot(df_plot_trim, aes(x = group, y = log2FoldChange)) +
  geom_boxplot() +
  theme_classic() +
  stat_compare_means(method = "wilcox.test",
                     label = "p.format")


############ Transcriptional reversal scatters #######################
### Classify all genes as chr19 near probe, chr19 others, others

probe_start <- 478637
probe_end   <- 702132

probe_chr19 <- GRanges(
  seqnames = "chr19",
  ranges = IRanges(start = probe_start, end = probe_end)
)
probe_width <- width(probe_chr19)

near_bp <- 1e6

df_lfc <- Gmin_G_df %>%
  transmute(
    gene = id,
    chr  = as.character(chromosome_name),
    start = as.numeric(start_position),
    end   = as.numeric(end_position),
    lfc   = as.numeric(log2FoldChange),
    padj  = as.numeric(padj)
  ) %>%
  # keep canonical chroms only; normalize to "chr*"
  mutate(chr = paste0("chr", sub("^chr", "", chr))) %>%
  filter(grepl("^chr([0-9]+|X|Y)$", chr)) %>%
  filter(!is.na(start), !is.na(end))

gr_all <- GRanges(
  seqnames = df_lfc$chr,
  ranges   = IRanges(start = df_lfc$start, end = df_lfc$end),
  gene     = df_lfc$gene,
  lfc      = df_lfc$lfc,
  padj     = df_lfc$padj
)

# de-duplicate repeated genes like ACTG1:
# keep the entry with smallest padj, then largest |lfc|
o <- order(mcols(gr_all)$padj, -abs(mcols(gr_all)$lfc))
gr_all <- gr_all[o]
gr_all <- gr_all[!duplicated(mcols(gr_all)$gene)]

# Compute distances 
gr19 <- gr_all[seqnames(gr_all) == "chr19"]
dtn  <- distanceToNearest(gr19, probe_chr19)
mcols(gr19)$dist_probe <- mcols(dtn)$distance

df_out <- bind_rows(
  as.data.frame(gr19) %>%
    transmute(
      gene = gene,
      chr  = as.character(seqnames),
      start, end,
      lfc, padj,
      dist_probe = dist_probe,
      group = ifelse(dist_probe <= near_bp, "chr19 near probe (±1Mb)", "chr19 other")
    ),
  as.data.frame(gr_all[seqnames(gr_all) != "chr19"]) %>%
    transmute(
      gene = gene,
      chr  = as.character(seqnames),
      start, end,
      lfc, padj,
      dist_probe = NA_real_,
      group = "other chromosomes"
    )
)

### Prepare side by side LF + distance groups
codf_GB <- as.data.frame(resG_B) %>%
  mutate(gene = rownames(resG_B)) %>%
  select(gene, lfc_GB = log2FoldChange)

df_GminG <- as.data.frame(resGmin_G) %>%
  mutate(gene = rownames(resGmin_G)) %>%
  select(gene, lfc_GminG = log2FoldChange)

df_merge2 <- inner_join(df_GB, df_GminG, by = "gene")

df_merge2 <- df_merge %>%
  left_join(
    df_out %>% select(gene, group),
    by = "gene"
  ) %>%
  filter(!is.na(group))

cor.test(df_merge2$lfc_GB, df_merge2$lfc_GminG, method = "spearman")
#-0.1884573, p-value < 2.2e-16

cor_by_group <- df_merge2 %>%
  filter(!is.na(group)) %>%
  group_by(group) %>%
  summarise(
    n = n(),
    rho = cor(lfc_GB, lfc_GminG, method = "spearman"),
    p_value = cor.test(lfc_GB, lfc_GminG, method = "spearman")$p.value
  )

cor_by_group
#group                       n    rho   p_value
#1 chr19 near probe (±1Mb)    47 -0.268 6.93e-  2
#2 chr19 other              1044 -0.221 5.69e- 13
#3 other chromosomes       13443 -0.186 1.47e-104

### Plot

ggplot() +
  # background genomewide points
  geom_point(
    data = df_merge2 %>% 
      filter(group == "other chromosomes"),
    aes(x = lfc_GB, y = lfc_GminG),
    color = "grey70",
    alpha = 0.2,
    size = 1,
    shape = 16
  ) +
  # chr19 other
  geom_point(
    data = df_merge2 %>% 
      filter(group == "chr19 other"),
    aes(x = lfc_GB, y = lfc_GminG),
    color = "#A374B2",
    alpha = 1,
    size = 1,
    shape = 16
  ) +
  # chr19 near probe
  geom_point(
    data = df_merge2 %>% 
      filter(group == "chr19 near probe (±1Mb)"),
    aes(x = lfc_GB, y = lfc_GminG),
    color = "#52335E",
    alpha = 1,
    size = 1.5,
    shape = 16
  ) +
  # regression lines
  geom_smooth(
    data = df_merge2,
    aes(x = lfc_GB, y = lfc_GminG),
    method = "lm",
    color = "black",
    se = FALSE
  ) +
  theme_classic() +
  labs(
    x = "log2FC (GRANTA vs BLAS)",
    y = "log2FC (GRANTA+Min vs GRANTA)"
  )

######### Enrichment analysis of transcriptional reversal #####################
df_GB <- as.data.frame(resG_B) %>%
  tibble::rownames_to_column("gene") %>%
  select(gene,
         lfc_GB = log2FoldChange,
         padj_GB = padj)

df_GminG <- as.data.frame(resGmin_G) %>%
  tibble::rownames_to_column("gene") %>%
  select(gene,
         lfc_GminG = log2FoldChange,
         padj_GminG = padj)

df_all <- full_join(df_GB, df_GminG, by = "gene")

lfc_thr  <- 1
padj_max <- 0.05

df_all <- df_all %>%
  mutate(
    GB_up = !is.na(padj_GB) & !is.na(lfc_GB) & padj_GB < padj_max & lfc_GB >  lfc_thr,
    Min_down = !is.na(padj_GminG) & !is.na(lfc_GminG) & padj_GminG < padj_max & lfc_GminG < -lfc_thr
  )

tab <- df_all %>%
  filter(GB_up) %>%
  summarise(
    n_GB_up = n(),
    n_reversed = sum(Min_down),
    frac_reversed = n_reversed / n_GB_up
  )

tab

# Rates (as fractions)
rate_bg   <- mean(df_all$Min_down)
rate_gbup <- mean(df_all$Min_down[df_all$GB_up])

# Fisher enrichment
mat <- table(df_all$GB_up, df_all$Min_down)
ft <- fisher.test(mat)

or   <- unname(ft$estimate)
pval <- ft$p.value

label_txt <- paste0(
  "OR = ", sprintf("%.2f", or),
  "\nFisher p = ", formatC(pval, format = "e", digits = 2)
)

df_bar <- tibble::tibble(
  group = factor(c("Genome-wide", "GRANTA-up genes"),
                 levels = c("Genome-wide", "GRANTA-up genes")),
  frac  = c(rate_bg, rate_gbup)
)

p_bar <- ggplot(df_bar, aes(x = group, y = frac)) +
  geom_col(width = 0.7, fill = "grey70") +
  geom_text(aes(label = sprintf("%.1f%%", 100 * frac)),
            vjust = -0.4, size = 4) +
  annotate("text",
           x = 1.5,
           y = max(df_bar$frac) * 1.15,
           label = label_txt,
           size = 4) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.2))
  ) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Fraction of genes downregulated by Minnelide\n(GRANTA+Min vs GRANTA)"
  )

p_bar

################## Selected genes expression ##########################

vst <- vst(dds, blind = TRUE)
expr <- assay(vst)["CCND1", ]

df_ccnd1 <- data.frame(
  expression = as.numeric(expr),
  condition  = colData(dds)$cond
) %>% filter(cond == "GRANTA" | cond == "GRANTA_min")


ggplot(df_ccnd1, aes(x = condition, y = expression, fill = condition)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 1) +
  theme_classic() + coord_cartesian(ylim = c(13, 14.5)) +
  labs(
    x = NULL,
    y = "CCND1 expression (VST)"
  ) + theme(legend.position = "none")



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

#time zone: Europe/Riga
#tzcode source: internal

#attached base packages:
#  [1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] DOSE_3.26.2                 enrichplot_1.20.0           msigdbr_7.5.1               fgsea_1.26.0                clusterProfiler_4.8.2       factoextra_1.0.7           
#[7] FactoMineR_2.11             latticeExtra_0.6-30         lattice_0.22-6              gtable_0.3.5                gridExtra_2.3               ggupset_0.4.0              
#[13] ggplotify_0.1.2             org.Hs.eg.db_3.17.0         AnnotationDbi_1.62.2        VennDiagram_1.7.3           futile.logger_1.4.3         ggvenn_0.1.10              
#[19] dplyr_1.1.4                 ggrepel_0.9.5               ggplot2_3.5.1               karyoploteR_1.26.0          regioneR_1.32.0             DESeq2_1.40.2              
#[25] SummarizedExperiment_1.30.2 Biobase_2.60.0              MatrixGenerics_1.12.3       matrixStats_1.3.0           GenomicRanges_1.52.1        GenomeInfoDb_1.36.4        
#[31] IRanges_2.34.1              S4Vectors_0.38.2            BiocGenerics_0.46.0         biomaRt_2.56.1             
#
#loaded via a namespace (and not attached):
#  [1] fs_1.6.4                 ProtGenerics_1.32.0      bitops_1.0-7             HDO.db_0.99.1            httr_1.4.7               RColorBrewer_1.1-3      
#[7] tools_4.3.0              backports_1.5.0          utf8_1.2.4               R6_2.5.1                 DT_0.33                  lazyeval_0.2.2          
#[13] withr_3.0.0              prettyunits_1.2.0        cli_3.6.3                formatR_1.14             flashClust_1.01-2        scatterpie_0.2.3        
#[19] sandwich_3.1-0           labeling_0.4.3           mvtnorm_1.2-5            Rsamtools_2.16.0         yulab.utils_0.1.4        gson_0.1.0              
#[25] foreign_0.8-87           dichromat_2.0-1          BSgenome_1.68.0          rstudioapi_0.16.0        RSQLite_2.3.7            generics_0.1.3          
#[31] gridGraphics_0.5-1       BiocIO_1.10.0            car_3.1-2                GO.db_3.17.0             leaps_3.2                Matrix_1.6-5            
#[37] interp_1.1-6             fansi_1.0.6              abind_1.4-7              lifecycle_1.0.4          scatterplot3d_0.3-44     multcomp_1.4-26         
#[43] yaml_2.3.9               carData_3.0-5            qvalue_2.32.0            BiocFileCache_2.13.0     blob_1.2.4               crayon_1.5.3            
#[49] cowplot_1.1.3            GenomicFeatures_1.52.2   KEGGREST_1.40.1          pillar_1.9.0             knitr_1.48               rjson_0.2.21            
#[55] estimability_1.5.1       codetools_0.2-20         fastmatch_1.1-4          glue_1.7.0               downloader_0.4           ggfun_0.1.5             
#[61] data.table_1.15.4        vctrs_0.6.5              png_0.1-8                treeio_1.24.3            cachem_1.1.0             xfun_0.46               
#[67] S4Arrays_1.0.6           tidygraph_1.3.1          coda_0.19-4.1            survival_3.7-0           TH.data_1.1-2            nlme_3.1-165            
#[73] ggtree_3.8.2             bit64_4.0.5              progress_1.2.3           filelock_1.0.3           rpart_4.1.23             colorspace_2.1-1        
#[79] DBI_1.2.3                Hmisc_5.1-3              nnet_7.3-19              tidyselect_1.2.1         emmeans_1.10.2           bit_4.0.5               
#[85] compiler_4.3.0           curl_5.2.1               htmlTable_2.4.3          bezier_1.1.2             xml2_1.3.6               DelayedArray_0.26.7     
#[91] shadowtext_0.1.4         rtracklayer_1.60.1       checkmate_2.3.1          scales_1.3.0             multcompView_0.1-10      rappdirs_0.3.3          
#[97] stringr_1.5.1            digest_0.6.36            rmarkdown_2.27           XVector_0.40.0           htmltools_0.5.8.1        pkgconfig_2.0.3         
#[103] jpeg_0.1-10              base64enc_0.1-3          dbplyr_2.5.0             fastmap_1.2.0            ensembldb_2.24.1         rlang_1.1.4             
#[109] htmlwidgets_1.6.4        farver_2.1.2             zoo_1.8-12               jsonlite_1.8.8           BiocParallel_1.34.2      GOSemSim_2.26.1         
#[115] VariantAnnotation_1.46.0 RCurl_1.98-1.16          magrittr_2.0.3           Formula_1.2-6            GenomeInfoDbData_1.2.10  patchwork_1.2.0         
#[121] munsell_0.5.1            Rcpp_1.0.13              babelgene_22.9           ape_5.8                  bamsignals_1.32.0        viridis_0.6.5           
#[127] stringi_1.8.4            ggraph_2.2.1             zlibbioc_1.46.0          MASS_7.3-60              plyr_1.8.9               parallel_4.3.0          
#[133] deldir_2.0-4             Biostrings_2.68.1        graphlayouts_1.1.1       splines_4.3.0            hms_1.1.3                locfit_1.5-9.10         
#[139] igraph_2.0.3             ggpubr_0.6.0             ggsignif_0.6.4           reshape2_1.4.4           futile.options_1.0.1     XML_3.99-0.17           
#[145] evaluate_0.24.0          biovizBase_1.48.0        lambda.r_1.2.4           tweenr_2.0.3             tidyr_1.3.1              purrr_1.0.2             
#[151] polyclip_1.10-6          ggforce_0.4.2            broom_1.0.6              xtable_1.8-4             restfulr_0.0.15          AnnotationFilter_1.24.0 
#[157] tidytree_0.4.6           rstatix_0.7.2            viridisLite_0.4.2        tibble_3.2.1             aplot_0.2.3              memoise_2.0.1           
#[163] GenomicAlignments_1.36.0 cluster_2.1.6  

