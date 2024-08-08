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
library(biomaRt)
library(FactoMineR)
library(factoextra)
library(clusterProfiler)
library(fgsea)
library(msigdbr)
library(enrichplot)
library(DOSE)
library(karyoploteR)
library(EnsDb.Hsapiens.v86)

################### Loading the count table ##############################
counts = read.table("salmon.merged.gene_counts.tsv", header = T)
rownames(counts) <- counts$gene_id
counts <- counts[,-1]
counts <- counts[,-1]

######################### Quality checks ################################
par(mfrow=c(3,3))
#Histograms of log transformed counts to see reads distribution
hist(log(as.numeric(counts[,1])+1), main = "MCL_A_1", xlab = "log(counts+1)")
hist(log(as.numeric(counts[,2])+1), main = "MCL_A_2", xlab = "log(counts+1)")
hist(log(as.numeric(counts[,3])+1), main = "MCL_A_3", xlab = "log(counts+1)")

hist(log(as.numeric(counts[,4])+1), main = "MCL_c3_1", xlab = "log(counts+1)")
hist(log(as.numeric(counts[,5])+1), main = "MCL_c3_2", xlab = "log(counts+1)")
hist(log(as.numeric(counts[,6])+1), main = "MCL_c3_3", xlab = "log(counts+1)")

hist(log(as.numeric(counts[,7])+1), main = "MCL_c7_1", xlab = "log(counts+1)")
hist(log(as.numeric(counts[,8])+1), main = "MCL_c7_2", xlab = "log(counts+1)")
hist(log(as.numeric(counts[,9])+1), main = "MCL_c7_3", xlab = "log(counts+1)")

hist(log(as.numeric(counts[,10])+1), main = "MCL_25_1", xlab = "log(counts+1)")
hist(log(as.numeric(counts[,11])+1), main = "MCL_25_2", xlab = "log(counts+1)")
hist(log(as.numeric(counts[,12])+1), main = "MCL_25_3", xlab = "log(counts+1)")

hist(log(as.numeric(counts[,13])+1), main = "MCL_50_2", xlab = "log(counts+1)")
hist(log(as.numeric(counts[,14])+1), main = "MCL_50_3", xlab = "log(counts+1)")
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
gsub(".[0-9]", "", counts)
dim(counts)

dds <- DESeqDataSetFromMatrix(counts, colData = DataFrame(cond), design = ~ cond)
colData(dds)
dds <- DESeq(dds)
rld <- rlog(dds)
plotPCA(rld, intgroup="cond")

resM25_c <- results(dds, contrast = c("cond", "MCL_M25", "MCL_c3"))
resM50_c <- results(dds, contrast = c("cond", "MCL_M50", "MCL_c3"))
resA_c <- results(dds, contrast = c("cond", "MCL_M50", "MCL_c7"))

summary(resM25_c)
summary(resM50_c)
summary(resMA_c)

resM25_c <- resM25_c[complete.cases(resM25_c),]
resM25_c <- resM25_c[order(resM25_c$padj),] 
m25 <- getBM(attributes = c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol', 'genecards', 'ensembl_gene_id', 'description'),
             filters = 'hgnc_symbol', values = row.names(resM25_c),
             mart = ensembl)
resM25_c_df <- as.data.frame(resM25_c)
resM25_c_df$id <- row.names(resM25_c)
resM25_c_df <- left_join(resM25_c_df, M25, by = join_by("id" == "hgnc_symbol"))
write.csv2(resM25_c_df, "results/res25Min_3d.csv")

resM50_c <- resM50_c[complete.cases(resM50_c),]
resM50_c <- resM50_c[order(resM50_c$padj),] 
m50 <- getBM(attributes = c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol', 'genecards', 'ensembl_gene_id', 'description'),
             filters = 'hgnc_symbol', values = row.names(resM50_c),
             mart = ensembl)
resM50_c_df <- as.data.frame(resM50_c)
resM50_c_df$id <- row.names(resM50_c)
resM50_c_df <- left_join(resM50_c_df, M50, by = join_by("id" == "hgnc_symbol"))
write.csv2(resM50_c_df, "results/res50Min_3d.csv")

resA_c <- resA_c[complete.cases(resA_c),]
resA_c <- resA_c[order(resA_c$padj),] 
a <- getBM(attributes = c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol', 'genecards', 'ensembl_gene_id', 'description'),
             filters = 'hgnc_symbol', values = row.names(resA_c),
             mart = ensembl)
resA_c_df <- as.data.frame(resA_c)
resA_c_df$id <- row.names(resA_c)
resA_c_df <- left_join(resA_c_df, a, by = join_by("id" == "hgnc_symbol"))
write.csv2(resA_c_df , "results/resAbe_7d.csv")

###Filtering
resM25_c_sign <- resM25_c[resM25_c[,'padj'] < 0.05, ]
resM25_c_up <- resM25_c_sign[resM25_c_sign[,'log2FoldChange'] > 1, ]
resM25_c_down <- resM25_c_sign[resM25_c_sign[,'log2FoldChange'] < -1, ]
summary(resM25_c_up)
summary(resM25_c_down)

resM50_c_sign <- resM50_c[resM50_c[,'padj'] < 0.05, ]
resM50_c_up <- resM50_c_sign[resM50_c_sign[,'log2FoldChange'] > 1, ]
resM50_c_down <- resM50_c_sign[resM50_c_sign[,'log2FoldChange'] < -1, ]

resA_c_sign <- resA_c[resA_c[,'padj'] < 0.05, ]
resA_c_up <- resA_c_sign[resA_c_sign[,'log2FoldChange'] > 1, ]
resA_c_down <- resA_c_sign[resA_c_sign[,'log2FoldChange'] < -1, ]
summary(resA_c_up)
summary(resA_c_down)

################################### Volcanos ##################################
volc_M50_c <- ggplot(data=as.data.frame(resM50_c)) + geom_point(aes(x=log2FoldChange, y=-log10(padj)), col="grey") +
  geom_point(data=as.data.frame(resM50_c_up), aes(x=log2FoldChange, y=-log10(padj)), col="indianred") +
  geom_point(data=as.data.frame(resM50_c_down), aes(x=log2FoldChange, y=-log10(padj)), col="skyblue4") +
  theme_minimal() +
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5)) +
  ggtitle("MCL patient cells + Minnelide 50 vs MCL patient cells")

volc_M25_c <- ggplot(data=as.data.frame(resM25_c)) + geom_point(aes(x=log2FoldChange, y=-log10(padj)), col="grey") +
  geom_point(data=as.data.frame(resM25_c_up), aes(x=log2FoldChange, y=-log10(padj)), col="indianred") +
  geom_point(data=as.data.frame(resM25_c_down), aes(x=log2FoldChange, y=-log10(padj)), col="skyblue4") +
  theme_minimal() +
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5)) +
  ggtitle("MCL patient cells + Minnelide 25 vs MCL patient cells")

volc_A_c <- ggplot(data=as.data.frame(resA_c)) + geom_point(aes(x=log2FoldChange, y=-log10(padj)), col="grey") +
  geom_point(data=as.data.frame(resA_c_up), aes(x=log2FoldChange, y=-log10(padj)), col="indianred") +
  geom_point(data=as.data.frame(resA_c_down), aes(x=log2FoldChange, y=-log10(padj)), col="skyblue4") +
  theme_minimal() +
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5)) +
  ggtitle("MCL patient cells + Abemaciclib vs MCL patient cells")

################################### Venns ##################################
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

#Min50
intersection_mcl_enh_min50 <- list(A= mcl_up$Gene_Name,
                                   B=genes_regulated_by_ehnancers$hgnc_symbol,
                                   C=rownames(resM50_c_down))
names(intersection_mcl_enh_min50) <- c("MCL_DAG_vs_naive_up",
                                        "regulated by enhancers, patients",
                                        "down by Min50, patient")

venn50 <- ggvenn(intersection_mcl_enh_min50,
               stroke_size = 0.5, set_name_size = 4,
               show_percentage = FALSE)

#Min25
intersection_mcl_enh_min25 <- list(A=mcl_up$Gene_Name,
                                   B=genes_regulated_by_ehnancers$hgnc_symbol,
                                   C=rownames(resM25_c_down))
names(intersection_mcl_enh_min25) <- c("MCL_DAG_vs_naive_up",
                                       "regulated by enhancers, patients",
                                       "down by Min25, patient")

venn25 <- ggvenn(intersection_mcl_enh_min25,
               stroke_size = 0.5, set_name_size = 4,
               show_percentage = FALSE)

#Abe
intersection_mcl_enh_abe <- list(A=mcl_up$Gene_Name,
                                 B=genes_regulated_by_ehnancers$hgnc_symbol,
                                 C=rownames(resA_c_down))
names(intersection_mcl_enh_abe) <- c("MCL_DAG_vs_naive_up",
                                       "regulated by enhancers, patients",
                                       "down by Abe, patient")

vennAbe <- ggvenn(intersection_mcl_enh_abe,
                 stroke_size = 0.5, set_name_size = 4,
                 show_percentage = FALSE)

############################## Enrichment analysis ################################
### GSEA ###
abe_gene_list <- resA_c_sign$log2FoldChange
names(abe_gene_list) <- row.names(resA_c_sign)
abe_gene_list = sort(abe_gene_list, decreasing = TRUE)

gse_abe_bp <- gseGO(geneList=abe_gene_list, 
                    ont ="BP", 
                    keyType = "SYMBOL", 
                    nPerm = 10000, 
                    minGSSize = 3, 
                    maxGSSize = 800, 
                    pvalueCutoff = "none", 
                    verbose = TRUE, 
                    OrgDb = org.Hs.eg.db, 
                    pAdjustMethod = "BH")
write.csv2(gse_abe_bp, "results/gse_abe_bp.csv")

gse_abe_bp <- pairwise_termsim(gse_abe_bp)
treeplot(gse_abe_bp, showCategory = 40, color = 'NES')
cnetplot(gse_abe_bp, foldChange=abe_gene_list)
cnetplot(gse_abe_bp, foldChange=abe_gene_list, circular = TRUE, colorEdge = TRUE)

min25_gene_list <- resM25_c_sign$log2FoldChange
names(min25_gene_list) <- row.names(resM25_c_sign)
min25_gene_list= sort(min25_gene_list, decreasing = TRUE)
gse_min25_bp <- gseGO(geneList=min25_gene_list, 
                    ont ="BP", 
                    keyType = "SYMBOL", 
                    nPerm = 10000, 
                    minGSSize = 3, 
                    maxGSSize = 800, 
                    pvalueCutoff = "none", 
                    verbose = TRUE, 
                    OrgDb = org.Hs.eg.db, 
                    pAdjustMethod = "BH")
write.csv2(gse_min25_bp, "results/gse_min25_bp.csv")

gse_min25_bp <- pairwise_termsim(gse_min25_bp)
treeplot(gse_min25_bp, showCategory = 40, color = 'NES')
cnetplot(gse_min25_bp, foldChange=min25_gene_list)
cnetplot(gse_min25_bp, foldChange=min25_gene_list, circular = TRUE, colorEdge = TRUE)

min50_gene_list <- resM50_c_sign$log2FoldChange
names(min50_gene_list) <- row.names(resM50_c_sign)
min50_gene_list= sort(min50_gene_list, decreasing = TRUE)
gse_min50_bp <- gseGO(geneList=min50_gene_list, 
                      ont ="BP", 
                      keyType = "SYMBOL", 
                      nPerm = 10000, 
                      minGSSize = 3, 
                      maxGSSize = 800, 
                      pvalueCutoff = "none", 
                      verbose = TRUE, 
                      OrgDb = org.Hs.eg.db, 
                      pAdjustMethod = "BH")
write.csv2(gse_min50_bp, "results/gse_min50_bp.csv")

gse_min50_bp <- pairwise_termsim(gse_min50_bp)
treeplot(gse_min50_bp, showCategory = 40, color = 'NES')
cnetplot(gse_min50_bp, foldChange=min50_gene_list)
cnetplot(gse_min50_bp, foldChange=min50_gene_list, circular = TRUE, colorEdge = TRUE)

### GO ###
listA <- list(up = rownames(resA_c_up), down = rownames(resA_c_down))
list25 <- list(up = rownames(resM25_c_up), down = rownames(resM25_c_down))
list50 <- list(up = rownames(resM50_c_up), down = rownames(resM50_c_down))

compGO_A <- compareCluster(geneCluster   = listA,
                          fun           = "enrichGO",
                          pvalueCutoff  = 0.05,
                          keyType       = "SYMBOL",
                          pAdjustMethod = "BH",
                          OrgDb='org.Hs.eg.db',
                          ont           = "BP",
                          universe = rownames(counts))
write.csv2(compGO_A, "./results/go_enrichment_BP_Abema.csv")
dotplot(compGO_A, showCategory = 10, title = "go_enrichment_BP_Abema")


compGO_25 <- compareCluster(geneCluster   = list25,
                           fun           = "enrichGO",
                           pvalueCutoff  = 0.05,
                           keyType       = "SYMBOL",
                           pAdjustMethod = "BH",
                           OrgDb='org.Hs.eg.db',
                           ont           = "BP",
                           universe = rownames(counts))
write.csv2(compGO_25, "./results/go_enrichment_BP_Min25.csv")
dotplot(compGO_25, showCategory = 10, title = "go_enrichment_BP_Min25")


compGO_50 <- compareCluster(geneCluster   = list50,
                            fun           = "enrichGO",
                            pvalueCutoff  = 0.05,
                            keyType       = "SYMBOL",
                            pAdjustMethod = "BH",
                            OrgDb='org.Hs.eg.db',
                            ont           = "BP",
                            universe = rownames(counts))
write.csv2(compGO_50, "./results/go_enrichment_BP_Min50.csv")
dotplot(compGO_50, showCategory = 10, title = "go_enrichment_BP_Min50")

compGO_A_MF <- compareCluster(geneCluster   = listA,
                           fun           = "enrichGO",
                           pvalueCutoff  = 0.05,
                           keyType       = "SYMBOL",
                           pAdjustMethod = "BH",
                           OrgDb='org.Hs.eg.db',
                           ont           = "MF",
                           universe = rownames(counts))
write.csv2(compGO_A_MF, "./results/go_enrichment_MF_Abema.csv")
dotplot(compGO_A_MF, showCategory = 10, title = "go_enrichment_MF_Abema")

compGO_25_MF <- compareCluster(geneCluster   = list25,
                            fun           = "enrichGO",
                            pvalueCutoff  = 0.05,
                            keyType       = "SYMBOL",
                            pAdjustMethod = "BH",
                            OrgDb='org.Hs.eg.db',
                            ont           = "MF",
                            universe = rownames(counts))
write.csv2(compGO_25_MF, "./results/go_enrichment_MF_Min25.csv")
dotplot(compGO_25_MF, showCategory = 10, title = "go_enrichment_MF_Min25")

### KEGG ###
#converting to entrez
universe_entrez = bitr(rownames(counts), 
                       OrgDb = org.Hs.eg.db, 
                       fromType = "SYMBOL", 
                       toType = "ENTREZID")
universe_entrez <- universe_entrez$ENTREZID

# Abemaciclib
resA_c_up_entrez = bitr(rownames(resA_c_up), 
                        OrgDb = org.Hs.eg.db, 
                        fromType = "SYMBOL", 
                        toType = "ENTREZID")
resA_c_up_entrez <- resA_c_up_entrez$ENTREZID
resA_c_down_entrez = bitr(rownames(resA_c_down), 
                          OrgDb = org.Hs.eg.db, 
                          fromType = "SYMBOL", 
                          toType = "ENTREZID")
resA_c_down_entrez <- resA_c_down_entrez$ENTREZID
listA_entrez <- list(up = resA_c_up_entrez, down = resA_c_down_entrez)

compKEGG_A <- compareCluster(geneCluster   = listA_entrez,
                           fun           = "enrichKEGG",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH",
                           universe = universe_entrez)
write.csv2(compGO_A, "./results/kegg_enrichment_Abema.csv")
dotplot(compKEGG_A, showCategory = 10, title = "kegg_enrichment_BP_Abema")

# Minnelide 25
res25_c_up_entrez = bitr(rownames(resM25_c_up), 
                        OrgDb = org.Hs.eg.db, 
                        fromType = "SYMBOL", 
                        toType = "ENTREZID")
res25_c_up_entrez <- res25_c_up_entrez$ENTREZID
res25_c_down_entrez = bitr(rownames(resM25_c_down), 
                          OrgDb = org.Hs.eg.db, 
                          fromType = "SYMBOL", 
                          toType = "ENTREZID")
res25_c_down_entrez <- res25_c_down_entrez$ENTREZID
list25_entrez <- list(up = res25_c_up_entrez, down = res25_c_down_entrez)

compKEGG_25 <- compareCluster(geneCluster   = list25_entrez,
                             fun           = "enrichKEGG",
                             pvalueCutoff  = 0.05,
                             pAdjustMethod = "BH",
                             universe = universe_entrez)
write.csv2(compKEGG_25, "./results/kegg_enrichment_Min25.csv")
dotplot(compKEGG_25, showCategory = 10, title = "kegg_enrichment_Min25")

# Minnelide 50
res50_c_up_entrez = bitr(rownames(resM50_c_up), 
                         OrgDb = org.Hs.eg.db, 
                         fromType = "SYMBOL", 
                         toType = "ENTREZID")
res50_c_up_entrez <- res50_c_up_entrez$ENTREZID
res50_c_down_entrez = bitr(rownames(resM50_c_down), 
                           OrgDb = org.Hs.eg.db, 
                           fromType = "SYMBOL", 
                           toType = "ENTREZID")
res50_c_down_entrez <- res50_c_down_entrez$ENTREZID
list50_entrez <- list(up = res50_c_up_entrez, down = res50_c_down_entrez)

compKEGG_50 <- compareCluster(geneCluster   = list50_entrez,
                              fun           = "enrichKEGG",
                              pvalueCutoff  = 0.05,
                              pAdjustMethod = "BH",
                              universe = universe_entrez)
write.csv2(compKEGG_50, "./results/kegg_enrichment_Min50.csv")
dotplot(compKEGG_50, showCategory = 10, title = "kegg_enrichment_Min50")

### Disease ontology enrichment 
compDO_50 <- compareCluster(geneCluster   = list50_entrez,
                              fun           = "enrichDO",
                              pvalueCutoff  = 0.05,
                              pAdjustMethod = "BH",
                              universe = universe_entrez,
                              readable = TRUE)
write.csv2(compDO_50, "./results/DO_enrichment_Min50.csv")
dotplot(compDO_50, showCategory = 10, title = "DO_enrichment_Min50")

compDO_25 <- compareCluster(geneCluster   = list25_entrez,
                            fun           = "enrichDO",
                            pvalueCutoff  = 0.05,
                            pAdjustMethod = "BH",
                            universe = universe_entrez,
                            readable = TRUE)
write.csv2(compDO_25, "./results/DO_enrichment_Min25.csv")
dotplot(compDO_25, showCategory = 10, title = "DO_enrichment_Min25")

compDO_A <- compareCluster(geneCluster   = listA_entrez,
                            fun           = "enrichDO",
                            pvalueCutoff  = 0.05,
                            pAdjustMethod = "BH",
                            universe = universe_entrez,
                            readable = TRUE)
write.csv2(compDO_A, "./results/DO_enrichment_A.csv")
dotplot(compDO_A, showCategory = 10, title = "DO_enrichment_A")

################## Enrichment analysis on intersection lists ###################
#Genes down by min and up in MCL GO
compGO_int1 <- enrichGO(intersect(rownames(resM25_c_down), rownames(mcl_up)),
                        pvalueCutoff  = 0.05,
                        keyType       = "SYMBOL",
                        pAdjustMethod = "BH",
                        OrgDb='org.Hs.eg.db',
                        ont           = "BP",
                        universe = rownames(counts))
write.csv2(compGO_int1, "./results/go_enrichment_BP_up_in_MCL_down_by_min.csv")
dotplot(compGO_int1, showCategory = 10, title = "go_enrichment_BP_up_in_MCL_down_by_min")

#Genes down by Abe and up in MCL GO
compGO_int3 <- enrichGO(intersect(rownames(resA_c_down), rownames(mcl_up)),
                        pvalueCutoff  = 0.05,
                        keyType       = "SYMBOL",
                        pAdjustMethod = "BH",
                        OrgDb='org.Hs.eg.db',
                        ont           = "BP",
                        universe = rownames(counts))
write.csv2(compGO_int3, "./results/go_enrichment_BP_up_in_MCL_down_by_abe.csv")
dotplot(compGO_int3, showCategory = 10, title = "go_enrichment_BP_up_in_MCL_down_by_abe")

#Genes down by min and regulated by enhancers GO
compGO_int2 <- enrichGO(intersect(rownames(resM25_c_down), genes_regulated_by_ehnancers$hgnc_symbol),
                        pvalueCutoff  = 0.05,
                        keyType       = "SYMBOL",
                        pAdjustMethod = "BH",
                        OrgDb='org.Hs.eg.db',
                        ont           = "BP",
                        universe = rownames(counts))
write.csv2(compGO_int2, "./results/go_enrichment_BP_reg_by_enh_down_by_min.csv")
dotplot(compGO_int2, showCategory = 10, title = "go_enrichment_BP_reg_by_enh_down_by_min")

#Genes down by Abe and regulated by enhancers GO
compGO_int4 <- enrichGO(intersect(rownames(resA_c_down), genes_regulated_by_ehnancers$hgnc_symbol),
                        pvalueCutoff  = 0.05,
                        keyType       = "SYMBOL",
                        pAdjustMethod = "BH",
                        OrgDb='org.Hs.eg.db',
                        ont           = "BP",
                        universe = rownames(counts))
write.csv2(compGO_int4, "./results/go_enrichment_BP_reg_by_enh_down_by_abe.csv")
dotplot(compGO_int4, showCategory = 10, title = "go_enrichment_BP_reg_by_enh_down_by_abe")

#Genes down by min and up in MCL KEGG
int1_entrez = bitr(intersect(rownames(resM25_c_down), rownames(mcl_up)), 
                   OrgDb = org.Hs.eg.db, 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID")
int1_entrez  <- int1_entrez$ENTREZID

compKEGG_int1 <- enrichKEGG(int1_entrez,
                            pvalueCutoff  = 0.05,
                            pAdjustMethod = "BH",
                            universe = universe_entrez)
write.csv2(compKEGG_int1, "./results/kegg_enrichment_up_in_MCL_down_by_min.csv")
dotplot(compKEGG_int1, showCategory = 10, title = "kegg_enrichment_up_in_MCL_down_by_min")

#Genes down by abe and up in MCL KEGG
int2_entrez = bitr(intersect(rownames(resA_c_down), rownames(mcl_up)), 
                   OrgDb = org.Hs.eg.db, 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID")
int2_entrez  <- int2_entrez$ENTREZID

compKEGG_int2 <- enrichKEGG(int2_entrez,
                            pvalueCutoff  = 0.05,
                            pAdjustMethod = "BH",
                            universe = universe_entrez)
write.csv2(compKEGG_int2, "./results/kegg_enrichment_up_in_MCL_down_by_abe.csv")
dotplot(compKEGG_int2, showCategory = 10, title = "kegg_enrichment_up_in_MCL_down_by_abe")

#Genes down by min and regulated by enhancers KEGG
int3_entrez = bitr(intersect(rownames(resM25_c_down), genes_regulated_by_ehnancers$hgnc_symbol), 
                   OrgDb = org.Hs.eg.db, 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID")
int3_entrez  <- int3_entrez$ENTREZID

compKEGG_int3 <- enrichKEGG(int3_entrez,
                        pvalueCutoff  = 0.05,
                        pAdjustMethod = "BH",
                        universe = universe_entrez)
write.csv2(compKEGG_int3, "./results/go_enrichment_kegg_reg_by_enh_down_by_min.csv")
dotplot(compKEGG_int3, showCategory = 10, title = "go_enrichment_kegg_reg_by_enh_down_by_min")

#Genes down by abe and regulated by enhancers KEGG
int4_entrez = bitr(intersect(rownames(resA_c_down), genes_regulated_by_ehnancers$hgnc_symbol), 
                   OrgDb = org.Hs.eg.db, 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID")
int4_entrez  <- int4_entrez$ENTREZID

compKEGG_int4 <- enrichKEGG(int4_entrez,
                            pvalueCutoff  = 0.05,
                            pAdjustMethod = "BH",
                            universe = universe_entrez)
write.csv2(compKEGG_int4, "./results/go_enrichment_kegg_reg_by_enh_down_by_abe.csv")
dotplot(compKEGG_int4, showCategory = 10, title = "go_enrichment_kegg_reg_by_enh_down_by_abe")

############## Barplots DEGs per chromosome ########################
### DEGs per chromosome, MCL
chromsizes <- read.csv2("./chromSizes.csv", header = TRUE)
chromsizes$normsizes <- chromsizes$size/min(chromsizes$size)
chromsizes$normnumber <- chromsizes$gene_number/min(chromsizes$gene_number)
chromsizes$chr <- gsub("chr", "", chromsizes$chr)

##  all DEGs Abema
df_a <- as.data.frame(table(resA_c_df$chromosome_name))
df_a <- df_a[1:22,]
colnames(df_a) <- c("chr", "freq")
df_a <- merge(df_a,chromsizes, by.x = "chr", by.y = "chr")
df_a$normfreq_size <- df_a$freq/df_a$normsizes #norm by size
df_a$normfreq_number <- df_a$freq/df_a$normnumber #norm by gene number
write.csv(df_a, "results/degs_per_chromosome_abe.csv")

##  all DEGs Min 25
df_25 <- as.data.frame(table(resM25_c_df$chromosome_name))
df_25 <- df_25[1:22,]
colnames(df_25) <- c("chr", "freq")
df_25 <- merge(df_25,chromsizes, by.x = "chr", by.y = "chr")
df_25$normfreq_size <- df_25$freq/df_25$normsizes #norm by size
df_25$normfreq_number <- df_25$freq/df_25$normnumber #norm by gene number
write.csv(df_25, "results/degs_per_chromosome_min25.csv")

##  all DEGs Min 50
df_50 <- as.data.frame(table(resM50_c_df$chromosome_name))
df_50 <- df_50[1:22,]
colnames(df_50) <- c("chr", "freq")
df_50 <- merge(df_50,chromsizes, by.x = "chr", by.y = "chr")
df_50$normfreq_size <- df_50$freq/df_50$normsizes #norm by size
df_50$normfreq_number <- df_50$freq/df_50$normnumber #norm by gene number
write.csv(df_50, "results/degs_per_chromosome_min50.csv")

## down DEGs Abe 
df_a_down <- as.data.frame(table(resA_c_df[resA_c_df$log2FoldChange < -1 & resA_c_df$padj < 0.05,]$chromosome_name))
df_a_down <- df_a_down[1:22,]
colnames(df_a_down) <- c("chr", "freq")
df_a_down <- merge(df_a_down,chromsizes, by.x = "chr", by.y = "chr")
df_a_down$normfreq_size <- df_a_down$freq/df_a_down$normsizes #norm by size
df_a_down$normfreq_number <- df_a_down$freq/df_a_down$normnumber #norm by gene number
write.csv(df_a_down, "results/down_degs_per_chromosome_abe.csv")

## down DEGs Min25
df_25_down <- as.data.frame(table(resM25_c_df[resM25_c_df$log2FoldChange < -1 & resM25_c_df$padj < 0.05,]$chromosome_name))
df_25_down <- df_25_down[1:22,]
colnames(df_25_down) <- c("chr", "freq")
df_25_down <- merge(df_25_down,chromsizes, by.x = "chr", by.y = "chr")
df_25_down$normfreq_size <- df_25_down$freq/df_25_down$normsizes #norm by size
df_25_down$normfreq_number <- df_25_down$freq/df_25_down$normnumber #norm by gene number
write.csv(df_25_down, "results/down_degs_per_chromosome_min25.csv")

## down DEGs Min50
df_50_down <- as.data.frame(table(resM50_c_df[resM50_c_df$log2FoldChange < -1 & resM50_c_df$padj < 0.05,]$chromosome_name))
df_50_down <- df_50_down[1:22,]
colnames(df_50_down) <- c("chr", "freq")
df_50_down <- merge(df_50_down,chromsizes, by.x = "chr", by.y = "chr")
df_50_down$normfreq_size <- df_50_down$freq/df_50_down$normsizes #norm by size
df_50_down$normfreq_number <- df_50_down$freq/df_50_down$normnumber #norm by gene number
write.csv(df_50_down, "results/down_degs_per_chromosome_min50.csv")

df_a_down <- df_a_down[order(df_a_down$freq,decreasing=TRUE),]
barplot(df_a_down$freq, names.arg = df_a_down$chr, ylab = "N down DEGs Abe")
df_a_down <- df_a_down[order(df_a_down$normfreq_size,decreasing=TRUE),]
barplot(df_a_down$normfreq_size, names.arg = df_a_down$chr, ylab = "N down DEGs Abe normalised to chromSize patients")
df_a_down <- df_a_down[order(df_a_down$normfreq_number,decreasing=TRUE),]
barplot(df_a_down$normfreq_number, names.arg = df_a_down$chr, ylab = "N down DEGs Abe normalised to gene number patients")

df_25_down <- df_25_down[order(df_25_down$freq,decreasing=TRUE),]
barplot(df_25_down$freq, names.arg = df_25_down$chr, ylab = "N down DEGs Min25")
df_25_down <- df_25_down[order(df_25_down$normfreq_size,decreasing=TRUE),]
barplot(df_25_down$normfreq_size, names.arg = df_25_down$chr, ylab = "N down DEGs Min25 normalised to chromSize patients")
df_25_down <- df_25_down[order(df_25_down$normfreq_number,decreasing=TRUE),]
barplot(df_25_down$normfreq_number, names.arg = df_25_down$chr, ylab = "N down DEGs Min25 normalised to gene number patients")

df_50_down <- df_50_down[order(df_50_down$freq,decreasing=TRUE),]
barplot(df_50_down$freq, names.arg = df_50_down$chr, ylab = "N down DEGs Min50") 
df_50_down <- df_50_down[order(df_50_down$normfreq_size,decreasing=TRUE),]
barplot(df_50_down$normfreq_size, names.arg = df_50_down$chr, ylab = "N down DEGs Min50 normalised to chromSize patients")
df_50_down <- df_50_down[order(df_50_down$normfreq_number,decreasing=TRUE),]
barplot(df_50_down$normfreq_number, names.arg = df_50_down$chr, ylab = "N down DEGs Min50 normalised to gene number patients")

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
#  [1] DOSE_3.26.2                 enrichplot_1.20.0           msigdbr_7.5.1               fgsea_1.26.0               
#[5] clusterProfiler_4.8.2       EnsDb.Hsapiens.v86_2.99.0   ensembldb_2.24.1            AnnotationFilter_1.24.0    
#[9] GenomicFeatures_1.52.2      biomaRt_2.56.1              ggupset_0.4.0               ggplotify_0.1.2            
#[13] org.Hs.eg.db_3.17.0         AnnotationDbi_1.62.2        VennDiagram_1.7.3           futile.logger_1.4.3        
#[17] ggvenn_0.1.10               dplyr_1.1.4                 ggrepel_0.9.5               ggplot2_3.5.1              
#[21] karyoploteR_1.26.0          regioneR_1.32.0             DESeq2_1.40.2               SummarizedExperiment_1.30.2
#[25] Biobase_2.60.0              MatrixGenerics_1.12.3       matrixStats_1.3.0           GenomicRanges_1.52.1       
#[29] GenomeInfoDb_1.36.4         IRanges_2.34.1              S4Vectors_0.38.2            BiocGenerics_0.46.0        

#loaded via a namespace (and not attached):
#  [1] splines_4.3.0            BiocIO_1.10.0            bitops_1.0-7             filelock_1.0.3           tibble_3.2.1            
#[6] polyclip_1.10-6          XML_3.99-0.17            rpart_4.1.23             lifecycle_1.0.4          MASS_7.3-60             
#[11] lattice_0.22-6           backports_1.5.0          magrittr_2.0.3           Hmisc_5.1-3              rmarkdown_2.27          
#[16] yaml_2.3.9               cowplot_1.1.3            DBI_1.2.3                RColorBrewer_1.1-3       abind_1.4-7             
#[21] zlibbioc_1.46.0          purrr_1.0.2              ggraph_2.2.1             biovizBase_1.48.0        RCurl_1.98-1.16         
#[26] yulab.utils_0.1.4        nnet_7.3-19              tweenr_2.0.3             VariantAnnotation_1.46.0 rappdirs_0.3.3          
#[31] GenomeInfoDbData_1.2.10  tidytree_0.4.6           codetools_0.2-20         DelayedArray_0.26.7      ggforce_0.4.2           
#[36] xml2_1.3.6               tidyselect_1.2.1         aplot_0.2.3              farver_2.1.2             viridis_0.6.5           
#[41] BiocFileCache_2.13.0     base64enc_0.1-3          bamsignals_1.32.0        GenomicAlignments_1.36.0 jsonlite_1.8.8          
#[46] tidygraph_1.3.1          Formula_1.2-6            ggnewscale_0.5.0         tools_4.3.0              progress_1.2.3          
#[51] treeio_1.24.3            Rcpp_1.0.13              glue_1.7.0               gridExtra_2.3            xfun_0.46               
#[56] qvalue_2.32.0            withr_3.0.0              formatR_1.14             fastmap_1.2.0            fansi_1.0.6             
#[61] digest_0.6.36            R6_2.5.1                 gridGraphics_0.5-1       colorspace_2.1-1         GO.db_3.17.0            
#[66] dichromat_2.0-1          RSQLite_2.3.7            tidyr_1.3.1              utf8_1.2.4               generics_0.1.3          
#[71] data.table_1.15.4        rtracklayer_1.60.1       graphlayouts_1.1.1       prettyunits_1.2.0        httr_1.4.7              
#[76] htmlwidgets_1.6.4        S4Arrays_1.0.6           scatterpie_0.2.3         pkgconfig_2.0.3          gtable_0.3.5            
#[81] blob_1.2.4               XVector_0.40.0           shadowtext_0.1.4         htmltools_0.5.8.1        ProtGenerics_1.32.0     
#[86] scales_1.3.0             png_0.1-8                ggfun_0.1.5              knitr_1.48               lambda.r_1.2.4          
#[91] rstudioapi_0.16.0        reshape2_1.4.4           rjson_0.2.21             nlme_3.1-165             checkmate_2.3.1         
#[96] curl_5.2.1               cachem_1.1.0             stringr_1.5.1            parallel_4.3.0           HDO.db_0.99.1           
#[101] foreign_0.8-87           restfulr_0.0.15          pillar_1.9.0             vctrs_0.6.5              dbplyr_2.5.0            
#[106] cluster_2.1.6            htmlTable_2.4.3          evaluate_0.24.0          cli_3.6.3                locfit_1.5-9.10         
#[111] compiler_4.3.0           bezier_1.1.2             futile.options_1.0.1     Rsamtools_2.16.0         rlang_1.1.4             
#[116] crayon_1.5.3             labeling_0.4.3           plyr_1.8.9               fs_1.6.4                 stringi_1.8.4           
#[121] viridisLite_0.4.2        BiocParallel_1.34.2      babelgene_22.9           munsell_0.5.1            Biostrings_2.68.1       
#[126] lazyeval_0.2.2           GOSemSim_2.26.1          Matrix_1.6-5             BSgenome_1.68.0          patchwork_1.2.0         
#[131] hms_1.1.3                bit64_4.0.5              KEGGREST_1.40.1          igraph_2.0.3             memoise_2.0.1           
#[136] ggtree_3.8.2             fastmatch_1.1-4          bit_4.0.5                downloader_0.4           gson_0.1.0              
#[141] ape_5.8  
