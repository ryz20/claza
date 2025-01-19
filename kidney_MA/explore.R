set.seed(0)
library(tidyverse)
library(ggrepel)
library(DESeq2)

load("tidied_data.RData")

# qc
boxplot <- countTable %>% 
  pivot_longer(names(countTable), names_to = "sample", values_to = "count")
ggplot(boxplot, aes(x = sample, y = count)) +
  geom_boxplot()
ggplot(boxplot, aes(x = reorder(sample, count, FUN = sum), y = count)) +
  geom_boxplot()
ggplot(boxplot, aes(x = reorder(sample, count, FUN = median), y = count)) +
  geom_boxplot()

# tissue composition
tmp <- as.data.frame(t(countTable[row.names(countTable) %in% c("ABCA13", "CRYAA"),]))
tmp$logratio <- log2(2^(tmp$ABCA13)/2^(tmp$CRYAA))
ggplot(tmp, aes(x=row.names(tmp), y=(logratio)))+
  geom_point()+
  geom_hline(yintercept = 0)
coldata <- merge(coldata, tmp, by=0)
row.names(coldata) <- coldata$Row.names
coldata$Row.names <- NULL
coldata <- coldata[,!colnames(coldata)%in%c("ABCA13", "CRYAA")] # remove ABCA13 and CRYAA cols
coldata$CvMlogRatio <- coldata$logratio
coldata$logratio <- NULL
genelist <- read.csv(file = "cortexmedulla_genelist.csv", header = TRUE, sep = ",")
genelist$ratio <- genelist$cortex_Median/genelist$medulla_Median
genelist <- genelist %>% 
  arrange(ratio)
genelist <- genelist[c(1:7, 48:54),]
sigGenes <- genelist$X
sigGenes <- sigGenes[sigGenes%in%rownames(countTable)]
plotDat <- countTable[sigGenes,]
anno_info <- coldata[,c("CvMlogRatio","sample")]
plotDat <- plotDat[,match(row.names(anno_info), colnames(plotDat))]
anno_info$CvMlogRatio <- NULL
summary(colnames(plotDat)==row.names(anno_info))
pheatmap::pheatmap(mat = plotDat,
                   scale="row", 
                   cluster_rows = FALSE, 
                   cluster_cols = TRUE, 
                   clustering_distance_cols="euclidean",
                   show_rownames = TRUE,
                   show_colnames = FALSE,
                   annotation_col = anno_info)
colnames(countTable[,out$tree_col[["order"]]])
tmp <- as.data.frame(sort(cutree(out$tree_col, k=2)))
tmp <- rownames_to_column(tmp, var = "sample")
tmp <- tmp[match(row.names(coldata), tmp$sample),]
tmp <- mutate(tmp, Cluster = tmp[,2])
tmp <-  tmp[, c(1, 3)]
coldata <- merge(coldata, tmp, by="sample")
coldata[,'Cluster']<-factor(coldata[,'Cluster'])
coldata <- droplevels(coldata)
ggplot(coldata, aes(x = Cluster, y = CvMlogRatio))+
  geom_jitter(aes(colour = Cluster))
coldata <- coldata %>% 
  mutate(Tissue = ifelse(Cluster == 1,
                         "cortex",
                         "medulla"))
rownames(coldata) <- coldata$Sample


# pca
library(PCAtools)
coldata <- cbind(coldata,
                 WGCNA::binarizeCategoricalVariable(coldata$biopsy,
                                                    includePairwise = TRUE,
                                                    includeLevelVsAll = FALSE))
row.names(coldata) <- coldata$sample
data <- as.matrix(countTable)
data <- data[order(rowVars(data), decreasing = TRUE),]
data <- data[c(1:500), ]
p <- pca(data, metadata = coldata)
screeplot(p, axisLabSize = 18, titleLabSize = 22)
biplot(p, showLoadings = TRUE,
       labSize = 5, pointSize = 5, sizeLoadingsNames = 5)
pairsplot(p,
          components = getComponents(p, seq_len(5)))
eigencorplot(p,
             metavars = c("sex", "week12.vs.pre",
                          "week52.vs.pre",
                          "week52.vs.week12", 
                          "batch", "DAGE", "CvMlogRatio", "Tissue"))
summary(row.names(coldata) == row.names(p[["rotated"]]))
pc_scores <- data.frame(p[["rotated"]], coldata[row.names(coldata) %in% row.names(p[["rotated"]]),])
PC1v2 <- ggplot(data = pc_scores, aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, color = 'grey70') +
  geom_vline(xintercept = 0, color = 'grey70') +
  xlab(paste("PC1 (", round(p[["variance"]][1], 2), "%)", sep = "")) +
  ylab(paste("PC2 (", round(p[["variance"]][2], 2), "%)", sep = "")) +
  ggtitle("PCA")
library(RColorBrewer)
PC1v2 + geom_point(aes(color = sex), size = 4)
PC1v2 + geom_point(aes(color = biopsy), size = 4)
PC1v2 + geom_point(aes(color = person), size = 4) +
  geom_line(aes(group = person), color = "black", linetype="dotted")
PC1v2 + geom_point(aes(color = years_to_visit1), size = 3)
PC1v2 + geom_point(aes(color = DAGE), size = 4)
PC1v2 + geom_point(aes(color = center), size = 4)
PC1v2 + geom_point(aes(color = Tissue), size = 4)
PC1v2 + geom_point(aes(color = CvMlogRatio), size = 4)
PC3v4 <- ggplot(data = pc_scores, aes(x = PC3, y = PC4)) +
  geom_hline(yintercept = 0, color = 'grey70') +
  geom_vline(xintercept = 0, color = 'grey70') +
  xlab(paste("PC3 (", round(p[["variance"]][3], 2), "%)", sep = "")) +
  ylab(paste("PC4 (", round(p[["variance"]][4], 2), "%)", sep = "")) +
  ggtitle("PCA")
PC3v4 + geom_point(aes(color = sex), size = 4)
PC3v4 + geom_point(aes(color = biopsy), size = 4)
PC3v4 + geom_point(aes(color = person), size = 4) +
  geom_line(aes(group = person), color = "black", linetype="dotted")
PC3v4 + geom_point(aes(color = years_to_visit1), size = 3)
PC3v4 + geom_point(aes(color = DAGE), size = 4)
PC3v4 + geom_point(aes(color = center), size = 4)
PC3v4 + geom_point(aes(color = Tissue), size = 4)
PC3v4 + geom_point(aes(color = CvMlogRatio), size = 4)

# cluster
data <- as.matrix(countTable)
if(all(row.names(coldata)==colnames(data))){
  tmp <- coldata$study_id
  tmp2 <- coldata$biopsy
  tmp3 <- ifelse(tmp2 == "pre",
                 "A",
                 ifelse(tmp2 == "week12",
                        "B",
                        "C"))
  tmp <- paste0(tmp, tmp3)
  colnames(data) <- tmp
}
mat <- dist(t(data))
sampleTree <- hclust(mat)
plot(sampleTree)
coldata$batch <- NULL

# remove 02_01_week52 as outlier (as on PCA and also WGCNA)
filt <- rownames(coldata)!="02_001_week52"
coldata <- coldata[filt,]
countTable <- countTable[,filt]

# pca2
data <- as.matrix(countTable)
data <- data[order(rowVars(data), decreasing = TRUE),]
data <- data[c(1:500), ]
p <- pca(data, metadata = coldata)
screeplot(p, axisLabSize = 18, titleLabSize = 22)
biplot(p, showLoadings = TRUE,
       labSize = 5, pointSize = 5, sizeLoadingsNames = 5)
pairsplot(p,
          components = getComponents(p, seq_len(5)))
eigencorplot(p,
             metavars = c("sex", "week12.vs.pre",
                          "week52.vs.pre",
                          "week52.vs.week12", 
                          "DAGE", "CvMlogRatio", "Tissue"))
summary(row.names(coldata) == row.names(p[["rotated"]]))
pc_scores <- data.frame(p[["rotated"]], coldata[row.names(coldata) %in% row.names(p[["rotated"]]),])
PC1v2 <- ggplot(data = pc_scores, aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, color = 'grey70') +
  geom_vline(xintercept = 0, color = 'grey70') +
  xlab(paste("PC1 (", round(p[["variance"]][1], 2), "%)", sep = "")) +
  ylab(paste("PC2 (", round(p[["variance"]][2], 2), "%)", sep = "")) +
  ggtitle("PCA")
PC1v2 + geom_point(aes(color = sex), size = 4)
PC1v2 + geom_point(aes(color = biopsy), size = 4)
PC1v2 + geom_point(aes(color = person), size = 4) +
  geom_line(aes(group = person), color = "black", linetype="dotted")
PC1v2 + geom_point(aes(color = years_to_visit1), size = 3)
PC1v2 + geom_point(aes(color = DAGE), size = 4)
PC1v2 + geom_point(aes(color = center), size = 4)
PC1v2 + geom_point(aes(color = Tissue), size = 4)
PC1v2 + geom_point(aes(color = CvMlogRatio), size = 4)
PC3v4 <- ggplot(data = pc_scores, aes(x = PC3, y = PC4)) +
  geom_hline(yintercept = 0, color = 'grey70') +
  geom_vline(xintercept = 0, color = 'grey70') +
  xlab(paste("PC3 (", round(p[["variance"]][3], 2), "%)", sep = "")) +
  ylab(paste("PC4 (", round(p[["variance"]][4], 2), "%)", sep = "")) +
  ggtitle("PCA")
PC3v4 + geom_point(aes(color = sex), size = 4)
PC3v4 + geom_point(aes(color = biopsy), size = 4)
PC3v4 + geom_point(aes(color = person), size = 4) +
  geom_line(aes(group = person), color = "black", linetype="dotted")
PC3v4 + geom_point(aes(color = years_to_visit1), size = 3)
PC3v4 + geom_point(aes(color = DAGE), size = 4)
PC3v4 + geom_point(aes(color = center), size = 4)
PC3v4 + geom_point(aes(color = Tissue), size = 4)
PC3v4 + geom_point(aes(color = CvMlogRatio), size = 4)

#~~~~~~~~~~~~~~~~~ EXPORT ~~~~~~~~~~~~~~~~~
coldata <- arrange(coldata, sample)
countTable <- countTable[,match(row.names(coldata), colnames(countTable))]

if (all(row.names(coldata)==colnames(countTable))){
  save(coldata, countTable, file = "explore_data.RData")
}
