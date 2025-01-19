set.seed(0)
library(tidyverse)
library(ggrepel)
library(DESeq2)

load("tidied_data.RData")

countTable <- countTable[rowSums(countTable) >58,]
tmp <- countTable>0
tmp2 <- rowSums(tmp)
tmp2 <- tmp2>29
countTable <- countTable[tmp2,]

# QC
boxplot <- as.data.frame(countTable) %>% 
  pivot_longer(names(countTable), names_to = "sample", values_to = "count")
ggplot(boxplot, aes(x = sample, y = count)) +
  geom_boxplot()+
  scale_y_log10()
ggplot(boxplot, aes(x = reorder(sample, count, FUN = sum), y = count)) +
  geom_boxplot()+
  scale_y_log10()
ggplot(boxplot, aes(x = reorder(sample, count, FUN = median), y = count)) +
  geom_boxplot()+
  scale_y_log10()

# remove the 1x low quality sample
tmp <- sort(rowSums(t(countTable)))
tmp <- names(tmp)[1]
filt <- !(colnames(countTable) %in% tmp)
countTable <- countTable[,filt]

boxplot <- countTable %>% 
  pivot_longer(names(countTable), names_to = "sample", values_to = "count")
ggplot(boxplot, aes(x = sample, y = count)) +
  geom_boxplot()+
  scale_y_log10()
ggplot(boxplot, aes(x = reorder(sample, count, FUN = sum), y = count)) +
  geom_boxplot()+
  scale_y_log10()
ggplot(boxplot, aes(x = reorder(sample, count, FUN = median), y = count)) +
  geom_boxplot()+
  scale_y_log10()

coldata <- coldata[filt,]

# sex
tmp <- as.data.frame(t(countTable[row.names(countTable) %in% c("XIST", "UTY"),]))
tmp$logratio <- log2((tmp$XIST+1)/(tmp$UTY+1))
tmp$sex <- as.factor(as.character(lapply(tmp$logratio, function(x) {
  if (x < -2)"male" 
  else if (x > 5)"female" 
  else "UNKNOWN"})))
hist(tmp$logratio)
ggplot(tmp, aes(x=row.names(tmp), y=(logratio), colour=sex))+
  geom_point()+
  geom_hline(yintercept = 0)
all(tmp$sex == coldata$sex)

# PCA
rm(boxplot, tmp, filt)
library(PCAtools)
coldata <- cbind(coldata,
                 WGCNA::binarizeCategoricalVariable(coldata$biopsy,
                                  includePairwise = TRUE,
                                  includeLevelVsAll = FALSE))
data <- as.matrix(countTable)
data <- DESeq2::vst(data)
data <- data[order(rowVars(data), decreasing = TRUE),]
data <- data[c(1:1000), ]
p <- pca(data, metadata = coldata)
screeplot(p, axisLabSize = 18, titleLabSize = 22)
biplot(p, showLoadings = TRUE,
       labSize = 5, pointSize = 5, sizeLoadingsNames = 5)
pairsplot(p,
          components = getComponents(p, seq_len(6)))
eigencorplot(p,
             metavars = c("sex", "week12.vs.pre",
                          "week52.vs.pre",
                          "week52.vs.week12", 
                          "batch", "age_biopsy"))
summary(row.names(coldata) == row.names(p[["rotated"]]))
pc_scores <- data.frame(p[["rotated"]], coldata[row.names(coldata) %in% row.names(p[["rotated"]]),])
PC1v2 <- ggplot(data = pc_scores, aes(x = PC1, y = PC2, label = sample)) +
  geom_hline(yintercept = 0, color = 'grey70') +
  geom_vline(xintercept = 0, color = 'grey70') +
  xlab(paste("PC1 (", round(p[["variance"]][1], 2), "%)", sep = "")) +
  ylab(paste("PC2 (", round(p[["variance"]][2], 2), "%)", sep = "")) +
  ggtitle("PCA")+
  geom_text_repel()
library(RColorBrewer)
PC1v2 + geom_point(aes(color = sex), size = 4)
PC1v2 + geom_point(aes(color = biopsy), size = 4)
PC1v2 + geom_point(aes(color = person), size = 4) +
  geom_line(aes(group = person), color = "black", linetype="dotted")
PC1v2 + geom_point(aes(color = batch), size = 4)
PC1v2 + geom_point(aes(color = years_to_visit1), size = 3)
PC1v2 + geom_point(aes(color = age_biopsy), size = 4)
PC1v2 + geom_point(aes(color = center), size = 4)
PC3v4 <- ggplot(data = pc_scores, aes(x = PC3, y = PC4, label=sample)) +
  geom_hline(yintercept = 0, color = 'grey70') +
  geom_vline(xintercept = 0, color = 'grey70') +
  xlab(paste("PC3 (", round(p[["variance"]][3], 2), "%)", sep = "")) +
  ylab(paste("PC4 (", round(p[["variance"]][4], 2), "%)", sep = "")) +
  ggtitle("PCA")+
  geom_text_repel()
PC3v4 + geom_point(aes(color = sex), size = 4)
PC3v4 + geom_point(aes(color = biopsy), size = 4)
PC3v4 + geom_point(aes(color = person), size = 4) +
  geom_line(aes(group = person), color = "black", linetype="dotted")
PC3v4 + geom_point(aes(color = batch), size = 4)
PC3v4 + geom_point(aes(color = years_to_visit1), size = 3)
PC3v4 + geom_point(aes(color = age_biopsy), size = 4)
PC3v4 + geom_point(aes(color = center), size = 4)

# remove 2nd outlier
coldata <- coldata[row.names(coldata)!="01_002_week52",]
countTable <- countTable[,colnames(countTable)!="01_002_week52"]

# PCA2
data <- as.matrix(countTable)
data <- DESeq2::vst(data)
data <- data[order(rowVars(data), decreasing = TRUE),]
data <- data[c(1:1000), ]
p <- pca(data, metadata = coldata)
screeplot(p, axisLabSize = 18, titleLabSize = 22)
biplot(p, showLoadings = TRUE,
       labSize = 5, pointSize = 5, sizeLoadingsNames = 5)
pairsplot(p,
          components = getComponents(p, seq_len(6)))
eigencorplot(p,
             metavars = c("sex", "week12.vs.pre",
                          "week52.vs.pre",
                          "week52.vs.week12", 
                          "batch", "age_biopsy"))
summary(row.names(coldata) == row.names(p[["rotated"]]))
pc_scores <- data.frame(p[["rotated"]], coldata[row.names(coldata) %in% row.names(p[["rotated"]]),])
PC1v2 <- ggplot(data = pc_scores, aes(x = PC1, y = PC2, label = sample)) +
  geom_hline(yintercept = 0, color = 'grey70') +
  geom_vline(xintercept = 0, color = 'grey70') +
  xlab(paste("PC1 (", round(p[["variance"]][1], 2), "%)", sep = "")) +
  ylab(paste("PC2 (", round(p[["variance"]][2], 2), "%)", sep = "")) +
  ggtitle("PCA")+
  geom_text_repel()
PC1v2 + geom_point(aes(color = sex), size = 4)
PC1v2 + geom_point(aes(color = biopsy), size = 4)
PC1v2 + geom_point(aes(color = person), size = 4) +
  geom_line(aes(group = person), color = "black", linetype="dotted")
PC1v2 + geom_point(aes(color = batch), size = 4)
PC1v2 + geom_point(aes(color = years_to_visit1), size = 3)
PC1v2 + geom_point(aes(color = age_biopsy), size = 4)
PC1v2 + geom_point(aes(color = center), size = 4)
PC3v4 <- ggplot(data = pc_scores, aes(x = PC3, y = PC4, label=sample)) +
  geom_hline(yintercept = 0, color = 'grey70') +
  geom_vline(xintercept = 0, color = 'grey70') +
  xlab(paste("PC3 (", round(p[["variance"]][3], 2), "%)", sep = "")) +
  ylab(paste("PC4 (", round(p[["variance"]][4], 2), "%)", sep = "")) +
  ggtitle("PCA")+
  geom_text_repel()
PC3v4 + geom_point(aes(color = sex), size = 4)
PC3v4 + geom_point(aes(color = biopsy), size = 4)
PC3v4 + geom_point(aes(color = person), size = 4) +
  geom_line(aes(group = person), color = "black", linetype="dotted")
PC3v4 + geom_point(aes(color = batch), size = 4)
PC3v4 + geom_point(aes(color = years_to_visit1), size = 3)
PC3v4 + geom_point(aes(color = age_biopsy), size = 4)
PC3v4 + geom_point(aes(color = center), size = 4)
pc1_genes <- as.data.frame(p[["loadings"]])
pc1_genes <- pc1_genes %>% 
  select(PC1) %>% 
  arrange(PC1)
pc1_genes <- c(row.names(pc1_genes)[1:40])
plotDat <- as.matrix(countTable)
plotDat <- plotDat[pc1_genes,]
plotDat <- DESeq2::varianceStabilizingTransformation(plotDat)
anno_info <- as.data.frame(coldata[,c("sequence", "biopsy","age_biopsy", "person")])
anno_info <- anno_info %>% 
  arrange(age_biopsy)
plotDat <- plotDat[,match(row.names(anno_info), colnames(plotDat))]
summary(colnames(plotDat)==row.names(anno_info))
pheatmap::pheatmap(mat = plotDat,
                   scale="row", 
                   cluster_rows = T, 
                   cluster_cols = F, 
                   clustering_distance_cols="euclidean",
                   show_rownames = TRUE,
                   show_colnames = FALSE,
                   annotation_col = anno_info)
pheatmap::pheatmap(mat = plotDat,
                   scale="row", 
                   cluster_rows = T, 
                   cluster_cols = T, 
                   clustering_distance_cols="euclidean",
                   show_rownames = TRUE,
                   show_colnames = FALSE,
                   annotation_col = anno_info)

# clustering
data <- as.matrix(countTable)
data <- DESeq2::rlog(data)
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

#~~~~~~~~~~~~~~~~~ EXPORT ~~~~~~~~~~~~~~~~~
coldata <- arrange(coldata, sample)
countTable <- countTable[,match(row.names(coldata), colnames(countTable))]

if (all(row.names(coldata)==colnames(countTable))){
  save(coldata, countTable, file = "explore_data.RData")
}

