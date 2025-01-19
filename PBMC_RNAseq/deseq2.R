set.seed(0)
library(tidyverse)
library(ggrepel)
library(DESeq2)
library("BiocParallel")
register(SnowParam(4))

load("explore_data.RData")

# tdy
coldata <- coldata %>% 
  select(sample, person, center, batch, sequence, sex, biopsy,
         maintenance, MI_tac, MI_CyA, MI_MPA, MI_GC, age_biopsy,
         DSA, ABMRabove0.2, ABMR_MMDx, ABMR_BANFF, GFR, GFR_percent, ABMR_MMDx)
coldata$biopsy <- factor(coldata$biopsy, levels = c("pre","week12", "week52"))
coldata$group <-  paste0(coldata$sequence, coldata$biopsy)

# placebo
coldata_placebo <- coldata[coldata$sequence=="Placebo",]
countTable_placebo <- countTable[,coldata$sequence=="Placebo"]
all(rownames(coldata)==colnames(countTable))
model <- as.formula( ~ person + batch + biopsy)
dds <- DESeqDataSetFromMatrix(countData = countTable_placebo,
                              colData = coldata_placebo,
                              design = model)
dds <- dds [rowSums(counts(dds)) >10 ,]
dds$batch <- relevel(dds$batch, ref = "1")
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds, maxit=500)
res <- results(dds, alpha = 0.05,
               contrast=c("biopsy", "week52", "week12"))
results <- na.omit(as.data.frame(res))
summary(res)
resultsNames(dds)
saveRDS(dds, file = "dds_placebo.rds")

# clz
coldata_clazakizumab <- coldata[coldata$sequence=="Clazakizumab",]
countTable_clazakizumab <- countTable[,coldata$sequence=="Clazakizumab"]
all(rownames(coldata)==colnames(countTable))
model <- as.formula( ~ person + batch + biopsy )
dds <- DESeqDataSetFromMatrix(countData = countTable_clazakizumab,
                              colData = coldata_clazakizumab,
                              design = model)
dds <- dds [rowSums(counts(dds)) >10 ,]
dds$batch <- relevel(dds$batch, ref = "1")
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds, maxit=500)
res <- results(dds, alpha = 0.05,
               contrast=c("biopsy", "week12", "pre"))
results <- na.omit(as.data.frame(res))
summary(res)
resultsNames(dds)
saveRDS(dds, file = "dds_clazakizumab.rds")

