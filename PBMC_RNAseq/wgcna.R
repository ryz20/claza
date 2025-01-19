set.seed(0)
library(WGCNA)
library(cluster)
library(flashClust)
library(tidyverse)
library(DESeq2)
library(limma)
options(stringsAsFactors = F)

load("explore_data.RData")

coldata <- coldata %>% 
  select(sample, person, batch, sequence, sex, biopsy,
         age_biopsy, DSA, ABMRabove0.2, ABMR_MMDx, ABMR_BANFF, GFR, IgG, GFR_percent, IgG_percent)
coldata$biopsy <- factor(coldata$biopsy, levels = c("pre","week12", "week52"))

coldata$group <-  paste0(coldata$sequence, coldata$biopsy)
coldata$batch <- as.factor(coldata$batch)
coldata$person <- as.factor(coldata$person)

tmp1= binarizeCategoricalVariable(coldata$sex, includePairwise = TRUE, includeLevelVsAll = FALSE)
tmp2= binarizeCategoricalVariable(coldata$group, includePairwise = TRUE, includeLevelVsAll = FALSE)
tmp3= binarizeCategoricalVariable(coldata$biopsy, includePairwise = TRUE, includeLevelVsAll = FALSE)

coldata <- data.frame(coldata, tmp1, tmp2, tmp3)
rm(tmp1, tmp2, tmp3)

# batch correction
model <- as.formula( ~ person + batch)
dds <- DESeqDataSetFromMatrix(countData = countTable,
                              colData = coldata,
                              design = model)
dds <- dds [rowSums(counts(dds)) > 20,]
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
mat <- assay(vsd)
mm <- model.matrix(~ person, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$batch, design=mm)
mat <- t(mat)
coldata <- coldata[match(row.names(mat), row.names(coldata)),]
table(rownames(coldata) == rownames(mat))
rm(vsd, mm, model, countTable)

A.sample <- adjacency(t(mat), type = "distance")
k.sample <- as.numeric(apply(A.sample, 2, sum)) - 1
Z.k.sample = scale(k.sample)
threshold_outlier <- -5
outlierC <- ifelse(Z.k.sample < threshold_outlier, "red", "black")
sampleTree <- flashClust(as.dist(1 - A.sample), method = "average")
coldata <- coldata %>% 
  select(c("DSA", 
           "ABMR_MMDx", 
           "ABMR_BANFF",
           "age_biopsy",
           "male.vs.female", 
           "Clazakizumabweek52.vs.Clazakizumabweek12",
           "Clazakizumabweek12.vs.Clazakizumabpre",
           "Placeboweek52.vs.Placeboweek12",
           "Placeboweek12.vs.Placebopre",
           "GFR",
           "IgG_percent"))
traitColors <- data.frame(numbers2colors(coldata, signed = FALSE))
dimnames(traitColors)[[2]] <- paste(names(coldata), "C", sep = "")
datColors <- data.frame(outlierC, traitColors)
pdf(file="wgcna/plots/sample_tree.pdf", width = 14, height = 7)
plotDendroAndColors(sampleTree, 
                    groupLabels = names(datColors),
                    colors = datColors,
                    main = "Sample dendrogram and trait heatmap") 
dev.off()
saveRDS(mat, file = "wgcna/objects/mat.rds")
rm(sampleTree, A.sample, datColors, outlierC, traitColors, Z.k.sample, k.sample, threshold_outlier)

# Choose power threshold
powers <- c(1:20) 
sft <- pickSoftThreshold(mat, powerVector = powers, networkType = "signed")
par(mfrow=c(1,2))
plot(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n",main=paste("Scale independence"))
text(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,col="red")
abline(h=0.90,col="red")
plot(sft$fitIndices[,1],sft$fitIndices[,5],type="n",
     xlab="Soft Threshold (power)",ylab="Mean Connectivity",main=paste("Mean connectivity"))
text(sft$fitIndices[,1],sft$fitIndices[,5],labels=powers,col="red") 
power = sft$fitIndices %>% 
  filter(SFT.R.sq >0.9) %>% 
  pull(Power) %>% 
  min()

A <-  adjacency(mat, power = power, type="signed", corFnc = "bicor") 
dissTOM <-  TOMdist(A, TOMType = "signed") 
geneTree <-  flashClust(as.dist(dissTOM), method="average") 
moduleLabels1 <- cutreeDynamic(dendro = geneTree,
                               distM = dissTOM,
                               method = "hybrid",
                               deepSplit = 1,
                               pamRespectsDendro = F,
                               minClusterSize = 30
)
moduleColors1 <- labels2colors(moduleLabels1) 
MEList1 <- moduleEigengenes(mat, colors = moduleColors1) 
MEs1 <-  MEList1$eigengenes
plotEigengeneNetworks(MEs1,
                      "",
                      marDendro=c(0,4,1,2),
                      marHeatmap=c(3,4,1,2),
                      cex.lab=0.8,
                      xLabelsAngle=90) 

mergingThresh <- dynamicMergeCut(nrow(coldata), mergeCor = 0.9, Zquantile = 2.35)
mergingThresh <- 0.3
merged <- mergeCloseModules(mat, moduleColors1, cutHeight = mergingThresh)
moduleColors2 <-  merged$colors
MEs2 <-  merged$newMEs
plotEigengeneNetworks(MEs2,
                      "",
                      marDendro=c(0,4,1,2),
                      marHeatmap=c(3,4,1,2),
                      cex.lab=0.8,
                      xLabelsAngle=90) 
variables_interest <- colnames(coldata)
varColors <- data.frame(matrix(vector(),ncol(mat),0),stringsAsFactors = FALSE)
for(i in 1:(length(variables_interest))){
  var.name <- as.character(variables_interest[i])
  var <- as.data.frame(select(coldata, var.name))
  names(var) <-  as.character(variables_interest[i])
  tmp <- as.numeric(cor(mat, var, use="p"))
  tmp <- numbers2colors(tmp,signed = T)
  varColors <- cbind(varColors, tmp)
}
datColors <- data.frame(moduleColors1, moduleColors2, varColors)
plotDendroAndColors(geneTree,
                    colors = datColors,
                    groupLabels=c("modules1", "modulesmerged", colnames(coldata)),
                    dendroLabels=FALSE,
                    hang=0.03,
                    addGuide=TRUE,
                    guideHang=0.05) 

mColors <- moduleColors2
nGenes <-  ncol(mat)
nSamples <-  nrow(mat)
MEs <-  moduleEigengenes(mat,mColors)$eigengenes
MEsOrdered = orderMEs(MEs)
modTraitCor = cor(MEsOrdered, coldata, use = "p")
modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(",
                   signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = modTraitCor, 
               xLabels = names(coldata),
               yLabels = names(MEsOrdered), 
               ySymbols = names(MEsOrdered),
               colorLabels =FALSE,
               colors=greenWhiteRed(50),
               textMatrix=textMatrix,
               setStdMargins = FALSE, 
               cex.text = 0.8, 
               zlim = c(-1,1),
               main = paste("Module-trait relationships")
) 

# MM
datKME <- signedKME(mat, MEs2)
colnames(modTraitCor) <- c("DSA",
                           "ABMR MMDx",
                           "ABMR BANFF",
                           "Donor Age",
                           "M vs. F",
                           "Claza B",
                           "Claza A",
                           "Placebo B",
                           "Placebo A",
                           "eGFR",
                           "Serum IgG %")

modTraitCor <- modTraitCor[1:(nrow(modTraitCor)-1),]
modTraitP <- modTraitP[1:(nrow(modTraitP)-1),]
textMatrix <- textMatrix[1:(nrow(textMatrix)-1),]
coldata <- coldata[1:(nrow(coldata)-1),]
modTraitCor <- modTraitCor[,c(1:10, 13)]
modTraitP <- modTraitP[,c(1:10, 13)]
textMatrix <- textMatrix[,c(1:10, 13)]

tmp <- row.names(modTraitCor)
tmp <- sub('..', '', tmp)
tmp <- str_to_title(tmp)
row.names(modTraitCor) <- tmp
row.names(modTraitP) <- tmp

labeledHeatmap(Matrix = modTraitCor, 
               xLabels = rep("", ncol(modTraitCor)),
               yLabels = rep("", nrow(modTraitCor)), 
               ySymbols = rownames(modTraitCor),
               colorLabels =FALSE,
               colors=blueWhiteRed(50),
               setStdMargins = FALSE, 
               cex.text = 0.6, 
               zlim = c(-1,1),
               plotLegend = T)