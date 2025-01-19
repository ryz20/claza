library(tidyverse)
library(WGCNA)
set.seed(0)

# read in adjacency matrices and module colours
datAdjKidney <- readRDS("A_kidney.rds")
datAdjBlood <- readRDS("A_PBMC.rds")
moduleColorsKidney <- readRDS("newmodulecolors.rds")

setLabels <- c("kidney", "blood")
multiExpr <- list(kidney = list(data = datAdjKidney),
                  blood = list(data = datAdjBlood))
moduleColorsKidney <- moduleColorsKidney
multiColor <- list(kidney = moduleColorsKidney)

nPermutations <- 200
set.seed(0)

system.time({
mp <- modulePreservation(multiExpr,
                         multiColor,
                         dataIsExpr = FALSE,
                         referenceNetworks = 1,
                         nPermutations = nPermutations,
                         randomSeed = 1,
                         quickCor = 0,
                         verbose = 3)
})

ref=1; test=2

obs.PreservationStats=mp$preservation$observed[[ref]][[test]]
Z.PreservationStats=mp$preservation$Z[[ref]][[test]]

obs.PreservationStats
Z.PreservationStats

modColors <- rownames(obs.PreservationStats)
moduleSize <- obs.PreservationStats$moduleSize
selectModules <- !(modColors %in% c("grey", "gold"))
point.label <- modColors[selectModules]
medianRank <- obs.PreservationStats$medianRank.pres
Zsummary <- Z.PreservationStats$Zsummary.pres

par(mfrow=c(1,2), mar=c(4.5,4.5,2.5,1))
plot(moduleSize[selectModules], medianRank[selectModules], col=1,
     bg=modColors[selectModules], pch=21, main = "medianRank Preservation",
     cex=2, ylab="medianRank", xlab="Module size", log="x")
labelPoints(moduleSize[selectModules],medianRank[selectModules],point.label,cex=1,offs=0.03)
plot(moduleSize[selectModules], Zsummary[selectModules], col=1,
     bg=modColors[selectModules], pch=21, main = "Zsummary Preservation",
     cex=2, ylab="Zsummary", xlab="Module size", log="x")
labelPoints(moduleSize[selectModules],Zsummary[selectModules],point.label,cex=1,offs=0.03)
abline(h=0); abline(h=2, col="blue", lty=2); abline(h=10, col="red", lty=2)

datExprKidney <- readRDS("C:/Users/royyn/Medicine/academic/projects/clazakizumab_2024/analysis_kidney_microarray_tissueadjusted/wgcna/objects/mat.rds")
MEskidney <- readRDS("wgcna/objects/mergedMEs.rds")
for (i in 1:ncol(MEskidney)){
  whichmodule <- sub('..', '', colnames(MEskidney)[i])
  Eigengene=MEskidney[,i]
  datExprModule <- datExprKidney[,moduleColorsKidney==whichmodule]
  par(mfrow=c(1,1), mar = c(0.3, 5.5, 3, 2))
  plotMat(t(scale(datExprModule)), cex.axis=2, nrgcols=30, rlabels=F,
          rcols=whichmodule, main = paste("heatmap", whichmodule, "module"))
}