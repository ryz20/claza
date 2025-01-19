set.seed(0)
library(tidyverse)
library(DESeq2)

dds <- readRDS(file = "dds_placebo.rds")
res <- results(dds, alpha = 0.05,
               contrast=c("biopsy","week12","pre"))
res.placebo.12vspre <- na.omit(as.data.frame(res))
res <- results(dds, alpha = 0.05,
               contrast=c("biopsy","week52","pre"))
res.placebo.52vspre <- na.omit(as.data.frame(res))
res <- results(dds, alpha = 0.05,
               contrast=c("biopsy","week52","week12"))
res.placebo.52vs12 <- na.omit(as.data.frame(res))

dds <- readRDS(file = "dds_clazakizumab.rds")
res <- results(dds, alpha = 0.05,
               contrast=c("biopsy","week12","pre"))
res.clazakizumab.12vspre <- na.omit(as.data.frame(res))
res <- results(dds, alpha = 0.05,
               contrast=c("biopsy","week52","pre"))
res.clazakizumab.52vspre <- na.omit(as.data.frame(res))
res <- results(dds, alpha = 0.05,
               contrast=c("biopsy","week52","week12"))
res.clazakizumab.52vs12 <- na.omit(as.data.frame(res))


make_preranked <- function(x, name, path){
  x <- as.data.frame(na.omit(x))
  x <- rownames_to_column(x, var = "gene")
  x <- x %>% 
    mutate(GSEAmetric = 1/x$pvalue * sign(x$log2FoldChange)) %>% 
    arrange(desc(GSEAmetric)) %>% 
    dplyr::select(c("gene", "GSEAmetric"))
  if (sum(duplicated(x$gene))!=0) {
    x <- x[complete.cases(x),]
  }
  filename <- paste0(path, name, ".rnk")
  write.table(x, file = filename, quote = F, sep = "\t", row.names = F)
}

path <- "ranked/"

make_preranked(res.placebo.12vspre, 
               "res.placebo.12vspre", path)
make_preranked(res.placebo.52vspre, 
               "res.placebo.52vspre", path)
make_preranked(res.placebo.52vs12, 
               "res.placebo.52vs12", path)
make_preranked(res.clazakizumab.12vspre, 
               "res.clazakizumab.12vspre", path)
make_preranked(res.clazakizumab.52vspre, 
               "res.clazakizumab.52vspre", path)
make_preranked(res.clazakizumab.52vs12, 
               "res.clazakizumab.52vs12", path)

# run in GSEA software
# h.all.v2023.2.Hs.symbols

gsea_analysis <- function(path1, path2, geneset){
  analyses <- list.files(path1)
  dir.create(paste0(path2, "enrichmentPlots"))
  dir.create(paste0(path2, "leadingEdge"))
  dir.create(paste0(path2, "gsea_plots"))
  for (y in 1:length(analyses)){
    analysis <- analyses[y]
    tmp <- list.files(paste0(path1,analysis,"/"))
    tmp <- grep("gsea_report", tmp, value = TRUE)
    tmp <- grep("tsv", tmp, value = TRUE)
    tmp <- paste0(path1,analysis,"/",tmp)
    neg <- read.table(tmp[1], fill = TRUE, sep="\t", header=TRUE)
    pos <- read.table(tmp[2], fill = TRUE, sep="\t", header=TRUE)
    gsea_output <- rbind(neg, pos)
    pathway_names <- gsea_output$NAME
    tmp <- pathway_names
    if (substr(pathway_names[1], 1, 8) == "HALLMARK"){
      tmp <- substr(gsea_output$NAME, 10,100)
    }
    tmp <- gsub("_", " ", tmp)
    row.names(gsea_output) <- tmp
    gsea_output <- gsea_output[,c(5,6,8)]
    tmp <- strsplit(analysis,".Gsea")
    tmp <-  tmp[[1]][1]
    gsea_output$Comparision <- tmp
    gsea_output$pathway <- row.names(gsea_output)
    gsea_output <- gsea_output [gsea_output$FDR.q.val< 0.05,]
    gsea_output <- gsea_output[order(gsea_output$NES, decreasing = FALSE),]
    gsea_output$Hallmark <- factor(row.names(gsea_output), levels = unique(row.names(gsea_output)))
    gsea_output <- droplevels(gsea_output)
    gsea_output$dir <- factor(gsea_output$NES > 0, levels = c("FALSE", "TRUE"))
    if (all(gsea_output$NES > 0)) {
      pdf(file=paste0(path2,"gsea_plots/",geneset,tmp,".pdf"), useDingbats = FALSE)
      print(ggplot( data=gsea_output, aes(x=Hallmark, y=NES,  size=-log10(FDR.q.val+0.00001))) +
              geom_point(colour="red") +
              geom_hline(yintercept=0, colour="grey50") +
              labs(y="NES", x=geneset, size="-log10(FDR+0.0001)", colour="")+
              theme_bw()+
              coord_flip()+
              ggtitle(tmp))
      dev.off()
    }  else if (all(gsea_output$NES < 0)){
      pdf(file=paste0(path2,"gsea_plots/",geneset,tmp,".pdf"), useDingbats = FALSE)
      print(ggplot( data=gsea_output, aes(x=Hallmark, y=NES, size=-log10(FDR.q.val+0.00001))) +
              geom_point(colour="dodgerblue2") +
              geom_hline(yintercept=0, colour="grey50") +
              labs(y="NES", x=geneset, size="-log10(FDR+0.0001)", colour="")+
              theme_bw()+
              coord_flip()+
              ggtitle(tmp)
      )
      dev.off()
    }  else {
      pdf(file=paste0(path2,"gsea_plots/",geneset,tmp,".pdf"), useDingbats = FALSE)
      print(ggplot( data=gsea_output, aes(x=Hallmark, y=NES, colour=dir, size=-log10(FDR.q.val+0.00001))) +
              geom_point() +
              geom_hline(yintercept=0, colour="grey50") +
              labs(y="NES", x=geneset, size="-log10(FDR+0.0001)", colour="")+
              scale_colour_manual(breaks = c("FALSE","TRUE"), values=c("dodgerblue2","red"), labels=c("Downregulated\nPathway", "Upregulated\nPathway"))+
              theme_bw()+
              coord_flip()+
              ggtitle(tmp))
      dev.off()
    }
    pathway_names <- pathway_names[str_detect(pathway_names, "/", negate = T)]
    for (z in 1:length(pathway_names)){
      tmp1 <- read.table(file = paste0(path1, analysis, "/", pathway_names[z], ".tsv"), sep = '\t', header = T)
      p1 <- ggplot(data=tmp1, aes(x=RANK.IN.GENE.LIST, y=RUNNING.ES))+
        geom_line()+
        geom_hline(yintercept=0, color='black')+
        geom_vline(xintercept = 0, color='black') +
        xlab("Gene Rank")+
        ylab("Running Enrichment Score")+
        theme_gray()+
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())
      tmp2 <- tmp1[c(-1, -nrow(tmp1)),]
      col<- if (mean(tmp2$RUNNING.ES)> 0)"red" else "dodgerblue2"
      p2 <- ggplot(data=tmp2, aes(y=RANK.IN.GENE.LIST, x=1))+
        geom_violin(fill=col)+
        theme(axis.text.x = element_blank())+
        ylab("Gene Rank")+
        geom_hline(yintercept = 0)+
        theme_gray()+
        expand_limits(y=c(0,nrow(tmp1)))+
        theme(axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank())+
        coord_flip()
      pdf(file = paste0(path2, "enrichmentPlots/EP_", tmp, "_", pathway_names[z], ".pdf", sep = ""), useDingbats = FALSE)
      print(cowplot::plot_grid(p1, p2, ncol=1, rel_heights = c(4,1), align = "v"))
      dev.off()
      tmp3 <- tmp1 %>% 
        dplyr::filter(CORE.ENRICHMENT == "Yes") %>% 
        dplyr::arrange(desc(abs(RANK.METRIC.SCORE)))
            write.table(tmp3, file = paste0(path2, "leadingEdge/LE_", tmp, "_", pathway_names[z], ".tsv", sep = ""), quote=FALSE, sep='\t', col.names = T)
    }
  }
}


gsea_analysis("pathway/hallmarks/output/",
              "pathway/hallmarks/",
              "hallmarks")
gsea_analysis("pathway/kegg/output/",
              "pathway/kegg/",
              "kegg")
gsea_analysis("pathway/pbmc162k/output/",
              "pathway/pbmc162k/",
              "pbmc162k")


