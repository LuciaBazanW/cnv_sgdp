library(CNVRanger)
library(tidyr)
library(Gviz)
library(ggplot2)

setwd("/Users/luciabazan/Documents/GitHub/cnv_sgdp")

# CNVs calls

calls <- read.csv("call_cnvr_cnvnator/individual_cnv_cnvr.txt", sep = "\t")
calls <-  drop_na(calls)
names(calls)[names(calls) == "CNV_Value"] <- "state"
#calls <- calls[calls$Sample_ID %in% c('LP6005441-DNA_H03', 'LP6005441-DNA_G07', 
#'LP6005441-DNA_B08','LP6005441-DNA_D10'),]
calls <- calls[ -c(4,6:7,9) ]
calls
grl <- GenomicRanges::makeGRangesListFromDataFrame(calls, 
                                                   split.field="Sample_ID", 
                                                   keep.extra.columns=TRUE)
grl


# SORT 
grl <- GenomicRanges::sort(grl)
grl

ra <- RaggedExperiment::RaggedExperiment(grl)
ra


## Identifying recurrent regions 
cnvrs <- populationRanges(grl, density=0.1, est.recur=TRUE)
cnvrs

subset(cnvrs, pvalue < 0.05)

plotRecurrentRegions(cnvrs, genome="hg19", chr="chr1")


## OVERLAP with functional genomic regions 
library(AnnotationHub)
library(BiocManager)
install("BiocFileCache")
ah <- AnnotationHub::AnnotationHub()
ahDb <- AnnotationHub::query(ah, pattern = c("Homo sapiens", "Ensembl"))
ahDb
ahEdb <- ahDb[["AH10684"]]
ahEdb 
hs.genes <- ensembldb::genes(ahEdb)
GenomeInfoDb::seqlevelsStyle(hs.genes) <- "UCSC"


BiocHubsShiny::BiocHubsShiny()


library(rtracklayer)

g <- readGFF("/Users/luciabazan/Downloads/gencode.v37lift37.annotation.gtf")
names(g)[names(g) == "seqid"] <- "chr"
g$chr<-gsub("chr","",as.character(g$chr))
g <- g[!duplicated(g$gene_id), ]
rownames(g) <- g$gene_id

granges_genes <- makeGRangesFromDataFrame(g, keep.extra.columns=TRUE)

granges_genes


## Finding and illustrating overlaps
olaps <- GenomicRanges::findOverlaps(granges_genes, cnvrs, ignore.strand=TRUE)
qh <- S4Vectors::queryHits(olaps)
sh <- S4Vectors::subjectHits(olaps)
cgenes <- granges_genes[qh]
cgenes$type <- cnvrs$type[sh]
subset(cgenes, select = "type")


cnvOncoPrint(grl, cgenes)


## CNV and RNA-seq 
rcounts  <- read.csv("data/rnaseq_counts_tpm.csv", row.names=1)
names(rcounts)[names(rcounts) == "LP6005441.DNA_D10"] <- "LP6005441-DNA_D10"
names(rcounts)[names(rcounts) == "LP6005441.DNA_B08"] <- "LP6005441-DNA_B08"
names(rcounts)[names(rcounts) == "LP6005441.DNA_G07"] <- "LP6005441-DNA_G07"
names(rcounts)[names(rcounts) == "LP6005441.DNA_H03"] <- "LP6005441-DNA_H03"
#rcounts <- rcounts[ -c(1:2) ]
rcounts <- as.matrix(rcounts)

rcounts

library(SummarizedExperiment)
rranges <- GenomicRanges::granges(granges_genes)[rownames(rcounts)]
rse <- SummarizedExperiment(assays=list(rcounts=rcounts), rowRanges=rranges)
rse



result_exp_filter <- cnvEQTL(cnvrs, grl, rse, verbose = TRUE, min.samples = 2, min.state.freq = 2, 
               de.method='edgeR') ##, filter.by.expr=FALSE

result_dt_exp_filter <-as.data.frame(result)




pdf(file= "degs_expression_filter.pdf" ) 
for (x in 1:75) {
  (r <- GRanges(names(result_exp_filter)[x]))
  print(x)
  try(
    plotEQTL(cnvr=r, genes=result_exp_filter[[x]], genome="hg19", cn="CN1"), 
    silent=TRUE)
  
}
dev.off()

result <- cnvEQTL(cnvrs, grl, rse, verbose = TRUE, min.samples = 2, min.state.freq = 2, 
                  de.method='edgeR', filter.by.expr=FALSE) ##, 

result_dt <-as.data.frame(result)


pdf(file= "degs_expression_filter_no.pdf" ) 
for (x in 1:91) {
  (r <- GRanges(names(result)[x]))
  print(x)
  try(
    plotEQTL(cnvr=r, genes=result[[x]], genome="hg19", cn="CN1"), 
    silent=TRUE)
  
}
dev.off()
