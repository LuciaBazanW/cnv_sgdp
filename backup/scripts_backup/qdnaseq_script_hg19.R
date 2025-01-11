library(QDNAseq)

library(Biobase)

library(BSgenome.Hsapiens.UCSC.hg38)
library("future")

bins <- getBinAnnotations(binSize=15)

readCounts <- binReadCounts(bins, bamfiles="/branchinecta/jbazanwilliamson/SGDP/mayan/SAMEA3302866/SAMEA3302866.sort.bam")

pdf(file="read_counts_SAME2866_hg_19_15k.pdf")

plot(readCounts, logTransform=FALSE, ylim=c(-5,200))
dev.off()

readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE)

readCountsFiltered <- estimateCorrection(readCountsFiltered)

readCountsFiltered <- applyFilters(readCountsFiltered, chromosomes=NA)

pdf(file="noise_plot_hg19_15k.pdf")

noisePlot(readCountsFiltered)

dev.off()

copyNumbers <- correctBins(readCountsFiltered)

copyNumbersNormalized <- normalizeBins(copyNumbers)

copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)

pdf(file="cnvs_same2866_hg19_15k.pdf")

plot(copyNumbersSmooth)

dev.off()

exportBins(copyNumbersSmooth, file="SAME866s_hg19_15k.txt")
exportBins(copyNumbersSmooth, file="SAME866s_hg19_15k.igv", format="igv")

copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")

copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)

pdf(file="cnvs_aftersegmentation_same2866_hg19_15k.pdf")
plot(copyNumbersSegmented)
dev.off()

copyNumbersCalled <- callBins(copyNumbersSegmented)

pdf(file="cnvs_called_hg19_15k.pdf")
plot(copyNumbersCalled)
dev.off()

#exportBins(copyNumbersCalled, file="SAME886cnvs.vcf", format="vcf")
#exportBins(copyNumbersCalled, file="SAME886cnvs.seg", format="seg")

exportBins(copyNumbersCalled, file="SAME886cnvs_hg19_15k.vcf", format="vcf")
exportBins(copyNumbersCalled, file="SAME886cnvs_hg19_15k.seg", format="seg")


