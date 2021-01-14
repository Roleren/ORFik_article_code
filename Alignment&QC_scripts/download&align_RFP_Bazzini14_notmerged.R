#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Ribo-seq Danio rerio stages (Bazzini et al, 2014)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Article: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4193932/
library(ORFik)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Settings (This is the only Custom part per user, rest you can just run)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# This is where you want your annotation and STAR index
annotation <- "/export/valenfs/data/references/Zv10_zebrafish/"
# Where to download fastq files
fastq.dir.rfp <- "/export/valenfs/data/raw_data/Ribo-Seq/bazzini_2014_zebrafish/"
fastq.dir.rna <- "/export/valenfs/data/raw_data/RNA-Seq/bazzini_zebrafish_2014/"
# Where you want mapped bam files
bam.dir.rfp <- "/export/valenfs/data/processed_data/Ribo-seq/bazzini_2014_zebrafish/"
bam.dir.rna <- "/export/valenfs/data/processed_data/RNA-seq/bazzini_2014_zebrafish/"

# SRA Meta data download (work for ERA and DRA too)
study <- download.SRA.metadata("SRP034750", fastq.dir.rfp)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Download fastq files for experiment and rename (Skip if you have the files already)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Remove CO libraries, we don't want those
# Remove 5h samples
study.rfp <- study[grep("mRNA|sucrose|5h ", sample_title, invert = TRUE),]
study.rna <- study[grep("RPFs|sucrose|5h ", sample_title, invert = TRUE),]

# Download fastq files (uses SRR numbers (Run column) from study)
download.SRA(study.rfp, fastq.dir.rfp)
download.SRA(study.rna, fastq.dir.rna)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Annotation
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Which organism is this, scientific name, like "Homo sapiens" or "Danio rerio"
organism <- study$ScientificName[1] # Usually you find organism here, else set it yourself


annotation <- getGenomeAndAnnotation(
  organism = organism,
  genome = TRUE, GTF = TRUE,
  phix = TRUE, ncRNA = TRUE, tRNA = TRUE, rRNA = TRUE,
  output.dir = annotation,
  assembly_type = "primary_assembly"
)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# STAR index
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
index <- STAR.index(annotation, wait = TRUE)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Alignment (with depletion of phix, rRNA, ncRNA and tRNAs) & (with MultiQC of final STAR alignment)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# All data is single end
paired.end <- study$LibraryLayout == "PAIRED"
paired.end.all <- all(paired.end)

alignment <-
  STAR.align.folder(fastq.dir.rfp, bam.dir.rfp, index,
                    paired.end = paired.end.all,
                    steps = "tr-co-ge", # (trim needed: adapters found, then genome)
                    adapter.sequence = "auto", # Adapters are auto detected
                    max.cpus = 80, trim.front = 0, min.length = 20)

alignment.rna <-
  STAR.align.folder(fastq.dir.rna, bam.dir.rna, index,
                    paired.end = paired.end.all,
                    steps = "tr-co-ge", # (trim needed: adapters found, then genome)
                    adapter.sequence = "auto", # Adapters are auto detected
                    max.cpus = 80, trim.front = 0, min.length = 20)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Create experiment (Starting point if alignment is finished)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
library(ORFik)

create.experiment(paste0(alignment, "/aligned/"), exper = "zf_baz14_RFP_all",
                  fa = annotation["genome"],
                  txdb = paste0(annotation["gtf"], ".db"),
                  organism = organism,
                  viewTemplate = FALSE)
create.experiment(paste0(alignment.rna, "/aligned/"), exper = "zf_baz14_RNA_all",
                  fa = annotation["genome"],
                  txdb = paste0(annotation["gtf"], ".db"),
                  organism = organism,
                  viewTemplate = FALSE)
# Now fix experiment non-unique rows in Excel, Libre office (I keep only rep 1 & 2 here)...
df <- read.experiment("zf_baz14_RFP_all")
df.rna <- read.experiment("zf_baz14_RNA_all")
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Convert files and run Annotation vs alignment QC
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
convertLibs(df, reassign.when.saving = TRUE)
ORFik:::countTable_regions(df, geneOrTxNames = "tx", longestPerGene = FALSE)

convertLibs(df.rna, reassign.when.saving = TRUE)
ORFik:::countTable_regions(df.rna, geneOrTxNames = "tx", longestPerGene = FALSE)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# P-shifting of Ribo-seq reads:
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Pre-pshifting heatmaps
heatMapRegion(df, region = c("TIS", "TTS"), shifting = "5prime", type = "ofst",
              outdir = paste0(dirname(df$filepath[1]), "/QC_STATS/heatmaps/pre-pshift/"))

shiftFootprintsByExperiment(df)
# Pshifted heatmaps
heatMapRegion(df, region = c("TIS", "TTS"), shifting = "5prime", type = "pshifted",
              outdir = paste0(dirname(df$filepath[1]), "/QC_STATS/heatmaps/pshifted/"))



counts_rfp <- countTable(df, "cds")
counts_rfp_m <- as.matrix(counts_rfp)
rownames(counts_rfp_m) <- rownames(counts_rfp)
counts_rna <- countTable(df.rna, "mrna")
counts_rna_m <- as.matrix(counts_rna[rownames(counts_rna) %in% rownames(counts_rfp),])
rownames(counts_rna_m) <- rownames(counts_rna)
conditions <- df$stage
ads2 <- anota2seqDataSetFromMatrix(dataP = counts_rfp_m,
                                  dataT = counts_rna_m,
                                  phenoVec = conditions,
                                  dataType = "RNAseq", transformation = "rlog",
                                  normalize = TRUE)
ads <- anota2seqRun(ads2, onlyGroup = TRUE, performQC = FALSE)
anota2seqPlotFC(ads, selContrast = 1, plotToFile = FALSE)
library(limma)

# DESEQ
df.rfp <- read.experiment("zf_baz14_RFP_all")
df.rna <- read.experiment("zf_baz14_RNA_all")

library(DESeq2); library(SummarizedExperiment); collapse <- FALSE
RFP_DESEQ <- countTable(df.rfp, "cds", type = "summarized", collapse = collapse)
colData(RFP_DESEQ)$stage <- bamVarName(df.rfp, T, T, F, T, T, T)
colData(RFP_DESEQ)$libtype <- bamVarName(df.rfp, T, T, T, T, T, F)
colData(RFP_DESEQ)$condition <- bamVarName(df.rfp, T, F, T, T, T, T)

RNA_DESEQ <- countTable(df.rna, "mrna", type = "summarized", collapse = collapse)
colData(RNA_DESEQ)$stage <- bamVarName(df.rna, T, T, F, T, T, T)
colData(RNA_DESEQ)$libtype <- bamVarName(df.rna, T, T, T, T, T, F)
colData(RNA_DESEQ)$condition <- bamVarName(df.rna, T, F, T, T, T, T)


se <- cbind(assay(RFP_DESEQ), assay(RNA_DESEQ))
colData <- rbind(colData(RFP_DESEQ), colData(RNA_DESEQ))
combined_se <- SummarizedExperiment(se,
                                    rowRanges = rowRanges(RFP_DESEQ),
                                    colData = colData)
combined_se <- combined_se[rowMeans(assay(combined_se)) > 1,]

combined_DESEQ <- DESeqDataSet(combined_se, design = ~ libtype + stage + libtype:stage)
dds <- DESeq2::DESeq(combined_DESEQ)
res <- DESeq2::results(dds, contrast = c("stage", "24hpf", "48hpf"))
res.final <- res[!is.na(res$padj),]
print(summary(res))
a <- rownames(res.final[(res.final$padj < 0.05) & (abs(res.final$log2FoldChange) > 1),])

plotMA(res)



