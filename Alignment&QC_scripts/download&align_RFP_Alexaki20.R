#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Ribo-seq HEK293 (2020) Investigative analysis of quality of new Ribo-seq data
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Article: https://f1000research.com/articles/9-174/v2
# Design: Wild type (WT) vs codon optimized (CO) (gene F9)
library(ORFik); library(data.table)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Config
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Specify paths wanted for NGS data, genome, annotation and STAR index
# If you use local files, make a conf variable with existing directories
conf <- config.exper(experiment = "Alexaki_Human",
                     assembly = "Homo_sapiens_GRCh38_110",
                     type = c("Ribo-Seq", "RNA-Seq"))

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Download fastq files for experiment and rename (Skip if you have the files already)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# SRA Meta data download (work for ERA and DRA too)
study <- download.SRA.metadata("PRJNA591214", outdir = conf["fastq Ribo-Seq"])
# Subset
study.rfp <- study[grep("mRNA-Seq", sample_title, invert = TRUE),]
study.rna <- study[grep("Ribo-Seq", sample_title, invert = TRUE),]
# Download fastq files (uses SRR numbers (RUN column) from study))
download.SRA(study.rfp, conf["fastq Ribo-Seq"],
             rename = study.rfp$sample_title)
download.SRA(study.rna, conf["fastq RNA-Seq"],
             rename = study.rna$sample_title)

# Which organism is this, scientific name, like "Homo sapiens" or "Danio rerio"
organism <- study$ScientificName[1] # Usually you find organism here, else set it yourself
paired.end.rfp <- study.rfp$LibraryLayout == "PAIRED"
paired.end.rna <- study.rna$LibraryLayout == "PAIRED"
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Annotation
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

annotation <- getGenomeAndAnnotation(
  organism = organism,
  genome = TRUE, GTF = TRUE,
  phix = TRUE, ncRNA = TRUE, tRNA = TRUE, rRNA = TRUE,
  output.dir = conf["ref"],
  assembly_type = "primary_assembly"
)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# STAR index
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Remove max.ram = 20 and SAsparse = 2, if you have more than 64GB ram
index <- STAR.index(annotation, wait = TRUE, max.ram = 20, SAsparse = 2)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Alignment (with depletion of phix, rRNA, ncRNA and tRNAs) & (with MultiQC of final STAR alignment)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

STAR.align.folder(conf["fastq Ribo-Seq"], conf["bam Ribo-Seq"], index,
                  paired.end = paired.end.rfp,
                  steps = "tr-co-ge", # (trim needed: adapters found, then genome)
                  adapter.sequence = "auto", # Adapters are auto detected
                  trim.front = 0, min.length = 20)

STAR.align.folder(conf["fastq RNA-Seq"], conf["bam RNA-Seq"], index,
                  paired.end = paired.end.rna,
                  steps = "tr-co-ge", # (trim needed: adapters found, then genome)
                  adapter.sequence = "auto", # Adapters are auto detected
                  trim.front = 0, min.length = 20)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Collapse data to send to ribotoolkit online tool
#   - original analysis done march 12th 2021
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
path <- file.path(conf["bam Ribo-Seq"], "trim") # trim folder
files <- dir(path, "trimmed", full.names = TRUE)
collapse.fastq(files, compress = TRUE)
# Settings for ribotoolkit: 
# Now upload these files, using collapsed fasta
# Homo sapiens, GRcH38, p-site read lengths: 20 -- 32 , 
# Allowed mis-matches in mapping: 3
# Allowed maximum multiple mapping: 10
# Rest of arguments are default: 
# Job IDs for finished jobs: url + id goes to page:
# http://rnabioinfor.tch.harvard.edu/RiboToolkit/batch_status.php?jobid=
# WT, replicates: gsTyI9huWQijr5vu, LrQ184ijgV2cw6qa, 3VaF8NplsslsKdjr
# CO, replicates: hgYLsqL508SVaSKr, Spa3mJWgYnDBXKD8, VPOKoMtXNFYdlbuY
# Ribotoolkit Group-run for DTEG analysis found further down
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Create experiment (Starting point if alignment is finished)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
library(ORFik)
create.experiment(paste0(conf["bam Ribo-Seq"], "/aligned/"),
                  exper = conf["exp Ribo-Seq"],
                  fa = annotation["genome"],
                  txdb = paste0(annotation["gtf"], ".db"),
                  organism = organism,
                  pairedEndBam = paired.end.rfp,
                  rep = c(1,2,3,1,2,3),
                  condition = c(rep("CO", 3), rep("WT", 3)),
                  viewTemplate = FALSE)
create.experiment(paste0(conf["bam RNA-Seq"], "/aligned/"),
                  exper = conf["exp RNA-Seq"],
                  fa = annotation["genome"],
                  txdb = paste0(annotation["gtf"], ".db"),
                  organism = organism,
                  pairedEndBam = paired.end.rna,
                  rep = c(1,2,3,1,2,3),
                  condition = c(rep("CO", 3), rep("WT", 3)),
                  viewTemplate = FALSE)

library(ORFik)
df.rfp <- read.experiment("Alexaki_Human_Ribo-Seq")
df.rna <- read.experiment("Alexaki_Human_RNA-Seq")
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Convert files and run Annotation vs alignment QC
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# General QC
# ORFik is fast because it uses a lot of memory
# 32 GB memory computers need to run QC single threaded for large genomes
BPPARAM <- BiocParallel::SerialParam()
ORFikQC(df.rfp)
ORFikQC(df.rna)

# Save countTable for Ribo-toolkit benchmark for DTEG:
# Upload this as the gene-counts
counts <- countTable(df.rna, region = "mrna")
counts <- data.table(Geneid = txNamesToGeneNames(rownames(counts), df.rna),
                     counts)
counts <- counts[!duplicated(Geneid),]
# Next line you download the valid genes 
valids <- fread("http://rnabioinfor.tch.harvard.edu/RiboToolkit/data/testdata/geneID/hg38.geneID.txt")
counts <- counts[counts$Geneid %in% valids$GeneID, ]
rna.ribotoolkit.path <- file.path(dirname(df.rna$filepath[1]), "QC_STATS","rna.counts.ribotoolkit.csv")
fwrite(counts, file = rna.ribotoolkit.path, row.names = FALSE, sep = " ")
# Now use this gene expression table in the group run for ribotoolkit
# Job id was: ID: zJkDBMhCpLOkJ7VB
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# P-shifting of Ribo-seq reads:
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# From ORFikQC it looks like 20, 21, 27 and 28 are candidates for Ribosomal footprints
shiftFootprintsByExperiment(df.rfp, accepted.lengths = c(20:21, 27:28))
# Ribo-seq specific QC
remove.experiments(df.rfp) # Remove loaded data (it is not pshifted)
ORFik:::RiboQC.plot(df.rfp)
remove.experiments(df.rfp)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Create heatmaps (Ribo-seq)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Pre-pshifting
heatMapRegion(df.rfp, region = c("TIS", "TTS"), shifting = "5prime", type = "ofst",
              outdir = paste0(dirname(df.rfp$filepath[1]), "/QC_STATS/heatmaps/pre-pshift/"))
heatMapRegion(df.rfp, region = c("TIS"), shifting = "5prime", type = "ofst", 
              acceptedLengths = c(21,27,28),
              outdir = paste0(dirname(df.rfp$filepath[1]), "/QC_STATS/heatmaps/pre-pshift_subset/"))
txdb <- loadTxdb(df.rfp)
txNames <- filterTranscripts(txdb, 51, 70, 1)
center <- loadRegion(txdb, "cds")[txNames]
mrna <- loadRegion(txdb, "mrna")[txNames]
ORFik:::heatMapL(center, mrna, df.rfp, 
                 outdir = paste0(dirname(df.rfp$filepath[1]), "/QC_STATS/heatmaps/pre-pshift_subset/"),
         scores = c("transcriptNormalized"), 
         upstream = c(50, 0), downstream = c(29, 0),
         addFracPlot = TRUE, location = "TIS",
         shifting = "5prime", skip.last = FALSE,
         acceptedLengths = c(21, 27,28), type = "ofst")

remove.experiments(df.rfp)
# After pshifting
heatMapRegion(df.rfp, region = c("TIS", "TTS"), shifting = "5prime", type = "pshifted",
              outdir = paste0(dirname(df.rfp$filepath[1]), "/QC_STATS/heatmaps/pshifted/"))
