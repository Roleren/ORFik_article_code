#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Ribo-seq Danio rerio stages (Bazzini et al, 2014) (Merge replicates)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Article: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4193932/
library(ORFik)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Config
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
conf <- config.exper(experiment = "bazzini_2014_zebrafish_test",
                     assembly = "Danio_rerio_GRCz10",
                     type = c("Ribo-Seq", "RNA-Seq"))

# SRA Meta data download (work for ERA and DRA too)
study <- download.SRA.metadata("SRP034750", conf["fastq Ribo-Seq"])
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Download fastq files for experiment and rename (Skip if you have the files already)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Remove CO libraries, we don't want those
# Remove 5h samples
study.rfp <- study[grep("mRNA|sucrose|5h ", sample_title, invert = TRUE),]
study.rna <- study[grep("RPFs|sucrose|5h ", sample_title, invert = TRUE),]

# Download fastq files (uses SRR numbers (Run column) from study)
#download.SRA(study.rfp, conf["fastq Ribo-Seq"]) # Data already downloaded
#download.SRA(study.rna, conf["fastq RNA-Seq"])
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Merge replicates (to get 1 file per stage)
# (warning: this will make you lose possibility for DESeq: translation efficiency)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Merge Ribo-seq

infiles <- dir(conf["fastq Ribo-Seq"], "*.fastq", full.names = TRUE)
if (length(infiles) == 0) stop("Run the download&align_RFP_Bazzini14.R first!")
in_files <- c(paste0(grep(infiles, pattern = paste0("2h_",  collapse = "|"), value = TRUE), collapse = " "),
              paste0(grep(infiles, pattern = paste0("12h_", collapse = "|"), value = TRUE), collapse = " "),
              paste0(grep(infiles, pattern = paste0("24h",  collapse = "|"), value = TRUE), collapse = " "),
              paste0(grep(infiles, pattern = paste0("48h",  collapse = "|"), value = TRUE), collapse = " "))
conf["fastq merged Ribo-Seq"] <- file.path(dirname(infiles)[1], "merged")
conf["bam merged Ribo-Seq"] <- file.path(conf["bam Ribo-Seq"], "merged/")
out_files <- file.path(conf["fastq merged Ribo-Seq"],
                       c("RPF_2h.fastq", "RPF_12h.fastq", "RPF_24h.fastq", "RPF_48h.fastq"))
mergeFastq(in_files, out_files)


# Merge RNA-seq
infiles <- dir(conf["fastq RNA-Seq"], "*.fastq", full.names = TRUE)
in_files <- c(paste0(grep(infiles, pattern = paste0("2h_",  collapse = "|"), value = TRUE), collapse = " "),
              paste0(grep(infiles, pattern = paste0("12h_", collapse = "|"), value = TRUE), collapse = " "),
              paste0(grep(infiles, pattern = paste0("24h",  collapse = "|"), value = TRUE), collapse = " "),
              paste0(grep(infiles, pattern = paste0("48h",  collapse = "|"), value = TRUE), collapse = " "))
conf["fastq merged RNA-Seq"] <- file.path(dirname(infiles)[1], "merged")
conf["bam merged RNA-Seq"] <- file.path(conf["bam RNA-Seq"], "merged/")
out_files <- file.path(dirname(infiles)[1], "merged",
                       c("mRNA_2h.fastq", "mRNA_12h.fastq", "mRNA_24h.fastq", "mRNA_48h.fastq"))
mergeFastq(in_files, out_files)

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
paired.end <- study.rna$LibraryLayout == "PAIRED"

STAR.align.folder(conf["fastq merged Ribo-Seq"], conf["bam merged Ribo-Seq"], index,
                  steps = "tr-co-ge", # (trim needed: adapters found, then genome)
                  adapter.sequence = "auto", # Adapters are auto detected
                  max.cpus = 80, trim.front = 0, min.length = 20)

STAR.align.folder(conf["fastq merged RNA-Seq"], conf["bam merged RNA-Seq"], index,
                  paired.end = paired.end,
                  steps = "tr-co-ge", # (trim needed: adapters found, then genome)
                  adapter.sequence = "auto", # Adapters are auto detected
                  max.cpus = 80, trim.front = 0, min.length = 20)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Create experiment (Starting point if alignment is finished)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
library(ORFik)

create.experiment(paste0(conf["bam merged Ribo-Seq"], "/aligned/"), 
                  exper = "zf_baz14_RFP",
                  fa = annotation["genome"],
                  txdb = paste0(annotation["gtf"], ".db"),
                  organism = organism,
                  viewTemplate = FALSE)
create.experiment(paste0(conf["bam merged RNA-Seq"], "/aligned/"), 
                  exper = "zf_baz14_RNA",
                  fa = annotation["genome"],
                  txdb = paste0(annotation["gtf"], ".db"),
                  organism = organism,
                  viewTemplate = FALSE)
# Now fix experiment non-unique rows in Excel, Libre office (I keep only rep 1 & 2 here)...
df <- read.experiment("zf_baz14_RFP")
df.rna <- read.experiment("zf_baz14_RNA")
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




# counts_rfp <- countTable(df, "cds")
# counts_rfp_m <- as.matrix(counts_rfp)
# rownames(counts_rfp_m) <- rownames(counts_rfp)
# counts_rna <- countTable(df.rna, "mrna")
# counts_rna_m <- as.matrix(counts_rna[rownames(counts_rna) %in% rownames(counts_rfp),])
# rownames(counts_rna_m) <- rownames(counts_rna)
# conditions <- df$stage
# ads2 <- anota2seqDataSetFromMatrix(dataP = counts_rfp_m,
#                                   dataT = counts_rna_m,
#                                   phenoVec = conditions,
#                                   dataType = "RNAseq", transformation = "rlog",
#                                   normalize = TRUE)
# ads <- anota2seqRun(ads2, onlyGroup = TRUE, performQC = FALSE)
# anota2seqPlotFC(ads, selContrast = 1, plotToFile = FALSE)
