#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# ORFik paper script, figure 2
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Analysis April/May 2020
# Contains data from Sel-TCP-seq Bohlen (HeLa cells) with CAGE from FANTOM 5
library(ORFik); library(data.table)

# Update this path (everything else should now run): ->
plotFolder <- "/export/valenfs/projects/Hakon/ORFik_paper/"

if (file.exists(file.path(plotFolder, "ORFik_paper_session.RData")))
  load(file.path(plotFolder, "ORFik_paper_session.RData"))
############################## CREATE Figure ###################################################
# Requires that Alignment script for Bohlen data is already done
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Load experiments 
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# TCP-seq (HeLa)
df <- read.experiment("sel_TCP_bohlen")
df <- df[df$rep == 1,] # We only need first replicate
df.SSU <- df[df$libtype == "SSU", ]
df.LSU <- df[df$libtype == "LSU", ]
# mRNA
df.rna <- read.experiment("bohlen_RNA")
# CAGE
df.cage <- read.experiment("cage_fantom_sub")

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Load NGS data
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# TCP-seq
outputLibs(df, type = "ofst") # Output rep 1 for LSU and SSU
SSU.5p <- convertToOneBasedRanges(SSU, addScoreColumn = TRUE, addSizeColumn = TRUE)
LSU.5p <- convertToOneBasedRanges(LSU, addScoreColumn = TRUE, addSizeColumn = TRUE)
# CAGE
cage <- filepath(df.cage, type = "ofst")
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Load annotation (Danio rerio)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# gtf
txdb <- loadTxdb(df)
# Filter and load 5' UTRs, cds and leaders to global environment
txNames100 <- filterTranscripts(txdb, 100, 100, 100)
loadRegions(txdb, names.keep = txNames100) 

# CAGE Annotation
leaders.hmm <- leaders[startSites(leaders, keep.names = FALSE, is.sorted = TRUE) > 52]
if (!file.exists(file.path(plotFolder, "TxDb_CAGE_Ovarian_cancer.db"))) {
  txdb.c <- reassignTxDbByCage(txdb, cage = fimport(cage, leaders))
  saveDb(txdb.c, file.path(plotFolder, "TxDb_CAGE_Ovarian_cancer.db"))
} else txdb.c <- loadTxdb(file.path(plotFolder, "TxDb_CAGE_Ovarian_cancer.db"))

t.cage <- filterTranscripts(txdb.c, 100, 100, 100)
loadRegions(txdb.c, c("leaders", "mrna"), t.cage, extension = ".cage")

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Count tables
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Get count tables
collapse <- FALSE
RNA_MRNA_FPKM <- countTable(df.rna, "mrna", type = "fpkm", collapse = collapse)
RNA_MRNA_FPKM <- data.table(id = rownames(RNA_MRNA_FPKM), RNA_MRNA_FPKM)
SSU_LEADERS_FPKM <- countTable(df.SSU, "leaders", type = "fpkm", collapse = collapse)
SSU_LEADERS_FPKM <- data.table(id = rownames(SSU_LEADERS_FPKM), SSU_LEADERS_FPKM)
LSU_CDS_FPKM <- countTable(df.LSU, "cds", type = "fpkm", collapse = collapse)
LSU_CDS_FPKM <- data.table(id = rownames(LSU_CDS_FPKM), LSU_CDS_FPKM)
dt <- data.table::merge.data.table(RNA_MRNA_FPKM, SSU_LEADERS_FPKM, by = "id")
dt <- data.table::merge.data.table(dt, LSU_CDS_FPKM, by = "id")
colnames(dt) <- c("txNames", "RNA_MRNA_FPKM", "SSU_LEADERS_FPKM", "LSU_CDS_FPKM")
dt[, `:=`(IR = (LSU_CDS_FPKM / SSU_LEADERS_FPKM),
          SE = (SSU_LEADERS_FPKM / RNA_MRNA_FPKM),
          TE = (LSU_CDS_FPKM / RNA_MRNA_FPKM))]
dt <- dt[RNA_MRNA_FPKM > 0 & SSU_LEADERS_FPKM > 0, ]


#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Heatmap (original & CAGE annotation)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

# Heatmaps
# SSU
cov.cage <- windowPerReadLength(leaders.cage, extendLeaders(mrna.cage, 31),
                                SSU.5p, upstream = 30, downstream = 30,
                                scoring = "transcriptNormalized", acceptedLengths = 21:75)
hm.cage <- coverageHeatMap(cov.cage, scoring = "transcriptNormalized", addFracPlot = TRUE, title = "CAGE updated annotation 40S",
                           xlab = "Position relative to Transcription start site")

cov <- windowPerReadLength(leaders.hmm, extendLeaders(mrna, extension = 51),
                           SSU.5p, upstream = 30, downstream = 30,
                           scoring = "transcriptNormalized", acceptedLengths = 21:75)
hm <- coverageHeatMap(cov, scoring = "transcriptNormalized", addFracPlot = TRUE, title = "Original annotation 40S",
                      xlab = "Position relative to Transcription start site",
                      gradient.max = max(cov.cage$score, cov$score))
# LSU
cov.cage.LSU <- windowPerReadLength(leaders.cage, extendLeaders(mrna.cage, 31),
                                    LSU.5p, upstream = 30, downstream = 30,
                                    scoring = "transcriptNormalized", acceptedLengths = 21:75)
hm.cage.LSU <- coverageHeatMap(cov.cage.LSU, scoring = "transcriptNormalized", addFracPlot = TRUE, title = "CAGE updated annotation 80S",
                               xlab = "Position relative to Transcription start site")

cov.LSU <- windowPerReadLength(leaders.hmm, extendLeaders(mrna, extension = 51),
                           LSU.5p, upstream = 30, downstream = 30,
                           scoring = "transcriptNormalized", acceptedLengths = 21:75)
hm.LSU <- coverageHeatMap(cov.LSU, scoring = "transcriptNormalized", addFracPlot = TRUE, title = "Original annotation 80S",
                      xlab = "Position relative to Transcription start site",
                      gradient.max = max(cov.cage.LSU$score, cov.LSU$score))

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Metacoverage, IR & TOP motif
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

# 1 Metacoverage
metaLSU <- ORFik:::splitIn3Tx(leaders, cds, trailers, LSU, fraction = "80S")
metaSSU <- ORFik:::splitIn3Tx(leaders, cds, trailers, SSU, fraction = "40S")

coveragePlot <- windowCoveragePlot(rbindlist(list(metaLSU, metaSSU)), colors = c('orange', 'skyblue4'))

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Kozak heatmap (IR)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
seqs <- startRegionString(cds, tx = mrna, faFile = df, upstream = 5, downstream = 4)
rates.heat <- dt[RNA_MRNA_FPKM > 10, ]
seqs <- seqs[names(seqs) %in% rates.heat$txNames]
kozakHeat <- kozakHeatmap(seqs = seqs, rate = rates.heat[txNames %in% names(seqs)]$IR,
                          start = 1, stop = max(nchar(seqs)), center = 6, type = "IR");kozakHeat
#ggslackR()
#ggsave(filename = p(plotFolder, "kozakHeatmap.png"), kozakHeat)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# TOP motif (Effect on SE FPKM)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

seqs <- startRegionString(leaders.cage[names(leaders.cage) %in% dt$txNames], NULL, df, 0, 4)
rate <- dt[txNames %in% names(seqs)]$SE
comb <- TOP.Motif.ecdf(seqs, rate, legend.position.1st = c(0.70, 0.28), legend.position.motif = c(0.70, 0.28))
#ggslackR(plot = comb)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Kozak ranking
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

rates <- dt[(txNames %in% names(cds)) & RNA_MRNA_FPKM > 10, ]
#rates <- dt[(txNames %in% names(cds)) & IR > 0.5, ]
cds_k <- cds[names(cds) %in% rates$txNames]

ranking <- ORFik:::kozak_IR_ranking(cds_k, mrna, rates, df, species = "human")
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Merge
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
library(gridExtra); library(ggplot2)
lay <- rbind(c(1,1),
             c(2,3),
             c(4, 5),
             c(6,7),
             c(8))
final <- grid.arrange(coveragePlot ,hm, hm.LSU, hm.cage, hm.cage.LSU,
                      kozakHeat, comb, ranking, layout_matrix = lay)
ggsave(paste0(plotFolder, "Figure_2_new",df@experiment,".png"), plot = final, width = 10, height = 15, dpi = 300)
save.image(file=file.path(plotFolder, "ORFik_paper_session.RData"))

library(cowplot)
lay_two <- rbind(c(1,2),
             c(3,4))
two <- grid.arrange(hm, hm.LSU, hm.cage, hm.cage.LSU, layout_matrix = lay_two)
three <- grid.arrange(kozakHeat, comb, nrow = 1)

final2 <- cowplot::plot_grid(coveragePlot, two, three, ranking, ncol = 1, labels = "AUTO", rel_heights = c(1,2,1,1), label_y = c(1, 1, 1, 1.1))
ggsave(paste0(plotFolder, "Figure_3_new",df@experiment,".pdf"), plot = final2, width = 10, height = 15, dpi = 300)
final2

# ranking2 <- kozak_IR_ranking(cds_k, mrna, rates, df, species = "human")
# coveragePlot2 <- windowCoveragePlot(rbindlist(list(metaLSU, metaSSU)), colors = c('orange', 'skyblue4'), scoring = "transcriptNormalized")
# coveragePlot2 <- windowCoveragePlot(rbindlist(list(metaLSU, metaSSU)), colors = c('orange', 'skyblue4'), scoring = "sum")
# 
# final3 <- cowplot::plot_grid(coveragePlot2, two, three, ranking, ncol = 1, labels = "AUTO", rel_heights = c(1,2,1,1), label_y = c(1, 1, 1, 1.1))
# ggsave(paste0(plotFolder, "Figure_3_new2",df@experiment,".pdf"), plot = final3, width = 10, height = 15, dpi = 300)
# final3