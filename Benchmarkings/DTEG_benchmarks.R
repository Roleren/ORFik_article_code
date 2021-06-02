#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Benchmark ORFik (deltaTE) vs anota2seq vs riborex
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
library(ORFik); library(anota2seq); library(riborex)
library(data.table)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Differential translation analysis, Alexaki et al (2020) 
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Data: 
df.rfp <- read.experiment("Alexaki_Human_Ribo-Seq")
df.rna <- read.experiment("Alexaki_Human_RNA-Seq")

RFP_counts <- countTable(df.rfp, "cds", type = "summarized")
RNA_counts <- countTable(df.rna, "mrna", type = "summarized")
# Go from all isoforms to longest per gene
single_gene_variant <- filterTranscripts(df.rfp, 0, minCDS = 1, 0)
RFP_counts <- RFP_counts[rownames(RFP_counts) %in% single_gene_variant,]
RNA_counts <- RNA_counts[rownames(RNA_counts) %in% single_gene_variant,]

## ORFik
# Run time for ORFik DTEG.analysis (54 seconds)
res.orfik_a <- DTEG.analysis(df.rfp, df.rna, design = "condition", output.dir = NULL,
                           RFP_counts = RFP_counts, RNA_counts = RNA_counts)
res.orfik_a.temp <- copy(res.orfik_a)
res.orfik_a.temp[abs(te) < 1.5, Regulation := "No change"]
DTEG.plot(res.orfik_a.temp, xlim = c(-5, 5), ylim = c(-5, 5), height = 4.8, width = 5, 
          output.dir = file.path(dirname(df.rfp$filepath[1]), "QC_STATS/"), 
          relative.name = "DTEG_Alexaki_ORFik.pdf", dot.size = 1.2, p.value = "")
## Anota
# Run time for anota2seqRun (40 minutes 31 seconds)

ads <- anota2seqDataSetFromMatrix(dataP = assay(RFP_counts),
                                  dataT = assay(RNA_counts),
                                  phenoVec = df.rna$condition,
                                  dataType = "RNAseq",
                                  normalize = TRUE)

alexaki.anota.file <- file.path(dirname(df.rfp$filepath[1]), "QC_STATS", "alexaki.anota.rds")
if (file.exists(alexaki.anota.file)) {
  ads <- readRDS(alexaki.anota.file)
} else {
  custom_contrast <- matrix(c(1, -1),
                            ncol = 1)
  rownames(custom_contrast) <- levels(as.factor(df.rna$condition))
  colnames(custom_contrast) <- c("CO vs WT")
  ads <- anota2seqRun(ads, performQC = FALSE, performROT = FALSE, 
                      contrasts = custom_contrast)
  saveRDS(ads, alexaki.anota.file)
}

anota2seqPlotFC(ads, selContrast = 1, plotToFile = FALSE, contrastName = "CO vs WT")
anota2seqPlotFC(ads, selContrast = 1, plotToFile = TRUE, contrastName = "CO vs WT",
                fileStem = file.path(dirname(df.rfp$filepath[1]), "QC_STATS/DTEG_Alexaki_"))
identical(rownames(ads@mRNAAbundance@totalmRNA[[1]]), 
          res.orfik.temp$id[res.orfik.temp$Regulation != "No change"])
# They match
## Riborex
RFP.cond <- colData(RFP_counts)$condition
RNA.cond <- colData(RNA_counts)$condition
RFP <- as.data.frame(assay(RFP_counts))
RNA <- as.data.frame(assay(RNA_counts))

res.deseq <- riborex(RNA, RFP, RNA.cond, RFP.cond)
# Since riborex normalizes a bit different, we need to set lfc to 1.2, the max lfc
F9_lfc <- res.deseq[rownames(res.deseq) == rownames(ads@mRNAAbundance@totalmRNA[[1]]),]
F9_lfc$log2FoldChange
summary(res.deseq[which(res.deseq$log2FoldChange >= F9_lfc$log2FoldChange & res.deseq$padj < 0.05),])
riborex_res_alexaki <- res.deseq[which(res.deseq$log2FoldChange >= F9_lfc$log2FoldChange & res.deseq$padj < 0.05),]
riborex_dt_alexaki <- copy(res.orfik_a.temp)
riborex_hits_dt <- which(rownames(RFP) %in% rownames(res.deseq)[which(res.deseq$log2FoldChange >= F9_lfc$log2FoldChange & res.deseq$padj < 0.05)])
riborex_dt_alexaki[, Regulation := "No change" ]
riborex_dt_alexaki[riborex_hits_dt, Regulation := "Translation"]
DTEG.plot(riborex_dt_alexaki, xlim = c(-5, 5), ylim = c(-5, 5), height = 4.8, width = 5, 
          output.dir = file.path(dirname(df.rfp$filepath[1]), "QC_STATS/"), 
          relative.name = "DTEG_Alexaki_riborex.pdf",dot.size = 1.2, p.value = "")

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Differential translation analysis, Bazzini et al (2014)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Data: 
df.rfp <- read.experiment("zf_baz14_RFP_all")
df.rna <- read.experiment("zf_baz14_RNA_all")
# Go from all isoforms to longest per gene
RFP_counts <- countTable(df.rfp, "cds", type = "summarized")
RNA_counts <- countTable(df.rna, "mrna", type = "summarized")
single_gene_variant_danio <- filterTranscripts(df.rfp, 0, minCDS = 1, 0)
RFP_counts <- RFP_counts[rownames(RFP_counts) %in% single_gene_variant_danio,]
RNA_counts <- RNA_counts[rownames(RNA_counts) %in% single_gene_variant_danio,]

# Run time for ORFik DTEG.analysis (61 seconds)
res.orfik <- DTEG.analysis(df.rfp, df.rna, design = "stage", output.dir = NULL,
                           RFP_counts = RFP_counts, RNA_counts = RNA_counts)
res.orfik.temp <- copy(res.orfik)
res.orfik.temp <- res.orfik.temp[variable == "Comparison: 24hpf vs 48hpf",]
res.orfik.temp[(abs(rfp) < 1.5 & abs(rna) < 1.5), Regulation := "No change"]
res.orfik.baz <- res.orfik.temp[Regulation != "No change",]
DTEG.plot(res.orfik.temp, xlim = c(-10, 10), ylim = c(-10, 10), 
          height = 4.8, width = 5,
          output.dir = file.path(dirname(df.rfp$filepath[1]), "QC_STATS/"),
          relative.name = "DTEG_Bazzini_ORFik.pdf", dot.size = 0.6, p.value = "")


# Run time for anota2seqRun (51 minutes 16 seconds)

ads2_data <- anota2seqDataSetFromMatrix(dataP = assay(RFP_counts),
                                  dataT = assay(RNA_counts),
                                  phenoVec = df.rna$stage,
                                  dataType = "RNAseq",
                                  normalize = TRUE, 
                                  filterZeroGenes = TRUE)
# Custom contrast (3 of them, 24hpf vs all three),
# only 24hpf vs 48hpf is important for analysis
custom_contrast <- matrix(c(c(0, 1, 0, -1), c(-1, 1, 0, 0), c(0, 1, -1, 0)),
                          ncol = 3)
rownames(custom_contrast) <- levels(as.factor(df.rna$stage))
colnames(custom_contrast) <- c("24hpf vs 48hpf", "24hpf vs 12hpf", "24hpf vs 2hpf")

bazzini.anota.file <- file.path(dirname(df.rfp$filepath[1]), "QC_STATS", "bazzini.anota.rds")
if (file.exists(bazzini.anota.file)) {
  ads2 <- readRDS(bazzini.anota.file)
} else {
  ads2 <- anota2seqAnalyze(ads2_data, contrasts = custom_contrast)
  saveRDS(ads2, bazzini.anota.file)
}

ads2_out <- anota2seqSelSigGenes(ads2, 
                                 selDeltaP = 1.5, selDeltaT = 1.5, 
                                 selDeltaPT = 1.5, selDeltaTP = 1.5,
                                 maxP = 0.1) # Dummy p-value for plotting
ads2_out <- anota2seqRegModes(ads2_out)
#summary(ads2_out@buffering@apvStats[[1]])
anota2seqPlotFC(ads2_out, selContrast = 1, plotToFile = FALSE, contrastName = "24hpf vs 48hpf")
anota2seqPlotFC(ads2_out, selContrast = 1, plotToFile = TRUE, contrastName = "24hpf vs 48hpf",
                fileStem = file.path(dirname(df.rfp$filepath[1]), "QC_STATS/DTEG_Bazzini_"))
a <- anota2seqGetOutput(ads2_out, output="singleDf", selContrast=1)
stopifnot(identical(nrow(a), nrow(res.orfik.temp)))
# Assign same filter as ORFik (alpha value 0.1 -> padj 0.1)
a <- a[a$translatedmRNA.apvRvmPAdj < 0.1 |
         a$totalmRNA.apvRvmPAdj < 0.1  |
         a$translation.apvRvmPAdj < 0.1 |
         abs(a$buffering.apvEff) > 1.5,] # No pvalue on buffering, since deltaTE uses LFC

## Riborex
RFP.cond <- colData(RFP_counts)$stage
RNA.cond <- colData(RNA_counts)$stage
RFP <- as.data.frame(assay(RFP_counts))
RNA <- as.data.frame(assay(RNA_counts))

res.deseq <- riborex(RNA, RFP, RNA.cond, RFP.cond)
# Since riborex normalizes a bit different, we need to set lfc to 1.2, the max lfc
summary(res.deseq)

riborex_res_bazzini <- res.deseq[which( abs(res.deseq$log2FoldChange) > 1.5  & res.deseq$padj < 0.05),]
riborex_dt_bazzini <- copy(res.orfik.temp)
riborex_hits_dt <- which(rownames(RFP) %in% rownames(res.deseq)[which(abs(res.deseq$log2FoldChange) > 1.5  & res.deseq$padj < 0.05)])
riborex_txnames_baz <- rownames(RFP)[riborex_hits_dt]
riborex_dt_bazzini[, Regulation := "No change" ]
riborex_dt_bazzini[riborex_hits_dt, Regulation := "Translation"]
riborex_dt_bazzini[(abs(rfp) < 1.5 & abs(rna) < 1.5), Regulation := "No change"]
DTEG.plot(riborex_dt_bazzini, xlim = c(-10, 10), ylim = c(-10, 10), height = 4.8, width = 5, 
          output.dir = file.path(dirname(df.rfp$filepath[1]), "QC_STATS/"), 
          relative.name = "DTEG_Bazzini_riborex.pdf", dot.size = 0.6, p.value = "")
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Compare results
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

{
  # Create a data.table with statistics
  stats <- data.table(Type = c("ORFik", "anota2seq", "Riborex"))
  #a_significant <- a[a$translation.apvRvmPAdj < 0.1 | a$buffering.apvRvmPAdj < 0.1 | a$totalmRNA.apvRvmPAdj < 0.1,]
  
  # Get anota values
  txNames_anota <- rownames(a)
  translated <- txNames_anota[a$singleRegMode == "translation"]
  mrna_abundance <- txNames_anota[a$singleRegMode == "abundance"]
  buffering <- txNames_anota[a$singleRegMode == "buffering"]
  all_active <- txNames_anota[a$singleRegMode != "background"]
  length(translated); length(mrna_abundance); length(buffering)
  length(all_active); nrow(res.orfik.baz)
  # Get Riborex value:
  txNames_riborex <- riborex_txnames_baz
  
  # Genes pre-excluded:
  # ORFik_filtered <- res.orfik.temp$id[is.na(res.orfik.temp$rna.padj) | is.na(res.orfik.temp$rfp.padj) | is.na(res.orfik.temp$te.padj)]
  # anota_filtered <- 0
  # stats[, pre_excluded := c(length(ORFik_filtered), length(anota_filtered))]
  # 
  # ORFik_filtered_over <- round(sum(ORFik_filtered %in% anota_filtered) / length(ORFik_filtered) * 100, 2)
  # anota_filtered_over <- round(sum(anota_filtered %in% ORFik_filtered) / length(anota_filtered) * 100, 2)
  # stats[, pre_excluded_overlap := c(ORFik_filtered_over, anota_filtered_over)]
  # DTEG genes
  stats[, DTEGs := c(nrow(res.orfik.baz), length(all_active), length(txNames_riborex))]
  
  ORFik_DTEG_over <- round(sum(res.orfik.baz$id %in% all_active) / nrow(res.orfik.baz) * 100, 2)
  anota_DTEG_over <- round(sum(all_active %in% res.orfik.baz$id) / length(all_active) * 100, 2)
  riborex_DTEG_over <- round(sum(txNames_riborex %in% res.orfik.baz$id) / length(txNames_riborex) * 100, 2)
  stats[, DTEGs_overlap := c(ORFik_DTEG_over, anota_DTEG_over, riborex_DTEG_over)]
  round(sum(txNames_riborex %in% all_active) / length(txNames_riborex) * 100, 2)
  # DTEG genes (Translating)
  ORFik.translated <- res.orfik.baz$Regulation == "Translation"
  stats[, Translation := c(sum(ORFik.translated), length(translated), NA)]
  
  ORFik.sub <- res.orfik.baz[ORFik.translated,]
  ORFik_DTEG_over <- round(sum(ORFik.sub$id %in% translated) / nrow(ORFik.sub) * 100, 2)
  anota_DTEG_over <- round(sum(translated %in% ORFik.sub$id) / length(translated) * 100, 2)
  stats[, Translation_overlap := c(ORFik_DTEG_over, anota_DTEG_over, NA)]
  # DTEG genes (mRNA abundance)
  ORFik.abundance <- res.orfik.baz$Regulation == "mRNA abundance"
  stats[, mRNA_abundance := c(sum(ORFik.abundance), length(mrna_abundance), NA)]
  
  ORFik.sub <- res.orfik.baz[ORFik.abundance,]
  ORFik_DTEG_over <- round(sum(ORFik.sub$id %in% mrna_abundance) / nrow(ORFik.sub) * 100, 2)
  anota_DTEG_over <- round(sum(mrna_abundance %in% ORFik.sub$id) / length(mrna_abundance) * 100, 2)
  stats[, mRNA_abundance_overlap := c(ORFik_DTEG_over, anota_DTEG_over, NA)]
  # DTEG genes (buffering)
  ORFik.buffering <- res.orfik.baz$Regulation == "Buffering"
  stats[, Buffering := c(sum(ORFik.buffering), length(buffering), NA)]
  
  ORFik.sub <- res.orfik.baz[ORFik.buffering,]
  ORFik_DTEG_over <- round(sum(ORFik.sub$id %in% buffering) / nrow(ORFik.sub) * 100, 2)
  anota_DTEG_over <- round(sum(buffering %in% ORFik.sub$id) / length(buffering) * 100, 2)
  stats[, Buffering_overlap := c(ORFik_DTEG_over, anota_DTEG_over, NA)]
  
  # Add row with differences:
  dif_row <- t(data.table(c("Difference", abs(stats[1,-1] - stats[2, -1]))))
  stats <- rbindlist(list(stats, dif_row))
}
stats

fwrite(stats, file = file.path(dirname(df.rfp$filepath[1]), "QC_STATS", "benchmark_stats_DTEG_all.csv"))
