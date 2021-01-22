#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Meta coverage sup. figure
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

library(ORFik)
plotFolder <- "/export/valenfs/projects/Hakon/ORFik_paper/"
df <- read.experiment("sel_TCP_bohlen")
loadRegions(df, names.keep = filterTranscripts(df, 100, 100, 100))

varNames <- bamVarName(df)
outputLibs(df, leaders, type = "ofst")
windowSize <- 100

coverage <- bplapply(varNames, function(x, leaders, cds, trailers,
                                        windowSize) {
  message(x)
  ORFik:::splitIn3Tx(leaders, cds, trailers,
             get(x), fraction = x,
             windowSize = windowSize)
}, leaders = leaders, cds = cds, trailers = trailers,
windowSize = windowSize)
coverage <- rbindlist(coverage)

coverage[, fraction := gsub("LSU", "80S", fraction)]
coverage[, fraction := gsub("SSU", "43S", fraction)]
coverage[, feature := gsub("leaders", "5´UTRs", feature)]
coverage[, feature := gsub("cds", "CDS", feature)]
coverage[, feature := gsub("trailers", "3´UTRs", feature)]

scores = c("sum", "zscore", "transcriptNormalized")
format <- ".png"
col <- c(rep('skyblue4', 3), rep('orange', 3))
title <- "Coverage metaplot";pasteDir <- ORFik:::pasteDir
a <- windowCoveragePlot(coverage, scoring = scores[1], title = title, colors = col) + scale_x_continuous(breaks = c(50, 100)) + theme(panel.spacing.x = unit(0.2, "lines"))
aa <- windowCoveragePlot(coverage, scoring = scores[2], title = title, colors = col) + scale_x_continuous(breaks = c(50, 100)) + theme(panel.spacing.x = unit(0.2, "lines"))
aaa <- windowCoveragePlot(coverage, scoring = scores[3], title = title, colors = col) + scale_x_continuous(breaks = c(50, 100)) + theme(panel.spacing.x = unit(0.2, "lines"))

final <- cowplot::plot_grid(a, aa, aaa, nrow = 1, labels = "AUTO")
ggsave(pasteDir(plotFolder, paste0(df@experiment,"_cp_all", format)),final,
       height = 10, width = 10, dpi = 300)
ggsave(pasteDir(plotFolder, paste0(df@experiment,"_cp_all", ".pdf")),final,
       height = 10, width = 10, dpi = 300)
