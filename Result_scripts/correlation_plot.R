#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Correlation plot (refined)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
library(ORFik); library(ggplot2); library(GGally)
plotFolder <- "/export/valenfs/projects/Hakon/ORFik_paper/"
pasteDir <- ORFik:::pasteDir
df.rfp <- read.experiment("zf_baz14_RFP")
df.rna <- read.experiment("zf_baz14_RNA")



data_for_pairs <- countTable(df.rfp, "mrna", type = "fpkm")
colnames(data_for_pairs) <- paste0("RFP_", c("2hpf", "12hpf", "24hpf", "48hpf"))

message("  - log2 scaled fpkm")
point_settings <- list(continuous = wrap("points", alpha = 0.3, size=0.1), 
                       combo = wrap("dot", alpha = 0.4, size=0.2))
paired_plot <- ggpairs(as.data.frame(log2(data_for_pairs + 1)),
                       columns = 1:ncol(data_for_pairs), 
                       lower = point_settings)

ggsave(pasteDir(plotFolder, "cor_plot_log2.png"), paired_plot,
       height = 200, width = 200, units = 'mm', dpi = 300)
