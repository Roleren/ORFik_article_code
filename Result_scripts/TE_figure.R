library(ORFik); library(ggplot2)
plotFolder <- "/export/valenfs/projects/Hakon/ORFik_paper/"

df.rfp <- read.experiment("zf_baz14_RFP_all")
df.rna <- read.experiment("zf_baz14_RNA_all")
# Between group differential analysis
dt <- DTEG.analysis(df.rfp, df.rna, output.dir = NULL)
# old: "Comparison: 24hpf vs 48hpf"
dt.subset <- dt[variable == "Comparison: 12hpf vs 24hpf",]
dt.subset[, variable := "Comparison: Sample 2 vs Sample 3"]
plot.between <- DTEG.plot(dt.subset, p.value = "", output.dir = NULL, xlim = c(-10, 10))

# Within group Translational efficiency
#df.rfp[5:8,], df.rna[5:8,]
te.within <- te.table(df.rfp[c(3:6),], df.rna[c(3:6),], collapse = TRUE)
te.within[variable == "12hpf", variable := "Sample 2"]
te.within[variable == "24hpf", variable := "Sample 3"]
plot.within <- te_rna.plot(te.within, filter.rfp = "")

# Put plots together
final.plot <- cowplot::plot_grid(plot.between, plot.within, ncol = 1,
                                 labels = "AUTO", rel_heights = c(1.5, 1))
final.plot
ggsave(file.path(plotFolder, "DTEG_plot.png"), final.plot,
       width = 6, height = 6, dpi = 300)

final.plot <- cowplot::plot_grid(plot.within, plot.between, ncol = 1,
                                 labels = "AUTO", rel_heights = c(1, 1.5))
final.plot
ggsave(file.path(plotFolder, "DTEG_plot_12_24.png"), final.plot,
       width = 6, height = 6, dpi = 300)
