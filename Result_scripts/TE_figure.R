library(ORFik); library(ggplot2)
plotFolder <- "/export/valenfs/projects/Hakon/ORFik_paper/"

df.rfp <- read.experiment("zf_baz14_RFP_all")
df.rna <- read.experiment("zf_baz14_RNA_all")
# Between group differential analysis
dt <- DTEG.analysis(df.rfp, df.rna)
plot.between <- DTEG.plot(dt[variable == "Comparison: 24hpf vs 48hpf"], p.value = "")

# Within group Translational efficiency
te.within <- te.table(df.rfp[5:8,], df.rna[5:8,], collapse = TRUE)
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
ggsave(file.path(plotFolder, "DTEG_plot2.png"), final.plot,
       width = 6, height = 6, dpi = 300)
