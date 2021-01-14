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

# DTEG.plot <- function(dt, output.dir = NULL,
#                       p.value = 0.05,
#                       plot.title = "", width = 6,
#                       height = 6, dot.size = 0.4,
#                       xlim = c(-5, 5), ylim = c(-10, 10)) {
#   color.values <- c("black", "orange4", "purple", "darkgreen")
#   p.caption <- if (p.value != "") {
#     labs(caption = paste("P-value <", p.value))
#   } else NULL
#   p.title <- if (plot.title != "") {
#     ggtitle(label = plot.title)
#   } else NULL
#   dt[, Regulation := Status]
#   dt$Status <- NULL
#   
#   dot.size <- rep(dot.size, nrow(dt))
#   dot.size[dt$Regulation != "No change"] <- dot.size[1]*2
#   plot.between <- ggplot(data = dt,
#                          aes(x = rna, y = rfp, color = Regulation)) +
#     geom_point(alpha = 0.6, size = dot.size) +
#     scale_color_manual(values = color.values) +
#     theme_minimal() +
#     geom_hline(aes(yintercept =  0), alpha = 0.2, color = "red") +
#     geom_vline(aes(xintercept =  0), alpha = 0.2, color = "red") +
#     xlab("mRNA (log2 fold change)") +
#     ylab("RFP (log2 fold change)") +
#     p.title +
#     p.caption +
#     facet_wrap(~ variable, ncol = 2) +
#     xlim(xlim) + ylim(ylim) +
#     guides(color = guide_legend(override.aes = list(alpha = 0.8, size = 1.3)))
#   plot(plot.between)
#   if (!is.null(output.dir)) {
#     ggsave(file.path(output.dir, "DTEG_plot.png"), plot.between,
#            width = width, height = height, dpi = 300)
#   }
#   return(plot.between)
# }
# 
# 
# te_rna.plot <- function(dt, output.dir = NULL,
#                         filter.rfp = 1, filter.rna = 1,
#                         plot.title = "",
#                         width = 6, height = "auto",
#                         dot.size = 0.4) {
#   
#   if (height == "auto") height <- 3+length(unique(dt$variable))
#   caption <- paste("Filter: RFP >", filter.rfp, " & mRNA >", filter.rna, "(FPKM)")
#   if (nrow(df.rfp) > 1 & nrow(df.rna) == 1)
#     caption <- paste(subtitle, "(Single mRNA sample)")
#   p.caption <- if (filter.rfp != "") {
#     labs(caption = caption)
#   } else NULL
#   p.title <- if (plot.title != "") {
#     ggtitle(label = plot.title)
#   } else NULL
#   
#   plot <- ggplot(data = dt, aes(x = rna_log10, y = LFC_TE)) +
#     geom_point(alpha = 0.2, size = dot.size) +
#     theme_minimal() +
#     geom_hline(aes(yintercept =  0), alpha = 0.2, color = "red") +
#     xlab("mRNA FPKM (log10)") +
#     ylab("TE (log2)") +
#     p.caption +
#     p.title +
#     ggtitle(label = plot.title) +
#     xlim(c(filter.rna, filter.rna + 2.5)) +
#     facet_wrap(~ variable, nrow = 1)
#   
#   plot(plot)
#   if (!is.null(output.dir)) {
#     ggsave(file.path(output.dir, "TE_within.png"), plot,
#            width = width, height = height, dpi = 300)
#   }
#   return(plot)
# }


