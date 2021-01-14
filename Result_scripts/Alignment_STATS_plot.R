#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Custom alignment plot Bohlen
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
library(ORFik); library(ggplot2); library(data.table)
df <- read.experiment("sel_TCP_bohlen")

# Update this path (everything else should now run): ->
plotFolder <- "/export/valenfs/projects/Hakon/ORFik_paper/"

dt_f <- fread(file.path(dirname(df$filepath[1]), "..", "full_process.csv"))
dt_f <- dt_f[c(1:2,4:5),]
dt_f$sample_id <- seq(4)
col_names <- colnames(dt_f)[-c(1,2)]
col_names[14] <- "trimmed reads"
col_names <- gsub("_", " ", col_names)
col_names <- gsub("of ", "", col_names)
col_names <- gsub("#", "", col_names)
col_names <- gsub("total mapped reads %", "% mapped reads", col_names)
col_names <- gsub("^ ", "", col_names)
col_names <- gsub("Uniquely mapped reads %", "% uniquely mapped reads", col_names)
col_names <- gsub("Uniquely mapped reads ", "uniquely mapped reads", col_names)
col_names <- gsub("total mapped reads ", "mapped reads", col_names)
col_names <- gsub("mapped reads", "mapped", col_names)
col_names <- gsub("reads multimapped", "multimapped", col_names)
col_names <- gsub("contamination", "contaminants", col_names)
colnames(dt_f)[-c(1,2)] <- col_names
temp_colnames <- col_names
col_names <- col_names[c(13,15,7, 1, seq(17)[-c(13,15,7, 1)])]



dt_plot <- melt(dt_f, id.vars = c("sample", "sample_id"))
dt_plot[, variable := factor(variable, 
                              levels = col_names, ordered = TRUE)]

gg_STAR <- ggplot(dt_plot,
                  aes(x=sample_id, y = value, fill = sample)) +
  geom_bar(aes(), stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  ylab("Value (log10)") +
  xlab("Samples") +
  facet_wrap(  ~ variable, scales = "free", ncol = 3) +
  scale_y_log10() +
  theme_minimal() + 
  theme(legend.position = "bottom", strip.text = element_text(size=8)) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  ggtitle("Alignment & trimming statistics") + 
  scale_y_continuous(n.breaks = 3)
  
gg_STAR
ggsave(file.path(plotFolder, "alignment_plot.png"), gg_STAR,
       width = 7, height = 9, dpi = 300)

gg_STAR_f <- gg_STAR + coord_flip()
ggsave(file.path(plotFolder, "alignment_plot_flipped.png"), gg_STAR_f,
       width = 18, height = 9)
