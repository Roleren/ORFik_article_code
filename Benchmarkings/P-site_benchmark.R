#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# P-site detection, benchmark ORFik vs ribotoolkit (which uses plastid)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Data: Alexaki et al 2020
# Steps are: 
# 1. Load shifts from ORFik and Ribotoolkit
# 2. Also pshift libraries to ribotoolkit shifts
# 3. Load both data and run the ORFik::orfScore function
# 4. plot the results as boxplots

library(ORFik); library(data.table); library(ggplot2)
df.rfp <- read.experiment("Alexaki_Human_Ribo-Seq")
shifts_orfik <- ORFik::shifts.load(df.rfp)
names(shifts_orfik) <- c(paste0("CO", 1:3), paste0("WT", 1:3))

shifts_orfik <- cbind(data.table::rbindlist(shifts_orfik),
                      library = c(rep("CO1", 4),
                                  rep("CO2", 4),
                                  rep("CO3", 4),
                                  rep("WT1", 4),
                                  rep("WT2", 4),
                                  rep("WT3", 4)))
colnames(shifts_orfik) <- c("length", "offset", "library")
shifts_orfik <- shifts_orfik[length %in% c(21, 27, 28),]

# ribotoolkit p-sites extracted from results explained in 
# alignment script of Alexaki data
shifts_ribotoolkit <- data.frame(length = c(21, 27, 28), 
                                 offset = c(c(3, 11, 12),
                                            c(3, 12, 12),
                                            c(NA, 11, 12),
                                            c(9, 11, 12), 
                                            c(3, 12, 12),
                                            c(3, 12, 12)),
                                 library = c(rep("CO1", 3),
                                             rep("CO2", 3),
                                             rep("CO3", 3),
                                             rep("WT1", 3),
                                             rep("WT2", 3),
                                             rep("WT3", 3)))
# Sanity checks 
all(shifts_orfik$length == shifts_ribotoolkit$length)
all(shifts_orfik$library == shifts_ribotoolkit$library)
# Now make difference
offset_difference <- data.frame(library = shifts_orfik$library,
                                length = shifts_orfik$length, 
                                offset_orfik = shifts_orfik$offset,
                                offset_ribotoolkit = -shifts_ribotoolkit$offset,
                                offset_diff = shifts_orfik$offset + shifts_ribotoolkit$offset)
offset_difference$in_frame <- offset_difference$offset_diff %% 3 == 0
fwrite(offset_difference, file.path(dirname(df.rfp$filepath[1]), "QC_STATS", "benchmark_pshifting.csv"))
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Make pshifted data for ribotoolkit offsets
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
outputLibs(df.rfp, type = "ofst") # non pshifted files
libs <- bamVarName(df.rfp)
save.dir <- file.path(dirname(filepath(df.rfp, 
                                       type = "default")[1]),
                      "pshifted_ribotoolkit")
if (!dir.exists(save.dir)) { # Do only once
  dir.create(save.dir)
  j = 1
  for(i in unique(shifts_ribotoolkit$library)) {
    shifts <- shifts_ribotoolkit[shifts_ribotoolkit$library %in% i,]
    colnames(shifts) <- c("fraction", "offsets_start", "lib")
    shifts <- shifts[!is.na(shifts$offsets_start),]
    shifts$offsets_start <- -shifts$offsets_start
    if(nrow(shifts) == 0) continue
    gr <- shiftFootprints(footprints = get(libs[j]), shifts)
    export.ofst(gr, file.path(save.dir, paste0(libs[j], ".ofst")))
    j = j + 1
  }
  remove.experiments(df.rfp)
}


#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Check ORFscores for the two ways to shift
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Filter out cds with fpkm < 1
counts <- countTable(df.rfp, "cds", type = "fpkm")
names <- rownames(counts)
names <- names[rowMins(as.matrix(counts)) >= 1] 
cds <- loadRegion(df.rfp, "cds", names.keep = names)

# orfScore for ORFik pshifts
outputLibs(df.rfp, type = "pshifted")
libs <- bamVarName(df.rfp)
res <- bplapply(libs, FUN = function(lib, cds) { 
  return(orfScore(cds, get(lib)[readWidths(get(lib)) %in% c(21, 27, 28)], is.sorted = TRUE)$ORFScores)
}, cds = cds)
remove.experiments(df.rfp) # unload data
# Add some names and convert
names(res) <- libs
data.table::setDT(res) # Will give 1 column per library
res # Now by columns
summary(res)

# orfScore for ribotoolkit pshifts
ofst.tool <- file.path(save.dir, paste0(libs, ".ofst"))
all(file.exists(ofst.tool))
res_tool <- bplapply(ofst.tool, FUN = function(lib, cds) { 
  return(orfScore(cds, fimport(lib), is.sorted = TRUE)$ORFScores)
}, cds = cds)
# Add some names and convert
names(res_tool) <- libs
data.table::setDT(res_tool) # Will give 1 column per library
res_tool # Now by columns
summary(res_tool)

orfik.melt <- melt(res)
orfik.melt$source <- "ORFik"
tool.melt <- melt(res_tool)
tool.melt$source <- "RiboToolkit"
plot.data <- rbindlist(list(orfik.melt, tool.melt))

plot <- ggplot(data = plot.data) + 
  geom_boxplot(aes(x = variable, y = value, fill = source)) + 
  ylab("ORFscore") + 
  xlab("Library")
plot
ggsave(file.path(save.dir, "ORFscore_comparison.png"), plot, dpi = 300)
ggsave(file.path(save.dir, "ORFscore_comparison.pdf"), plot, dpi = 300)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#  ORFscore without TIS and TTS region
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

# Remove 15 bases from 5' end and 3' end of CDS
cds_sub <- cds[widthPerGroup(cds, F) > 100]
start_transform <- startRegion(cds_sub,
                               upstream = -15, downstream = 15)
stop_transform <- stopRegion(cds_sub,
                               upstream = 15, downstream = -15)
startSites(start_transform)[1]
stopSites(stop_transform)[1]
cds_sub_done <- extendLeaders(cds_sub, start_transform)
cds_sub_done <- extendTrailers(cds_sub_done, stop_transform)

# orfScore for ORFik pshifts
outputLibs(df.rfp, type = "pshifted")
libs <- bamVarName(df.rfp)
res <- bplapply(libs, FUN = function(lib, cds) { 
  return(orfScore(cds, get(lib)[readWidths(get(lib)) %in% c(21, 27, 28)], is.sorted = TRUE)$ORFScores)
}, cds = cds_sub_done)
remove.experiments(df.rfp) # unload data
# Add some names and convert
names(res) <- libs
data.table::setDT(res) # Will give 1 column per library
res # Now by columns
summary(res)

# orfScore for ribotoolkit pshifts
ofst.tool <- file.path(save.dir, paste0(libs, ".ofst"))
all(file.exists(ofst.tool))
res_tool <- bplapply(ofst.tool, FUN = function(lib, cds) { 
  return(orfScore(cds, fimport(lib), is.sorted = TRUE)$ORFScores)
}, cds = cds_sub_done)
# Add some names and convert
names(res_tool) <- libs
data.table::setDT(res_tool) # Will give 1 column per library
res_tool # Now by columns
summary(res_tool)

orfik.melt <- melt(res)
orfik.melt$source <- "ORFik*"
tool.melt <- melt(res_tool)
tool.melt$source <- "RiboToolkit*"
plot.data.sub <- rbindlist(list(orfik.melt, tool.melt))

plot <- ggplot(data = plot.data.sub) + 
  geom_boxplot(aes(x = variable, y = value, fill = source)) + 
  ylab("ORFscore") + 
  xlab("Library")
plot
ggsave(file.path(save.dir, "ORFscore_comparison_noTIS&TTS.png"), plot, dpi = 300)
ggsave(file.path(save.dir, "ORFscore_comparison_noTIS&TTS.pdf"), plot, dpi = 300)

# Merged plot
plot.data.both <- rbindlist(list(plot.data, plot.data.sub))
plot.data.both[, variable := gsub("RFP_", "", variable)]
fwrite(plot.data.both, file = file.path(save.dir, "Comparison_table.csv"))
plot.data.both[, variable := factor(variable, levels = c("WT_r1", "WT_r2", "WT_r3", "CO_r1", "CO_r2", "CO_r3"), ordered = TRUE)]
# plot.data.both <- fread(file.path(save.dir, "Comparison_table.csv"))
plot <- ggplot(data = plot.data.both) + 
  geom_boxplot(aes(x = variable, y = value, fill = source),
               outlier.size = 0.1) + 
  ylab("ORFscore") + 
  xlab("Library") + 
  ylim(c(-12, 12)) + 
  scale_fill_manual(values = c("#e56598", "#e56571", "#7da1d4" ,"#7da1a1"))
plot

ggsave(file.path(save.dir, "ORFscore_comparison_both.png"), plot, dpi = 300, width = 6, height = 4)
ggsave(file.path(save.dir, "ORFscore_comparison_both.pdf"), plot, dpi = 300, width = 6, height = 4)
