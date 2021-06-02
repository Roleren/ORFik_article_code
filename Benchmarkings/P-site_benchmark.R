#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# P-site detection, benchmark ORFik, shoelaces, Ribotoolkit (which uses plastid)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Data: Alexaki et al 2020
# Steps are:
# 1. Make shift in RiboProfiling
# 2. Load shifts from ORFik, Ribotoolkit, shoelaces and Riboprofiling
# 3. Also pshift libraries from Ribotoolkit, shoelaces and RiboProfiling
# 4. Load data and run the ORFik::orfScore function
# 5. plot the results as boxplots

library(ORFik); library(data.table); library(ggplot2)
df.rfp <- read.experiment("Alexaki_Human_Ribo-Seq")
libs <- bamVarName(df.rfp)
save.dir <- file.path(dirname(df.rfp$filepath[1]), "QC_STATS")
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Make RiboProfiling shifts
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
library(RiboProfiling)
rprofiling.shift.path <- file.path(dirname(df.rfp$filepath[1]), 
                                   "QC_STATS", "RiboProfiling_shifts.csv")
if (!file.exists(rprofiling.shift.path)) {
  outputLibs(df.rfp, type = "default")
  txdb <- loadTxdb(df.rfp) #txdb object with annotations
  shifts.rprofiling1 <- c()
  for(i in libs) {
    aln <- get(i)
    
    alnGRanges <- readsToStartOrEnd(aln, what="start")
    
    oneBinRanges <- aroundPromoter(txdb, alnGRanges, percBestExpressed=0.001)
    #the coverage in the TSS flanking region for the reads with match sizes 29:31
    listPromoterCov <-readStartCov(alnGRanges,
                                   oneBinRanges,
                                   matchSize=c(21, 27,28),
                                   fixedInterval=c(-20, 20),
                                   renameChr="aroundTSS",
                                   charPerc="perc")
    #plotSummarizedCov(listPromoterCov)
    res <- (unlist(listPromoterCov))[2:4]
    max.pos <- lapply(res, function(iSumCov) { # From RiboProfiling internal
      maxPeak <- max(iSumCov$values)
      maxPeakPos <- start(iSumCov)[which(iSumCov$values == 
                                           maxPeak)][1]
    })
    max.pos <- unlist(max.pos, use.names = TRUE)
    shifts.rprofiling1 <- c(shifts.rprofiling1, max.pos)
  }
  stopifnot(length(shifts.rprofiling1) == 18) # 3*6
  shifts_rprofiling <- data.frame(length = c(21, 27, 28),
                             offset = shifts.rprofiling1,
                             library = c(rep("CO1", 3),
                                         rep("CO2", 3),
                                         rep("CO3", 3),
                                         rep("WT1", 3),
                                         rep("WT2", 3),
                                         rep("WT3", 3)))
  fwrite(shifts_rprofiling, rprofiling.shift.path)
} else shifts_rprofiling <- fread(rprofiling.shift.path, header = TRUE)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Get all offsets
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
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

# shoelaces p-sites extracted from running shoelaces GUI version (Linux)
# The offsets were extracted and written down and added here:
shifts_shoelaces <- data.frame(length = c(21, 27, 28),
                                 offset = c(c(9, 9, 12),
                                            c(9, 9, 12),
                                            c(9, 11, 12),
                                            c(9, 11, 12),
                                            c(9, 11, 12),
                                            c(9, 11, 12)),
                                 library = c(rep("CO1", 3),
                                             rep("CO2", 3),
                                             rep("CO3", 3),
                                             rep("WT1", 3),
                                             rep("WT2", 3),
                                             rep("WT3", 3)))

# Sanity checks
all(shifts_orfik$length == shifts_ribotoolkit$length)
all(shifts_orfik$library == shifts_ribotoolkit$library)
all(shifts_orfik$length == shifts_shoelaces$length)
all(shifts_orfik$library == shifts_shoelaces$library)
all(shifts_orfik$length == shifts_rprofiling$length)
all(shifts_orfik$library == shifts_rprofiling$library)
# Now make difference
offset_difference <- data.frame(library = shifts_orfik$library,
                                length = shifts_orfik$length,
                                offset_orfik = shifts_orfik$offset,
                                offset_ribotoolkit = -shifts_ribotoolkit$offset,
                                offset_shoelaces = -shifts_shoelaces$offset,
                                offset_rprofiling = shifts_rprofiling$offset)
offset_difference$All_In_Frame <- ((shifts_orfik$offset + shifts_ribotoolkit$offset) %% 3 == 0) & 
                                    ((shifts_orfik$offset + shifts_shoelaces$offset) %% 3 == 0) & 
                                    ((shifts_orfik$offset + shifts_rprofiling$offset) %% 3 == 0)
  
fwrite(offset_difference, file.path(dirname(df.rfp$filepath[1]), "QC_STATS", "benchmark_pshifting.csv"))
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Make pshifted data for ribotoolkit, Shoelaces, RiboProfiling offsets
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
outputLibs(df.rfp, type = "ofst") # non pshifted files

save.dir.ribo <- file.path(dirname(filepath(df.rfp,
                                       type = "default")[1]),
                      "pshifted_ribotoolkit")
if (!dir.exists(save.dir.ribo)) { # Do only once
  dir.create(save.dir.ribo)
  j = 1
  for(i in unique(shifts_ribotoolkit$library)) {
    shifts <- shifts_ribotoolkit[shifts_ribotoolkit$library %in% i,]
    colnames(shifts) <- c("fraction", "offsets_start", "lib")
    shifts <- shifts[!is.na(shifts$offsets_start),]
    shifts$offsets_start <- -shifts$offsets_start
    if(nrow(shifts) == 0) continue
    gr <- shiftFootprints(footprints = get(libs[j]), shifts)
    export.ofst(gr, file.path(save.dir.ribo, paste0(libs[j], ".ofst")))
    j = j + 1
  }
}
remove.experiments(df.rfp)

# Shoelaces
outputLibs(df.rfp, type = "ofst") # non pshifted files
save.dir.shoe <- file.path(dirname(filepath(df.rfp,
                                            type = "default")[1]),
                           "pshifted_shoelaces")
if (!dir.exists(save.dir.shoe)) { # Do only once
  dir.create(save.dir.shoe)
  j = 1
  for(i in unique(shifts_shoelaces$library)) {
    shifts <- shifts_shoelaces[shifts_shoelaces$library %in% i,]
    colnames(shifts) <- c("fraction", "offsets_start", "lib")
    shifts <- shifts[!is.na(shifts$offsets_start),]
    shifts$offsets_start <- -shifts$offsets_start
    if(nrow(shifts) == 0) continue
    gr <- shiftFootprints(footprints = get(libs[j]), shifts)
    export.ofst(gr, file.path(save.dir.shoe, paste0(libs[j], ".ofst")))
    j = j + 1
  }
  remove.experiments(df.rfp)
}

# RiboProfiling
outputLibs(df.rfp, type = "ofst") # non pshifted files
save.dir.rprof <- file.path(dirname(filepath(df.rfp,
                                            type = "default")[1]),
                           "pshifted_RiboProfiling")
if (!dir.exists(save.dir.rprof)) { # Do only once
  dir.create(save.dir.rprof)
  j = 1
  for(i in unique(shifts_rprofiling$library)) {
    shifts <- shifts_rprofiling[shifts_rprofiling$library %in% i,]
    colnames(shifts) <- c("fraction", "offsets_start", "lib")
    shifts <- shifts[!is.na(shifts$offsets_start),]
    shifts$offsets_start <- shifts$offsets_start
    if(nrow(shifts) == 0) continue
    gr <- shiftFootprints(footprints = get(libs[j]), shifts)
    export.ofst(gr, file.path(save.dir.rprof, paste0(libs[j], ".ofst")))
    j = j + 1
  }
  remove.experiments(df.rfp)
}


#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Check ORFscores for the four ways to shift
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Filter out cds with fpkm < 1
counts <- countTable(df.rfp, "cds", type = "fpkm")
names <- rownames(counts)
names <- names[rowMins(as.matrix(counts)) >= 1]
cds <- loadRegion(df.rfp, "cds", names.keep = names)

# orfScore for ORFik pshifts
outputLibs(df.rfp, type = "pshifted")
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
ofst.tool <- file.path(save.dir.ribo, paste0(libs, ".ofst"))
stopifnot(all(file.exists(ofst.tool)))
res_tool <- bplapply(ofst.tool, FUN = function(lib, cds) {
  return(orfScore(cds, fimport(lib), is.sorted = TRUE)$ORFScores)
}, cds = cds)
# Add some names and convert
names(res_tool) <- libs
data.table::setDT(res_tool) # Will give 1 column per library
res_tool # Now by columns
summary(res_tool)

# orfScore for shoelaces pshifts
ofst.shoe <- file.path(save.dir.shoe, paste0(libs, ".ofst"))
stopifnot(all(file.exists(ofst.shoe)))
res_shoe <- bplapply(ofst.shoe, FUN = function(lib, cds) {
  return(orfScore(cds, fimport(lib), is.sorted = TRUE)$ORFScores)
}, cds = cds)
# Add some names and convert
names(res_shoe) <- libs
data.table::setDT(res_shoe) # Will give 1 column per library
res_shoe # Now by columns
summary(res_shoe)

# orfScore for RiboProfiling pshifts
ofst.rprof <- file.path(save.dir.rprof, paste0(libs, ".ofst"))
stopifnot(all(file.exists(ofst.rprof)))
res_rprof <- bplapply(ofst.rprof, FUN = function(lib, cds) {
  return(orfScore(cds, fimport(lib), is.sorted = TRUE)$ORFScores)
}, cds = cds)
# Add some names and convert
names(res_rprof) <- libs
data.table::setDT(res_rprof) # Will give 1 column per library
res_rprof # Now by columns
summary(res_rprof)

orfik.melt <- melt(res)
orfik.melt$source <- "ORFik"
tool.melt <- melt(res_tool)
tool.melt$source <- "RiboToolkit"
shoe.melt <- melt(res_shoe)
shoe.melt$source <- "Shoelaces"
rprof.melt <- melt(res_rprof)
rprof.melt$source <- "RiboProfiling"
plot.data <- rbindlist(list(orfik.melt, tool.melt, shoe.melt, rprof.melt))

plot.data[, variable := gsub("RFP_", "", variable)]
fwrite(plot.data, file = file.path(save.dir, "Comparison_table_full.csv"))
# plot.data <- fread(file.path(save.dir, "Comparison_table_full.csv"))
plot.data[, variable := factor(variable, levels = c("WT_r1", "WT_r2", "WT_r3", "CO_r1", "CO_r2", "CO_r3"), ordered = TRUE)]

plot <- ggplot(data = plot.data) +
  geom_boxplot(aes(x = variable, y = value, fill = source),
               outlier.size = 0.05) +
  ylab("ORFscore") +
  xlab("Library") +
  ylim(c(-12, 12)) +
  scale_fill_manual(values = c("#e56598", "#7da1d4" ,
                               "#FF9100",  "#31906D"))
plot
ggsave(file.path(save.dir, "ORFscore_comparison.png"), plot, dpi = 300, width = 6, height = 4)
ggsave(file.path(save.dir, "ORFscore_comparison.pdf"), plot, dpi = 300, width = 6, height = 4)

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
res; summary(res) # Now by columns

# orfScore for ribotoolkit pshifts
ofst.tool <- file.path(save.dir.ribo, paste0(libs, ".ofst"))
all(file.exists(ofst.tool))
res_tool <- bplapply(ofst.tool, FUN = function(lib, cds) {
  return(orfScore(cds, fimport(lib), is.sorted = TRUE)$ORFScores)
}, cds = cds_sub_done)
# Add some names and convert
names(res_tool) <- libs
data.table::setDT(res_tool) # Will give 1 column per library
res_tool; summary(res_tool) # Now by columns

# orfScore for shoelaces pshifts
ofst.shoe <- file.path(save.dir.shoe, paste0(libs, ".ofst"))
stopifnot(all(file.exists(ofst.shoe)))
res_shoe <- bplapply(ofst.shoe, FUN = function(lib, cds) {
  return(orfScore(cds, fimport(lib), is.sorted = TRUE)$ORFScores)
}, cds = cds_sub_done)
# Add some names and convert
names(res_shoe) <- libs
data.table::setDT(res_shoe) # Will give 1 column per library
res_shoe; summary(res_shoe) # Now by columns


# orfScore for RiboProfiling pshifts
ofst.rprof <- file.path(save.dir.rprof, paste0(libs, ".ofst"))
stopifnot(all(file.exists(ofst.rprof)))
res_rprof <- bplapply(ofst.rprof, FUN = function(lib, cds) {
  return(orfScore(cds, fimport(lib), is.sorted = TRUE)$ORFScores)
}, cds = cds_sub_done)
# Add some names and convert
names(res_rprof) <- libs
data.table::setDT(res_rprof) # Will give 1 column per library
res_rprof; summary(res_rprof) # Now by columns


orfik.melt <- melt(res)
orfik.melt$source <- "ORFik*"
tool.melt <- melt(res_tool)
tool.melt$source <- "RiboToolkit*"
shoe.melt <- melt(res_shoe)
shoe.melt$source <- "Shoelaces*"
rprof.melt <- melt(res_rprof)
rprof.melt$source <- "RiboProfiling*"
plot.data.sub <- rbindlist(list(orfik.melt, tool.melt, shoe.melt, rprof.melt))

plot <- ggplot(data = plot.data.sub) +
  geom_boxplot(aes(x = variable, y = value, fill = source)) +
  ylab("ORFscore") +
  xlab("Library")
plot
ggsave(file.path(save.dir, "ORFscore_comparison_noTIS&TTS.png"), plot, dpi = 300)
ggsave(file.path(save.dir, "ORFscore_comparison_noTIS&TTS.pdf"), plot, dpi = 300)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Merge and plot
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Merged plot
plot.data.both <- rbindlist(list(plot.data, plot.data.sub))
plot.data.both[, variable := gsub("RFP_", "", variable)]
fwrite(plot.data.both, file = file.path(save.dir, "Comparison_table.csv"))

plot.data.both[, variable := factor(variable, levels = c("WT_r1", "WT_r2", "WT_r3", "CO_r1", "CO_r2", "CO_r3"), ordered = TRUE)]
# plot.data.both <- fread(file.path(save.dir, "Comparison_table.csv"))
plot <- ggplot(data = plot.data.both) +
  geom_boxplot(aes(x = variable, y = value, fill = source),
               outlier.size = 0.05) +
  ylab("ORFscore") +
  xlab("Library") +
  ylim(c(-12, 12)) +
  scale_fill_manual(values = c("#e56598", "#e56571", "#7da1d4" ,"#7da1a1",
                               "#FF9100", "#FF9F50", "#31906D", "#3190AD"))
plot

ggsave(file.path(save.dir, "ORFscore_comparison_both.png"), plot, dpi = 300, width = 6, height = 4)
ggsave(file.path(save.dir, "ORFscore_comparison_both.pdf"), plot, dpi = 300, width = 6, height = 4)

