
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Libraries needed
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
library(uORFomePipe)
# devtools::document(pkg = "/export/valenfs/projects/Hakon/uORFomePipe/")
#devtools::install("/export/valenfs/projects/Hakon/uORFomePipe/")

{ # This part will vary according to what your experiments looks like, here I pick 4 stages to use
  # Load experiments
  df.cage <- NULL
  df.rfp  <- read.experiment("zf_baz14_RFP") # RNA-seq is optional, but makes results better
  df.rna  <- read.experiment("zf_baz14_RNA")
  mainPath_aCDS <-"/export/valenfs/projects/Hakon/ORFik_paper/artificial_500_new"
  startCodons = c("ATG", "CTG", "TTG", "AAG", "AGG")
  stopCodons = c("TAA", "TGA", "TAG")
}

# getCDSTraining <- uORFomePipe:::getCDSTraining; getSpecialThreeUTRs <- uORFomePipe:::getSpecialThreeUTRs; getThreeUTRs <- uORFomePipe:::getThreeUTRs; getCDSFiltered <- uORFomePipe:::getCDSFiltered; getGeneralRiboFeatures <- uORFomePipe:::getGeneralRiboFeatures; mainPath = mainPath_aCDS; organism = "Danio rerio"; startCodons.cds.allowed = startCodons; stopCodons.cds.allowed = stopCodons; mode = "aCDS"; max.artificial.length = 100

find_uORFome(mainPath = mainPath_aCDS,
             organism = "Danio rerio", # <- scientific name for organism, will let you know if you misspelled
             df.rfp, df.rna, df.cage,
             startCodons, stopCodons,
             mode = "aCDS")

find_uORFome(mainPath = "/export/valenfs/projects/Hakon/ORFik_paper/artificial_100_new",
             organism = "Danio rerio", # <- scientific name for organism, will let you know if you misspelled
             df.rfp, df.rna, df.cage,
             startCodons, stopCodons,
             mode = "aCDS")
find_uORFome(mainPath = "/export/valenfs/projects/Hakon/ORFik_paper/artificial_300_new",
             organism = "Danio rerio", # <- scientific name for organism, will let you know if you misspelled
             df.rfp, df.rna, df.cage,
             startCodons, stopCodons, max.artificial.length = 300,
             mode = "aCDS")

dt2 <- verification.conf.matrix(artificial = "/export/valenfs/projects/Hakon/ORFik_paper/artificial_100_new",
                                output = 1)


verify <- rbind(cbind(y = dt2$pred_full, model = "Full CDS", readRDS(paste0("features/PredictionData_","combined", "verify",".rds"))),
                cbind(y = dt2$pred_art, model = "art. CDS", readRDS(paste0("features/PredictionData_","combined",".rds"))))
t <- rbind(cbind(model = "Training", readRDS(paste0("features/TrainingData_","combined",".rds"))), verify)
verify$model.fn <- verify$model
verify[verify$model == "art. CDS" & dt2$false.negative,]$model.fn <- "art. CDS false negative"
verify.fn <- verify[rep(dt2$false.negative, 2), ]
verify.tp <- verify[rep(dt2$true.positive, 2), ]

# Positive model for all 3
gg <- ggdensity(t[y == 1 | y == TRUE,], x = colnames(t)[-c(1,2,3)],
                add = "median", rug = TRUE,
                color = "model", fill = "model",
                palette = c("#00AFBB", "#E7B800", "green"), combine = TRUE,
                repel = TRUE, xlim = c(0.01, 5000), title = "Prediction model features: positive",
                subtitle = p("samples: ", nrow(t))) + xscale("log10", .format = TRUE)
ggsave("/export/valenfs/projects/Hakon/ORFik_paper/artificial_100_new/training_model_pos.png", gg, width = 7.64, height = 5.51)
# Negative model for all 3
gg <- ggdensity(t[y == 0 | y == FALSE,], x = colnames(t)[-c(1,2)],
                add = "median", rug = TRUE,
                color = "model", fill = "model",
                palette = c("#00AFBB", "#E7B800", "green"), combine = TRUE,
                repel = TRUE, xlim = c(0.01, 5000), title = "Prediction model features: negative",
                subtitle = p("samples: ", nrow(t))) + xscale("log10", .format = TRUE)
ggsave("/export/valenfs/projects/Hakon/ORFik_paper/artificial_100_new/training_model_neg.png", gg, width = 7.64, height = 5.51)
# FALSE negative of aCDS vs same CDS
gg <- ggdensity(verify.fn, x = colnames(t)[!(colnames(t) %in% c("y", "model", "model.fn"))],
                add = "median", rug = TRUE,
                color = "model", fill = "model",
                palette = c("#00AFBB", "#E7B800"), combine = TRUE,
                repel = TRUE, xlim = c(0.01, 5000), title = "Prediction model features: false.negative",
                subtitle = p("samples: ", nrow(t))) + xscale("log10", .format = TRUE)
ggsave("/export/valenfs/projects/Hakon/ORFik_paper/artificial_100_new/training_model_fn.png", gg, width = 7.64, height = 5.51)
# TRUE positive of aCDS vs same CDS
gg <- ggdensity(verify.tp, x = colnames(t)[!(colnames(t) %in% c("y", "model", "model.fn"))],
                add = "median", rug = TRUE,
                color = "model", fill = "model",
                palette = c("#00AFBB", "#E7B800"), combine = TRUE,
                repel = TRUE, xlim = c(0.01, 5000), title = "Prediction model features: true.positive",
                subtitle = p("samples: ", nrow(t))) + xscale("log10", .format = TRUE)
ggsave("/export/valenfs/projects/Hakon/ORFik_paper/artificial_100_new/training_model_tp.png", gg, width = 7.64, height = 5.51)

gg <- ggdensity(t[y == 1,], x = colnames(t)[-1],
                add = "median", rug = TRUE,
                color = "y", fill = "y",
                palette = c("#00AFBB", "#E7B800"), combine = TRUE,
                repel = TRUE, xlim = c(0.01, 5000), title = "Full length CDS model features",
                subtitle = p("samples: ", nrow(t), " , TRUE: ", sum(t$y))) + xscale("log10", .format = TRUE)
ggsave("/export/valenfs/projects/Hakon/ORFik_paper/artificial_100_new/CDS_model.png", gg, width = 7.64, height = 5.51)


gg <- ggdensity(t[y == 1,], x = colnames(t)[-1],
                add = "median", rug = TRUE,
                color = "y", fill = "y",
                palette = c("#00AFBB", "#E7B800"), combine = TRUE,
                repel = TRUE, xlim = c(0.01, 5000), title = "Artificial length CDS model features",
                subtitle = p("samples: ", nrow(t), " , TRUE: ", sum(t$y))) + xscale("log10", .format = TRUE)
ggsave("/export/valenfs/projects/Hakon/ORFik_paper/artificial_100_new/aCDS_model.png", gg, width = 7.64, height = 5.51)
# uORF test after CDS validation
find_uORFome(mainPath = "/export/valenfs/projects/Hakon/ORFik_paper/uorfs_bazzini",
             organism = "Danio rerio", # <- scientific name for organism, will let you know if you misspelled
             df.rfp, df.rna, df.cage,
             startCodons, stopCodons)

