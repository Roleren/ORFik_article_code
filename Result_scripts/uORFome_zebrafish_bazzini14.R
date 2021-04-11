#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# uORFome pipeline for Bazzini et al 2014, with CAGE from nepal et al
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
### This script uses the extension package for ORFik called uORFomePipe, install if wanted.
#devtools::install_github("Roleren/uORFomePipe")
library(uORFomePipe)

# Output folder
mainPath <- "/export/valenfs/projects/Hakon/uORFomes/uORFome_zebrafish_bazzini_test4"
# Input data experiment creation
exp.name.CAGE <- "zf_nepal"
exp.name.RFP <- "zf_baz14_RFP"
exp.name.RNA <- "zf_baz14_RNA"

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# INIT (START HERE)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
{ # Load experiments (subset to 3 stages)
  df.cage <- read.experiment(exp.name.CAGE)
  df.rfp  <- read.experiment(exp.name.RFP)
  df.rna  <- read.experiment(exp.name.RNA)
  # Subset experiments to 12hpf, 24hpf and 48hpf
  conditions <- c("", NA) # Only empty conditions allowed (no mutants etc.)
  # Create the 3 stages: 12hpf, 24hpf & 48hpf
  stages <- c("Somite", "12hpf","24hpf", "prim6", "48hpf", "prim20") # 3 stages total
  df.rfp <- df.rfp[df.rfp$stage %in% stages & df.rfp$condition %in% conditions,]
  df.rna <- df.rna[df.rna$stage %in% stages & df.rna$condition %in% conditions,]
  df.cage <- df.cage[df.cage$stage %in% stages & df.cage$condition %in% conditions,]
  df.cage[3,2] <- df.rna$stage[2]; df.cage[4:5,2] <- df.rna$stage[3]
  df.rna[1,2] <- "12hpf"; df.rfp[1,2] <- "12hpf"; df.cage[1:2,2] <-  df.rna$stage[1]


  # Start and stop codons used
  startCodons = "ATG|CTG|TTG|GTG|AAG|AGG|ACG|ATC|ATA|ATT"
  stopCodons = "TAA|TAG|TGA"

  # Run
  find_uORFome(mainPath = mainPath,
               organism = organism.df(df.rfp),
               df.rfp, df.rna, df.cage,
               startCodons, stopCodons)
}
# Run this line to reactivate database if you restarted session:
uORFomePipe:::createDataBase(file.path(dataBaseFolder, "uorfCatalogue.sqlite"))

prediction <- predictUorfs(mode = "uORFs")
predictionVsCageHits(saveName = paste0("differential_", "uORFs", "_usage.pdf"), prediction = prediction)

# Load results and do validations if you want now
uORFs <- getUorfsInDb()
listTables() # the data in database

pred <- readTable("tissueAtlasByCageAndPred")
colSums(pred)
pred$total <- NULL # Remove total prediction column

grl[rowSums(pred) == 4,] # uORF in all 4 tissues



# Venn diagram custom
uORFomePipe:::venn.diagram.uORFs(width = 4, height = 3)
predictions = readTable("tissueAtlasByCageAndPred")
preds <- copy(predictions)
preds$total <- NULL
preds <- preds[rowSums(preds) > 0, ]

my_list <- list()
for(i in seq_along(preds)) {
  my_list[[length(my_list) + 1]] <- which(preds[, i, with = F] == 1)
}

VennDiagram::venn.diagram(imagetype = "png",
                          x = my_list,
                          category.names = colnames(preds),
                          filename = 'venn_diagramm_uORFs.png',
                          output=TRUE, units = "in",
                          width = 3.3, height = 3.3,
                          cat.dist =c(0.05,0.05,0.05),
                          print.mode = c("raw","percent"))

VennDiagram::venn.diagram(imagetype = "pdf",
                          x = my_list,
                          category.names = colnames(preds),
                          filename = 'venn_diagramm_uORFs.pdf',
                          output=TRUE, units = "in",
                          width = 3.3, height = 3.3,
                          cat.dist =c(0.05,0.05,0.05),
                          print.mode = c("raw","percent"))
