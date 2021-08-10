library("SingleCellExperiment")
library("scRNAseq")

# Mini muestreo del set de datos usado en: https://bioconductor.org/books/release/OSCA/zeisel-mouse-brain-strt-seq.html#introduction-5

archivo_cuentas <- "https://raw.githubusercontent.com/emarquezz/minidataset_osca/main/min_sce.csv"
archivo_rowData <- "https://raw.githubusercontent.com/emarquezz/minidataset_osca/main/rowD.csv"
archivo_colData <- "https://raw.githubusercontent.com/emarquezz/minidataset_osca/main/colD.csv"


counts <- read.csv(archivo_cuentas, row.names = 1, header = TRUE, check.names = F)
col.data <- DataFrame(read.csv(archivo_colData, row.names = 1, header = TRUE, check.names = F))
row.data <- read.csv(archivo_rowData, row.names = 1, header = TRUE, check.names = F)


sce <- SingleCellExperiment(
      assays = list(counts = counts),
      colData = col.data,
      rowData = row.data)
sce

int_gen <- c("Angpt1", "Chic2", "Mir503", "Magee2", "Nenf", "Eps15l1", "Hsf2bp", "Gnptg", "Vegfb", "Atmin", "Gad1", "Gad2", "Slc32a1", "Dner", "Slc2a13", "Slc6a1", "Nrxn3")

bool_data <- rownames(sce) %in% int_gen
min_sce <- sce[bool_data,]

tej_int <- min_sce$level1class == "interneurons"
tej_pyr <- min_sce$level1class == "pyramidal CA1"
tej_min_sce <- min_sce[,tej_int | tej_pyr]



library("scater")
tej_min_sce <- scater::logNormCounts(tej_min_sce)
plotHeatmap(object = tej_min_sce, features = rownames(tej_min_sce), order_columns_by = "level1class")
