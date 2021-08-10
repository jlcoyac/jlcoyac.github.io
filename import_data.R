# Descarga datos de ejemplo procesados con CellRanger
# Paréntesis: al usar BiocFileCache solo tenemos que descargar
#             los datos una vez.
library("BiocFileCache")
bfc <- BiocFileCache()
pbmc.url <-
  paste0(
    "http://cf.10xgenomics.com/samples/cell-vdj/",
    "3.1.0/vdj_v1_hs_pbmc3/",
    "vdj_v1_hs_pbmc3_filtered_feature_bc_matrix.tar.gz"
  )
pbmc.data <- bfcrpath(bfc, pbmc.url)


## adding rname 'http://cf.10xgenomics.com/samples/cell-vdj/3.1.0/vdj_v1_hs_pbmc3/vdj_v1_hs_pbmc3_filtered_feature_bc_matrix.tar.gz'
# Extrae los archivos en un directorio temporal
untar(pbmc.data, exdir = tempdir())

# Enumera los archivos que descargamos y que extrajimos
# Estos son los archivos típicos de CellRanger
pbmc.dir <- file.path(
  tempdir(),
  "filtered_feature_bc_matrix"
)
list.files(pbmc.dir)

# Importa los datos como un objeto de tipo SingleCellExperiment
library("DropletUtils")
sce.pbmc <- read10xCounts(pbmc.dir)
# Revisa el objeto que acabamos de construir
sce.pbmc

## ¿Qué tan grande es el objeto de R?
lobstr::obj_size(sce.pbmc) / 1024^2 ## En MB

# Almacena la información de CITE-seq como un experimento alternativo
sce.pbmc <- splitAltExps(sce.pbmc, rowData(sce.pbmc)$Type)
# Revisa el objeto que acabamos de actualizar
sce.pbmc

## ¿Qué tan grande es el objeto de R?
lobstr::obj_size(sce.pbmc) / 1024^2 ## En MB

# Descarga datos de ejemplo procesados con scPipe
library("BiocFileCache")
bfc <- BiocFileCache()
sis_seq.url <-
  "https://github.com/LuyiTian/SIS-seq_script/archive/master.zip"
sis_seq.data <- bfcrpath(bfc, sis_seq.url)
## adding rname 'https://github.com/LuyiTian/SIS-seq_script/archive/master.zip'

# Extrae los archivos en un directorio temporal
unzip(sis_seq.data, exdir = tempdir())

# Enumera (algunos de) los archivos que descargamos y extrajimos
# Estos son los archivos típicos de scPipe
sis_seq.dir <- file.path(
  tempdir(),
  "SIS-seq_script-master",
  "data",
  "BcorKO_scRNAseq",
  "RPI10"
)
list.files(sis_seq.dir)

# Importa los datos como un objeto de tipo SingleCellExperiment
library("scPipe")


sce.sis_seq <- create_sce_by_dir(sis_seq.dir)

# Revisa el objeto que acabamos de construir
sce.sis_seq

## ¿Qué tan grande es el objeto de R?
lobstr::obj_size(sce.sis_seq) / 1024^2 ## En MB

# Descarga un ejemplo de un montón de archivos
library("BiocFileCache")
bfc <- BiocFileCache()
lun_counts.url <-
  paste0(
    "https://www.ebi.ac.uk/arrayexpress/files/",
    "E-MTAB-5522/E-MTAB-5522.processed.1.zip"
  )
lun_counts.data <- bfcrpath(bfc, lun_counts.url)

lun_coldata.url <-
  paste0(
    "https://www.ebi.ac.uk/arrayexpress/files/",
    "E-MTAB-5522/E-MTAB-5522.sdrf.txt"
  )
lun_coldata.data <- bfcrpath(bfc, lun_coldata.url)

# Extrae los archivos en un directorio temporal
lun_counts.dir <- tempfile("lun_counts.")
unzip(lun_counts.data, exdir = lun_counts.dir)

# Enumera los archivos que descargamos y extrajimos
list.files(lun_counts.dir)

# Lee la matriz de cuentas (para una placa)
lun.counts <- read.delim(
  file.path(lun_counts.dir, "counts_Calero_20160113.tsv"),
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)
# Almacena la información de la longitud de los genes para después
gene.lengths <- lun.counts$Length
# Convierte los datos de cuentas de genez a una matriz (quitamos las longitudes)
lun.counts <- as.matrix(lun.counts[, -1])

# Lee la información de las muestras (células)
lun.coldata <- read.delim(lun_coldata.data,
                          check.names = FALSE,
                          stringsAsFactors = FALSE
)
library("S4Vectors")
lun.coldata <- as(lun.coldata, "DataFrame")

# Pon en orden la información de las muestras para que
# sea idéntico al orden en la matriz de cuentas
m <- match(
  colnames(lun.counts),
  lun.coldata$`Source Name`
)
lun.coldata <- lun.coldata[m, ]

# Construye la tabla de información de los genes
lun.rowdata <- DataFrame(Length = gene.lengths)

# Construye el objeto de SingleCellExperiment
lun.sce <- SingleCellExperiment(
  assays = list(assays = lun.counts),
  colData = lun.coldata,
  rowData = lun.rowdata
)
# Revisa el objeto que acabamos de construir
lun.sce

## ¿Qué tan grande es el objeto de R?
lobstr::obj_size(lun.sce) / 1024^2 ## En MB

