## Paquetes de este capítulo
library("scRNAseq") ## para descargar datos de ejemplo
library("AnnotationHub") ## para obtener información de genes
library("scater") ## para gráficas y control de calidad
library("BiocFileCache") ## para descargar datos
library("DropletUtils") ## para detectar droplets
library("Matrix") ## para leer datos en formatos comprimidos


## Datos
library("scRNAseq")
sce.416b <- LunSpikeInData(which = "416b")
library("SingleCellExperiment")

sce.416b$block <- factor(sce.416b$block)

# Descarga los archivos de anotación de la base de datos de Ensembl
# correspondientes usando los recursos disponibles vía AnnotationHub
library("AnnotationHub")
ah <- AnnotationHub()


query(ah, c("Mus musculus", "Ensembl", "v97"))
# Obtén la posición del cromosoma para cada gen
ens.mm.v97 <- ah[["AH73905"]]
location <- mapIds(
  ens.mm.v97,
  keys = rownames(sce.416b),
  keytype = "GENEID",
  column = "SEQNAME"
)


# Identifica los genes mitocondriales
is.mito <- which(location == "MT")
library("scater")
sce.416b <- addPerCellQC(sce.416b,
                         subsets = list(Mito = is.mito)
)
plotColData(sce.416b, x = "block", y = "detected")


