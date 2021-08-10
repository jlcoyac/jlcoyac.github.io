library("scRNAseq")
sce.416b <- LunSpikeInData(which = "416b")

# Carga el paquete SingleCellExperiment
library("SingleCellExperiment")

"Primera parte
aquí checamos el slot assays"
# Extrae la matriz de cuentas del set de datos de 416b
counts.416b <- counts(sce.416b)
# CHEQUEMOS clase y dimensiones
class(counts.416b) # es matriz
dim(counts.416b) # indicará genes y células

# CONSTRUIR un nuevo objeto SCE de la matriz de cuentas !!!!!!
sce <- SingleCellExperiment(assays = list(counts = counts.416b))

# Revisa el objeto que acabamos de crear
sce

## ¿Qué tan grande es el objeto de R?
lobstr::obj_size(sce) / 1024^2 ## En MB


# Accesa la matriz de cuenta del compartimento (slot) "assays"
# assays(sce, "counts")
# OJO: ¡esto puede inundar tu sesión de R!

# 1. El método general
assay(sce, "counts")[110:115, 1:3] # gene, cell


# 2. El método específico para accesar la matriz de cuentas "counts"
counts(sce)[110:115, 1:3]


# AGREGAR MAS ASSAYS
sce <- scater::logNormCounts(sce)
# Revisa el objeto que acabamos de actualizar
sce


## ¿Qué tan grande es el objeto de R?
lobstr::obj_size(sce) / 1024^2 ## En MB
# 1. El método general
assay(sce, "logcounts")[110:115, 1:3]

#2. El método específico para accesar la matriz de cuentas
# transformadas "logcounts"
logcounts(sce)[110:115, 1:3]

# agregemos un assay mas, esta vez de manera manual
assay(sce, "counts_100") <- assay(sce, "counts") + 100 # suma 100 a counts assay
# Enumera los "assays" en el objeto
assays(sce) # indica num y nombre de assays

assayNames(sce) # solo nos dará los nombres de los assays


#      assay(sce, "counts_100")[110:115, 1:3]
## ¿Qué tan grande es el objeto de R?
lobstr::obj_size(sce) / 1024^2 ## En MB

"segunda parte:
aquí checaremos metadata de las células"
## [1] "segunda parte:\naquí checaremos metadata de las células"

# Extrae la información de las muestras (metadata) del set de datos de 416b
colData.416b <- colData(sce.416b) # podemos checar objeto en la cajita de environment de RStudio!!

# explorar datooos
table(colData.416b$phenotype)

table(colData.416b$block) # fue en varios dias?

# Agrega algo de esa información a nuestro objeto de SCE
colData(sce) <- colData.416b[, c("phenotype", "block")]
# Revisa el objeto que acabamos de actualizar
sce


# Accesa a la información de las muestras (metadata) en nuestro SCE
colData(sce) # usar head?


# Accesa una columna específica de la información de las muestras (metadata)
table(sce$block)

table(colData(sce)$block) # otra manera
# Ejemplo de una función que agrega columnas nuevas al colData
sce <- scater::addPerCellQC(sce.416b) # añade datos de control de calidad
# Accesa a la información de las muestras (metadata) en nuestro SCE actualizado
colData(sce)

# Revisa el objeto que acabamos de actualizar
sce
## ¿Qué tan grande es el objeto de R?
lobstr::obj_size(sce) / 1024^2 ## En MB

## Agrega las cuentas normalizadas (lognorm) de nuevo
sce <- scater::logNormCounts(sce)

## ¿Qué tan grande es el objeto de R?
lobstr::obj_size(sce) / 1024^2 ## En MB

# Ejemplo: obtén el subconjunto de células de fenotipo "wild type"
# Acuérdate que las células son columnas del SCE !!!!
sce[, sce$phenotype == "wild type phenotype"]


"Tercera parte:
examinaremos metadata de features (rowData)"
## [1] "Tercera parte:\nexaminaremos metadata de features (rowData)"

# Accesa la información de los genes de nuestro SCE
# ¡Está vació actualmente!
rowData(sce)

# Ejemplo de una función que agrega campos nuevos en el rowData
sce <- scater::addPerFeatureQC(sce)
# Accesa a la información de las muestras (metadata) en nuestro SCE actualizado
rowData(sce)

## ¿Qué tan grande es el objeto de R?
lobstr::obj_size(sce) / 1024^2 ## En MB


# Descarga los archivos de anotación de la base de datos de Ensembl
# correspondientes usando los recursos disponibles vía AnnotationHub
library("AnnotationHub")
ah <- AnnotationHub()
query(ah, c("Mus musculus", "Ensembl", "v97"))


# Obtén la posición del cromosoma para cada gen
ensdb <- ah[["AH73905"]]
chromosome <- mapIds(ensdb,
                     keys = rownames(sce),
                     keytype = "GENEID",
                     column = "SEQNAME"
)
rowData(sce)$chromosome <- chromosome

# Accesa a la información de las muestras (metadata) en nuestro SCE actualizado
rowData(sce)

## ¿Qué tan grande es el objeto de R?
lobstr::obj_size(sce) / 1024^2 ## En MB


# Ejemplo: obtén el subconjunto de datos donde los genes están en el
# cromosoma 3
# NOTA: which() fue necesario para lidear con los nombres de cromosoma
# que son NA
sce[which(rowData(sce)$chromosome == "3"), ]

"Cuarta parte:
examinamos slot metadata"
## [1] "Cuarta parte:\nexaminamos slot metadata"

# Accesa la información de nuestro experimento usando metadata()
# ¡Está vació actualmente!
metadata(sce)

# La información en el metadata() es como Vegas - todo se vale
metadata(sce) <- list(
  favourite_genes = c("Shh", "Nck1", "Diablo"),
  analyst = c("Pete")
)

# Accesa la información de nuestro experimento usando metadata() de
# nuestro objeto actualizado
metadata(sce)


"Quinta parte:
examinamos slot de reducción de dimensiones"
## [1] "Quinta parte:\nexaminamos slot de reducción de dimensiones"

# Ejemplo: agrega los componentes principales (PCs) de las logcounts
# NOTA: aprenderemos más sobre análisis de componentes principales (PCA) después
sce <- scater::runPCA(sce)
# Revisa el objeto que acabamos de actualizar
sce

# Accesa la matriz de PCA del componente (slot) reducedDims
reducedDim(sce, "PCA")[1:6, 1:3]


# Ejemplo, agrega una representación de los logcounts en t-SNE
# NOTA: aprenderemos más sobre t-SNE después
sce <- scater::runTSNE(sce)
# Revisa el objeto que acabamos de actualizar
sce


# Accesa a la matriz de t-SNE en el componente (slot) de reducedDims
head(reducedDim(sce, "TSNE"))


# Ejemplo: agrega una representación 'manual' de los logcounts en UMAP
# NOTA: aprenderemos más sobre UMAP después y de una forma más sencilla de
#       calcularla
u <- uwot::umap(t(logcounts(sce)), n_components = 2)

# Agrega la matriz de UMAP al componente (slot) reducedDims
reducedDim(sce, "UMAP") <- u

# Accesa a la matriz de UMAP desde el componente (slot) reducedDims
head(reducedDim(sce, "UMAP"))

# Enumera los resultados de reducción de dimensiones en nuestro objeto SCE
reducedDims(sce)

"Sexta parte:
experimentos alternativos"
## [1] "Sexta parte:\nexperimentos alternativos"

# Extrae la información de ERCC de nuestro SCE para el set de datos de 416b
ercc.sce.416b <- altExp(sce.416b, "ERCC")
# Inspecciona el SCE para los datos de ERCC
ercc.sce.416b

# Agrega el SCE de ERCC como un experimento alternativo a nuestro SCE
altExp(sce, "ERCC") <- ercc.sce.416b
# Revisa el objeto que acabamos de actualizar
sce

## ¿Qué tan grande es el objeto de R?
lobstr::obj_size(sce) / 1024^2 ## En MB


# Enumera los experimentos alternativos almacenados en nuestro objeto
altExps(sce)


# El crear un subconjunto del SCE por muestra (célula) automáticamente
# obtiene el subconjunto de los experimentos alternativos
sce.subset <- sce[, 1:10]
ncol(sce.subset)

ncol(altExp(sce.subset))

## ¿Qué tan grande es el objeto de R?
lobstr::obj_size(sce.subset) / 1024^2 ## En MB

"Septima parte:
sizefactores"
## [1] "Septima parte:\nsizefactores"

# Extrae los factores de tamaño (size factors)
# Estos fueron añadidos a nuestro objeto cuando corrimos
# scater::logNormCounts(sce)
head(sizeFactors(sce))

# "Automáticamente" reemplaza los factores de tamaño
sce <- scran::computeSumFactors(sce)
head(sizeFactors(sce))

# "Manualmente" reemplaza los factores de tamaño
sizeFactors(sce) <- scater::librarySizeFactors(sce)
head(sizeFactors(sce))

