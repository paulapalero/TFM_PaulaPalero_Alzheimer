################################################################################
# SCRIPT PARA INTEGRAR LOS METADATOS 
# DEL ESTUDIO SYN18485171
################################################################################

################################################################################
# LIBRERIAS
################################################################################

library(Seurat)
library(dplyr)
library(S4Vectors)
library(Matrix)
library(tibble)

################################################################################
# CARGAR DATOS
################################################################################

setwd('/clinicfs/userhomes/ppalero/SYN52293442/syn52293442_pfc')

cortex <- readRDS('data/Prefrontal_cortex.rds')
class(cortex)
dim(cortex)

################################################################################
# GENERAR METADATOS
# A PARTIR DE LAS TABLAS SUPLEMENTARIAS DEL ESTUDIO
################################################################################

metadata <- read.table("metadata/Supplementary_Table_1_sample_metadata.txt", 
                       sep = "\t", header = T)
metadata <- metadata[metadata$region == "PFC",]

metadata7 <- read.table("metadata/MIT_ROSMAP_Multiomics_individual_metadata.csv", 
                        sep = ",", header = T)
rosmap_metadata <- read.table("metadata/ROSMAP_clinical.csv", sep = ",", 
                              header = T)

metadata <- merge(metadata, unique(metadata7[ ,c(1,19)]),
                  by.x = 'subject', by.y = 'subject', all.x = TRUE)

metadata <- merge(metadata, rosmap_metadata, 
                  by.x = 'individualID', by.y = 'individualID')

write.table(metadata, file = "metadata/final_metadata.txt")

#----------------------- Guardamos los metadatos -------------------------------

coldata <- cortex@meta.data %>%
  rownames_to_column("cell_id")   # guardamos IDs de célula (los que aparece
                                  # en el seurat 'cortex' como rownames)

coldata <- coldata %>%
  left_join(metadata, by = "projid")

#----------------------- Creamos el colData ------------------------------------

# Crear columna 'sex'
coldata <- coldata %>%
  mutate(sex = recode(msex.x,
                      '0' = 'Mujer',
                      '1' = 'Hombre'))

# Crear columna 'condition'
coldata <- coldata %>%
  mutate(condition = recode(pathAD,
                            'non-AD' = 'Control',
                            'AD' = 'AD'))
# Crear columna de covariable
coldata <- coldata %>%
  mutate(group.sex = paste(condition, sex, sep = '_'))

#----------------------- Restaurar rownames ------------------------------------

coldata <- coldata %>%
  column_to_rownames("cell_id")

# Reordenar correctamente
coldata <- coldata[colnames(cortex), ]

# Comprobación
stopifnot(all(rownames(coldata) == colnames(cortex)))

#----------------------- Asignamos de nuevo al Seurat --------------------------

cortex@meta.data <- coldata

#------------------------------- Guardamos -------------------------------------

saveRDS(cortex, 'results/seurat.raw.rds')

