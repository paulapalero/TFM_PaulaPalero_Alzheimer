################################################################################
# RNAAgeCalc EN DATOS SINGLE-CELL
# Prueba en los datos del estudio GSE157827
################################################################################

################################################################################
# LIBRERÍAS
################################################################################

library(recount)
library(SummarizedExperiment)
library(RNAAgeCalc)
library(ggplot2)
library(Matrix)
library(dplyr)
library(broom)

################################################################################
# CARGA DE DATOS
################################################################################

sce.annot <- readRDS(
  "/clinicfs/userhomes/ppalero/RScripts/sce.annot.clean.rds")

sce <- sce.annot

################################################################################
# EXPLORACIÓN INICIAL DEL OBJETO
################################################################################

assayNames(sce.annot)        # matrices disponibles
dim(sce.annot)               # genes x células
colnames(colData(sce.annot)) # metadatos
rownames(sce.annot)[1:10]    # genes
head(colData(sce.annot))

# PatientID      -> identificador del paciente
# AGE            -> edad cronológica
# SingleR.labels -> tipo celular
# cell.type      -> utilizar si se emplea el objeto "clean"

################################################################################
# FILTRADO DE CÉLULAS
################################################################################

# Eliminar células descartadas
sce <- sce[, sce$discard == FALSE]

# Eliminar células sin anotación
sce <- sce[, sce$SingleR.labels != "unknown"]

################################################################################
# GENERACIÓN DE PSEUDOBULKS
################################################################################

# Identificador único:
# paciente + tipo celular

group_id <- paste(
  sce$PatientID,
  sce$SingleR.labels,
  sep = "_"
)

# Matriz de conteos
counts_mat <- counts(sce)

# Agregación pseudobulk
pb <- rowsum(t(counts_mat), group_id)
pb <- t(pb)

dim(pb)
head(colnames(pb))
head(rownames(pb))

# rownames(pb) = IDs ENSEMBL

################################################################################
# CONSTRUCCIÓN DEL DATAFRAME DE EDAD CRONOLÓGICA
################################################################################

# Extraer edad por paciente
meta_df <- as.data.frame(colData(sce))

patient_age <- unique(meta_df[, c("PatientID", "AGE")])

# Crear tabla de muestras pseudobulk
pb_samples <- data.frame(
  sample_id = colnames(pb)
)

pb_samples$PatientID <- sub(
  "_.*",
  "",
  pb_samples$sample_id
)

# Asociar edad cronológica

chronage <- merge(
  pb_samples,
  patient_age,
  by = "PatientID",
  all.x = TRUE
)

# Formato requerido por RNAAgeCalc

chronage_final <- chronage[, c("sample_id", "AGE")]

colnames(chronage_final) <- c(
  "sample_id",
  "age"
)

all(chronage_final$sample_id == colnames(pb))

################################################################################
# PREDICCIÓN DE EDAD BIOLÓGICA
################################################################################

res <- predict_age(
  exprdata  = pb,
  tissue    = "brain",
  exprtype  = "counts",
  idtype    = "ENSEMBL",
  signature = "all",
  chronage  = chronage_final
)

################################################################################
# VISUALIZACIÓN GENERAL DE RNAAgeCalc
################################################################################

makeplot(res)

################################################################################
# PREPARACIÓN DE TABLA DE RESULTADOS
################################################################################

res_df <- data.frame(
  sample_id     = rownames(res),
  predicted_age = res$RNAAge,
  chron_age     = res$ChronAge,
  age_accel     = res$AgeAccelResid
)

# Extraer paciente

res_df$PatientID <- sub(
  "_.*",
  "",
  res_df$sample_id
)

# Extraer tipo celular

res_df$celltype <- sub(
  ".*_",
  "",
  res_df$sample_id
)

################################################################################
# AÑADIR METADATOS CLÍNICOS
################################################################################

meta_patient <- unique(
  as.data.frame(colData(sce))[
    ,
    c("PatientID", "CONDITION", "SEX")
  ]
)

res_df <- merge(
  res_df,
  meta_patient,
  by = "PatientID"
)

################################################################################
# BOXPLOTS DE ACELERACIÓN DE EDAD
################################################################################

# AD vs Control

ggplot(
  res_df,
  aes(
    x = CONDITION,
    y = age_accel,
    fill = CONDITION
  )
) +
  geom_boxplot() +
  geom_jitter(width = 0.15) +
  facet_wrap(~celltype, scales = "free_y") +
  theme_bw()

################################################################################
# DIFERENCIAS POR SEXO
################################################################################

ggplot(
  res_df,
  aes(
    x = SEX,
    y = age_accel,
    fill = SEX
  )
) +
  geom_boxplot() +
  facet_wrap(~celltype, scales = "free_y") +
  theme_bw()

################################################################################
# INTERACCIÓN CONDICIÓN × SEXO
################################################################################

ggplot(
  res_df,
  aes(
    x = CONDITION,
    y = age_accel,
    fill = SEX
  )
) +
  geom_boxplot(
    position = position_dodge(0.8)
  ) +
  facet_wrap(~celltype) +
  theme_bw()

################################################################################
# ANÁLISIS ESTADÍSTICO
################################################################################

# Efecto de la condición

res_df %>%
  group_by(celltype) %>%
  do(
    tidy(
      lm(age_accel ~ CONDITION, data = .)
    )
  )

################################################################################
# EFECTO CONDICIÓN × SEXO
################################################################################

res_df %>%
  group_by(celltype) %>%
  do(
    tidy(
      lm(age_accel ~ CONDITION * SEX, data = .)
    )
  )

################################################################################
# EDAD CRONOLÓGICA VS EDAD BIOLÓGICA
################################################################################

ggplot(
  res_df,
  aes(
    x = chron_age,
    y = predicted_age,
    color = celltype
  )
) +
  geom_point(
    alpha = 1,
    size = 3
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    color = "black"
  ) +
  theme_classic() +
  labs(
    x = "Edad cronológica",
    y = "Edad biológica (RNAAgeCalc)",
    title = "Edad biológica vs cronológica por tipo celular"
  )

################################################################################
# ACELERACIÓN DE EDAD (ANÁLISIS ALTERNATIVO)
################################################################################

rna_age_results$age_acceleration <-
  rna_age_results$predicted_age -
  rna_age_results$chronological_age

ggplot(
  rna_age_results,
  aes(
    x = cell_type,
    y = age_acceleration,
    fill = diagnosis
  )
) +
  geom_boxplot() +
  theme_classic() +
  labs(
    y = "Aceleración de edad (RNAAgeCalc)"
  )


