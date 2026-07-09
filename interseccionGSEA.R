################################################################################
# SCRIPTS AUTOMATIZADO PARA ANALIZAR
# LA INTERSECCIÓN DE LOS RESULTADOS DEL ANALISIS FUNCIONAL
# PARA CADA DATABASE ENTRE GSE157827, SYN18485171 Y SYN52293442
################################################################################

################################################################################
# LIBRERÍAS
################################################################################

library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tidyr)

################################################################################
# CONFIGURAR RUTAS A DIRECTORIOS/CARPETAS
################################################################################
rutas <- list(
  GSE157827_GO       = "/clinicfs/userhomes/ppalero/GSE157827/GSEA/GO",
  GSE157827_Reactome = "/clinicfs/userhomes/ppalero/GSE157827/GSEA/Reactome",
  GSE157827_Aging = "/clinicfs/userhomes/ppalero/GSE157827/GSEA/Aging",
  SYN18485171_GO     = "/clinicfs/userhomes/ppalero/SYN18485171/GSEA/GO",
  SYN18485171_Reactome = "/clinicfs/userhomes/ppalero/SYN18485171/GSEA/Reactome",
  SYN18485171_Aging = "/clinicfs/userhomes/ppalero/SYN18485171/GSEA/Aging",
  SYN52293442_GO     = "/clinicfs/userhomes/ppalero/SYN52293442/GSEA/GO",
  SYN52293442_Reactome = "/clinicfs/userhomes/ppalero/SYN52293442/GSEA/Reactome",
  SYN52293442_Aging = "/clinicfs/userhomes/ppalero/SYN52293442/GSEA/Aging"
)

################################################################################
# FUNCIÓN PRINCIPAL
# CREAR UNA TABLA RESUMEN CON TODA LA INFORMACIÓN DE GSEA
# PARA CADA ESTUDIO SEPARANDO POR DATABASE
################################################################################

crear_tabla <- function(ruta, etiqueta) {
  
  db_type <- case_when(
    grepl("GO", ruta)       ~ "GO",
    grepl("Reactome", ruta) ~ "Reactome",
    grepl("Aging", ruta)    ~ "Aging",
  )
  estudio <- str_extract(etiqueta, "GSE157827|SYN18485171|SYN52293442")
  
  # Listamos todos los archivos que terminan en .csv
  all_files <- list.files(ruta, pattern = "\\.csv$", full.names = TRUE)
  
  if(length(all_files) == 0) return(NULL)
  
  # Combinar
  resumen <- all_files %>% map_df(function(f) {
    df <- read_csv(f, show_col_types = FALSE)
    
    filename <- basename(f)
    
    # Extracción del CellType basada en nombre de archivo:
    # Formato: gseaGO_TipoCelular_Condicion_...
    # Separamos por "_" y tomamos el segundo elemento
    partes <- str_split(filename, "_")[[1]]
    tipo_celular <- partes[2]
    condicion    <- partes[3]
    
    df %>% 
      mutate(
        CellType = tipo_celular,
        Condition = condicion,
        Database = db_type,
        Study = estudio,
      ) %>%
      select(ID, Description, NES, pvalue, p.adjust,
             CellType, Condition, Database, Study, everything())
  })
  
  # Guardar
  output_dir <- "/clinicfs/userhomes/ppalero/ALL"
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  output_name <- paste0("Total_GSEA_", db_type, "_", estudio, ".csv")
  write_csv(resumen, file.path(output_dir, output_name))

  return(resumen)
}

# Ejecución
tablas_finales <- imap(rutas, ~crear_tabla(.x, .y))


################################################################################
# FUNCIÓN PARA RECOPILAR TODOS LOS DATOS DE UN TIPO CELULAR ESPECÍFICO
################################################################################

all_info <- function(rutas, celula_target) {
  map_df(rutas, function(ruta) {
    db_type <- case_when(
      grepl("GO", ruta)       ~ "GO",
      grepl("Reactome", ruta) ~ "Reactome",
      grepl("Aging", ruta)    ~ "Aging",
    )
    files <- list.files(ruta, pattern = paste0("_", celula_target, "_"), 
                        full.names = TRUE)
    
    map_df(files, function(f) {
      df <- read_csv(f, show_col_types = FALSE)
      condicion <- str_extract(basename(f), "AD|Control")
      estudio <- case_when(
        grepl("GSE157827", f) ~ "GSE157827",
        grepl("SYN18485171", f) ~ "SYN18485171",
        grepl("SYN52293442", f) ~ "SYN52293442",
      )
      
      df %>% mutate(Cell_Type = celula_target, Condition = condicion, 
                    Database = db_type, Study = estudio)
    })
  })
}

################################################################################
# FUNCIÓN PARA INTEGRAR UN TIPO CELULAR ESPECÍFICO
# ENTRE LOS DATASETS
################################################################################

interseccion <- function(tipo_celular) {
  
  rutas_todos <- c(
    "/clinicfs/userhomes/ppalero/GSE157827/GSEA/GO",
    "/clinicfs/userhomes/ppalero/GSE157827/GSEA/Reactome",
    "/clinicfs/userhomes/ppalero/GSE157827/GSEA/Aging",
    "/clinicfs/userhomes/ppalero/SYN18485171/GSEA/GO",
    "/clinicfs/userhomes/ppalero/SYN18485171/GSEA/Reactome",
    "/clinicfs/userhomes/ppalero/SYN18485171/GSEA/Aging",
    "/clinicfs/userhomes/ppalero/SYN52293442/GSEA/GO",
    "/clinicfs/userhomes/ppalero/SYN52293442/GSEA/Reactome",
    "/clinicfs/userhomes/ppalero/SYN52293442/GSEA/Aging"
  )
  
  data_total <- all_info(rutas_todos, tipo_celular)
  bases_datos <- unique(data_total$Database) 
  
  # Iteración interna por cada DB
  for (db in bases_datos) {
    
    data_db <- data_total %>% filter(Database == db)
    
    sig_terms <- data_db %>%
      filter(p.adjust <= 0.05) %>%
      pull(ID) %>%
      unique()
    
    tabla_integrada <- data_db %>%
      filter(ID %in% sig_terms) %>%
      select(
        ID,
        Description,
        setSize,
        enrichmentScore,
        NES,
        pvalue,
        p.adjust,
        qvalue,
        rank,
        leading_edge,
        core_enrichment,
        Condition,
        Study,
        Database,
        Cell_Type,
      ) %>%
      arrange(ID, Study, Condition)
    
    output_dir <- "/clinicfs/userhomes/ppalero/ALL/Integracion"
    if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    file_name <- paste0("Integracion_", db, "_", tipo_celular, ".csv")
    write_csv(tabla_integrada, file.path(output_dir, file_name))
    
    message(paste("Guardada intersección", db, "para:", tipo_celular))
  }
}

################################################################################
# BUCLE PARA GUARDAR TODAS LAS INTERSECCIONES
# POR TIPO CELULAR
################################################################################

tipos <- c("Microglía", "Neuronas", "Astrocitos", "Oligodendrocitos")

for (t in tipos) {
  interseccion(t)
}

################################################################################
# TABLAS INTERSECCIÓN POR BASE DE DATOS
################################################################################

dir_integracion <- "/clinicfs/userhomes/ppalero/ALL/Integracion"

# ----------------------------------- GO ---------------------------------------
cell_types <- c("Microglía", "Astrocitos", "Neuronas", "Oligodendrocitos")

tabla_integrada_GO <- map_dfr(cell_types, function(ct) {
  read_csv(file.path(dir_integracion, paste0("Integracion_GO_", ct, ".csv"))) %>%
    mutate(Cell_Type = ct)
})

write_csv(tabla_integrada_GO, file.path(dir_integracion, "GO_Integracion_Total.csv"))

# --------------------------------- Reactome -----------------------------------

tabla_integrada_React <- map_dfr(cell_types, function(ct) {
  read_csv(file.path(dir_integracion, paste0("Integracion_Reactome_", ct, ".csv"))) %>%
    mutate(Cell_Type = ct)
})

write_csv(tabla_integrada_React, file.path(dir_integracion, "Reactome_Integracion_Total.csv"))

# ---------------------------------- Aging -------------------------------------

tabla_integrada_Aging <- map_dfr(cell_types, function(ct) {
  read_csv(file.path(dir_integracion, paste0("Integracion_Aging_", ct, ".csv"))) %>%
    mutate(Cell_Type = ct)
})

write_csv(tabla_integrada_Aging, file.path(dir_integracion, "Aging_Integracion_Total.csv"))



