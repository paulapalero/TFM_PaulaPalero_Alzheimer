################################################################################
# SCRIPT AUTOMATIZADO PARA
# CREAR TABLA RESUMEN DE GENES SIGNIFICATIVOS (UP/DOWN)
# SOBRE LOS RESULTADOS DEA DE SYN18485171
################################################################################

################################################################################
# LIBRERIAS
################################################################################

library(dplyr)
library(tidyr)

################################################################################
# DIRECTORIO RESULTADOS
################################################################################

setwd("/clinicfs/userhomes/ppalero/SYN18485171/")
results_dir <- "DEA"

methods <- list.dirs(results_dir, recursive = FALSE, full.names = FALSE)

################################################################################
# FUNCION CREAR TABLA RESUMEN
################################################################################

create_summary <- function(method){
  
  cat("\n=================================\n")
  cat("Procesando Método:", method, "\n")
  cat("=================================\n")
  
  method_dir <- file.path(results_dir, method)
  
  files <- list.files(method_dir, pattern="csv$", full.names=TRUE)
  
  if(length(files)==0){
    cat("No se han encontrado los archivos\n")
    summarise(n_genes = n_distinct(Gene)) %>%
      arrange(cell_group, desc(n_genes))
    return(NULL)
  }
  
  ####################################################
  # LEER RESULTADOS
  ####################################################
  
  res_list <- lapply(files, read.csv)
  
  all_res <- bind_rows(res_list)
  
  ####################################################
  # CLASIFICAR GENES
  ####################################################
  
  res_sig <- all_res %>%
    filter(!is.na(padj)) %>%
    mutate(direction = case_when(
      padj < 0.05 & log2FC > 0  ~ "UP",
      padj < 0.05 & log2FC < 0 ~ "DOWN",
      TRUE ~ "NS"
    )) %>%
    filter(direction != "NS")
  
  ####################################################
  # CONTAR GENES
  ####################################################
  
  counts <- res_sig %>%
    group_by(CellType, Comparison, direction) %>%
    summarise(n = n_distinct(Gene), .groups="drop")
  
  ####################################################
  # CREAR TABLA MATRIZ
  ####################################################
  
  summary_table <- counts %>%
    tidyr::unite("Comparison_direction", Comparison, direction) %>%
    pivot_wider(
      names_from = Comparison_direction,
      values_from = n,
      values_fill = 0
    )
  
  ####################################################
  # EXPORTAR
  ####################################################
  
  write.csv(
    summary_table,
    paste0("syn18_DEA_summary.csv"),
    row.names=FALSE
  )
  
  cat("Resumen guardado:", paste0("syn18_DEA_summary.csv"), "\n")
  
}

################################################################################
# LOOP SOBRE MÉTODO
################################################################################

for(m in methods){
  create_summary(m)
}

cat("\nLa tabla ha sido creada correctamente\n")
