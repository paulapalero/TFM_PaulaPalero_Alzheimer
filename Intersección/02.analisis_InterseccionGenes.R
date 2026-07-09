################################################################################
# ANÁLISIS INDIVIDUAL DE LAS RUTAS COMPARTIDAS ENTRE TIPOS CELULARES
# UTILIZANDO LAS INTERSECCIONES DE LAS BASES DE DATOS
# ENTRE GSE157827, SYN18485171 Y SYN52293442
# PARA OBTENER LA INTERSECCIÓN DE GENES 
################################################################################

################################################################################
# LIBRERÍAS
################################################################################

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(readr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)

################################################################################
# CONFIGURACIÓN
################################################################################

base_path <- "/clinicfs/userhomes/ppalero/ALL/Integracion"

rutas_dea <- list(
  "GSE157827"   = "/clinicfs/userhomes/ppalero/GSE157827/DEA/Wilcoxon",
  "SYN18485171" = "/clinicfs/userhomes/ppalero/SYN18485171/DEA/Wilcoxon",
  "SYN52293442" = "/clinicfs/userhomes/ppalero/SYN52293442/DEA/Wilcoxon"
)

comparaciones <- c(
  "GSE157827_AD"        = "AD1",
  "SYN18485171_AD"      = "AD2",
  "SYN52293442_AD"      = "AD3",
  "GSE157827_Control"   = "Control1",
  "SYN18485171_Control" = "Control2",
  "SYN52293442_Control" = "Control3"
)

abreviaturas <- c(
  "Astrocitos"       = "Ast",
  "Neuronas"         = "Neu",
  "Microglía"        = "Mic",
  "Oligodendrocitos" = "Oli"
)

orden_columnas_base <- c(
  "Ast_AD1",  "Ast_AD2",  "Ast_AD3",  "Ast_Control1",  "Ast_Control2",  "Ast_Control3",
  "Neu_AD1",  "Neu_AD2",  "Neu_AD3",  "Neu_Control1",  "Neu_Control2",  "Neu_Control3",
  "Mic_AD1",  "Mic_AD2",  "Mic_AD3",  "Mic_Control1",  "Mic_Control2",  "Mic_Control3",
  "Oli_AD1",  "Oli_AD2",  "Oli_AD3",  "Oli_Control1",  "Oli_Control2",  "Oli_Control3"
)

# Umbral para mostrar nombres de genes
UMBRAL_NOMBRES_GENES <- 40

################################################################################
# FUNCIONES DE INTERÉS POR BASE DE DATOS
################################################################################

funciones_por_bd <- list(
  
  "Aging" = c(
    "GeneAge_UP",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE_UP",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_ALLOGRAFT_REJECTION_UP",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE_UP",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING_UP",
    "CellAge_UP",
    "GeneAge_DOWN"
  ),
  
  "GO" = c(
    "protein refolding",
    "response to unfolded protein",
    "cellular response to heat",
    "chaperone cofactor-dependent protein refolding",
    "response to virus",
    "cholesterol biosynthetic process",
    "positive regulation of cell population proliferation",
    "protein folding",
    "proton transmembrane transport",
    "skeletal muscle cell differentiation"
  ),
  
  "Reactome" = c(
    "Nervous system development",
    "Neuronal System",
    "L1CAM interactions",
    "Axon guidance"
  )
)

archivos_gsea <- list(
  "Aging"    = file.path(base_path, "Aging_Integracion_Total.csv"),
  "GO"       = file.path(base_path, "GO_Integracion_Total.csv"),
  "Reactome" = file.path(base_path, "Reactome_Integracion_Total.csv")
)

################################################################################
# FUNCIÓN: LEER TODOS LOS DEA
################################################################################

leer_dea_completo <- function(rutas_dea, comparaciones, abreviaturas) {
  
  map_dfr(names(rutas_dea), function(estudio) {
    
    archivos <- list.files(rutas_dea[[estudio]], pattern = "\\.csv$", full.names = TRUE)
    
    map_dfr(archivos, function(f) {
      
      nombre      <- basename(f) %>% str_remove("\\.csv$")
      partes      <- str_split(nombre, "_")[[1]]
      cell_type   <- partes[1]
      condition   <- partes[2]
      clave       <- paste(estudio, condition, sep = "_")
      comparacion <- comparaciones[clave]
      
      if (is.na(comparacion)) {
        message("Clave no encontrada: ", clave)
        return(NULL)
      }
      
      read_csv(f, show_col_types = FALSE) %>%
        select(Gene, log2FC, padj) %>%
        mutate(
          Cell_Type   = recode(cell_type, !!!abreviaturas),
          Comparacion = comparacion
        )
    })
  })
}

################################################################################
# FUNCIÓN AUXILIAR: CONSTRUIR MATRICES A PARTIR DEL DEA FILTRADO
################################################################################

construir_matrices <- function(dea_completo, genes_funcion, orden_columnas_base) {
  
  dea_filtrado <- dea_completo %>%
    filter(Gene %in% genes_funcion) %>%
    mutate(Columna = paste(Cell_Type, Comparacion, sep = "_")) %>%
    distinct(Gene, Columna, .keep_all = TRUE)
  
  if (nrow(dea_filtrado) == 0) return(NULL)
  
  mat <- dea_filtrado %>%
    select(Gene, Columna, log2FC) %>%
    pivot_wider(names_from = Columna, values_from = log2FC) %>%
    column_to_rownames("Gene") %>%
    as.matrix()
  
  mat_p <- dea_filtrado %>%
    select(Gene, Columna, padj) %>%
    pivot_wider(names_from = Columna, values_from = padj) %>%
    column_to_rownames("Gene") %>%
    as.matrix()
  
  cols_presentes <- intersect(orden_columnas_base, colnames(mat))
  mat   <- mat[, cols_presentes, drop = FALSE]
  mat_p <- mat_p[, cols_presentes, drop = FALSE]
  
  list(mat = mat, mat_p = mat_p)
}

################################################################################
# FILTRO 1: Significativos en ≥1 casilla
################################################################################

filtro_significativos <- function(mat, mat_p) {
  sig <- apply(mat_p, 1, function(x) any(!is.na(x) & x < 0.05))
  list(mat = mat[sig, , drop = FALSE], mat_p = mat_p[sig, , drop = FALSE])
}

################################################################################
# EXTRAER UN SOLO TIPO CELULAR
################################################################################

extraer_tipo_celular <- function(mat, mat_p, tipo) {
  
  cols <- grep(
    paste0("^", tipo, "_"),
    colnames(mat),
    value = TRUE
  )
  
  list(
    mat   = mat[, cols, drop = FALSE],
    mat_p = mat_p[, cols, drop = FALSE]
  )
}

################################################################################
# FILTRO 2: Tendencias OPUESTAS estricta entre AD y Control en ≥3 tipos celulares
#
# Para cada tipo celular se evalúa si:
#   - Todas las 3 réplicas de AD son positivas (> 0) Y todas las 3 de Control son negativas (< 0)
#   - O BIEN: Todas las de AD son negativas (< 0) Y todas las de Control son positivas (> 0)
################################################################################

filtro_tendencias_opuestas <- function(mat, mat_p,
                                       tipos     = c("Ast", "Neu", "Mic", "Oli"),
                                       min_tipos = 2) {
  
  pasa <- apply(mat, 1, function(gen) {
    
    tipos_ok <- sapply(tipos, function(ct) {
      
      ad_vals  <- gen[grep(paste0("^", ct, "_AD"),      names(gen))]
      ctl_vals <- gen[grep(paste0("^", ct, "_Control"), names(gen))]
      
      # Si hay algún NA en estas casillas, no podemos garantizar el bloque estricto
      if (any(is.na(ad_vals)) || any(is.na(ctl_vals))) return(FALSE)
      
      # Comprobamos homogeneidad de signos por bloque
      ad_pos  <- all(ad_vals > 0)
      ad_neg  <- all(ad_vals < 0)
      ctl_pos <- all(ctl_vals > 0)
      ctl_neg <- all(ctl_vals < 0)
      
      # Patrón opuesto estricto: AD todo (+) y Control todo (-), o viceversa
      (ad_pos && ctl_neg) || (ad_neg && ctl_pos)
    })
    
    sum(tipos_ok) >= min_tipos
  })
  
  list(mat = mat[pasa, , drop = FALSE], mat_p = mat_p[pasa, , drop = FALSE])
}

################################################################################
# FILTRO OPUESTOS EN UN SOLO TIPO CELULAR
################################################################################

filtro_opuestos_tipo <- function(mat, mat_p, tipo) {
  
  pasa <- apply(mat, 1, function(gen) {
    
    ad_vals <- gen[grep("_AD", names(gen))]
    ctl_vals <- gen[grep("_Control", names(gen))]
    
    if (any(is.na(ad_vals)) || any(is.na(ctl_vals))) {
      return(FALSE)
    }
    
    ad_pos  <- all(ad_vals > 0)
    ad_neg  <- all(ad_vals < 0)
    
    ctl_pos <- all(ctl_vals > 0)
    ctl_neg <- all(ctl_vals < 0)
    
    (ad_pos && ctl_neg) ||
      (ad_neg && ctl_pos)
  })
  
  list(
    mat   = mat[pasa, , drop = FALSE],
    mat_p = mat_p[pasa, , drop = FALSE]
  )
}

################################################################################
# FILTRO 3: Tendencias IGUALES estricta entre AD y Control en ≥3 tipos celulares
# Mismo criterio pero las 6 casillas del tipo celular deben compartir el mismo signo.
################################################################################

filtro_tendencias_iguales <- function(mat, mat_p,
                                      tipos     = c("Ast", "Neu", "Mic", "Oli"),
                                      min_tipos = 2) {
  
  pasa <- apply(mat, 1, function(gen) {
    
    tipos_ok <- sapply(tipos, function(ct) {
      
      ad_vals  <- gen[grep(paste0("^", ct, "_AD"),      names(gen))]
      ctl_vals <- gen[grep(paste0("^", ct, "_Control"), names(gen))]
      
      # Si hay algún NA en estas casillas, no podemos garantizar el bloque estricto
      if (any(is.na(ad_vals)) || any(is.na(ctl_vals))) return(FALSE)
      
      # Comprobamos homogeneidad de signos por bloque
      ad_pos  <- all(ad_vals > 0)
      ad_neg  <- all(ad_vals < 0)
      ctl_pos <- all(ctl_vals > 0)
      ctl_neg <- all(ctl_vals < 0)
      
      # Patrón igual estricto: Las 6 casillas son (+) o las 6 casillas son (-)
      (ad_pos && ctl_pos) || (ad_neg && ctl_neg)
    })
    
    sum(tipos_ok) >= min_tipos
  })
  
  list(mat = mat[pasa, , drop = FALSE], mat_p = mat_p[pasa, , drop = FALSE])
}

################################################################################
# FILTRO IGUALES EN UN SOLO TIPO CELULAR
################################################################################

filtro_iguales_tipo <- function(mat, mat_p, tipo) {
  
  pasa <- apply(mat, 1, function(gen) {
    
    ad_vals <- gen[grep("_AD", names(gen))]
    ctl_vals <- gen[grep("_Control", names(gen))]
    
    if (any(is.na(ad_vals)) || any(is.na(ctl_vals))) {
      return(FALSE)
    }
    
    ad_pos  <- all(ad_vals > 0)
    ad_neg  <- all(ad_vals < 0)
    
    ctl_pos <- all(ctl_vals > 0)
    ctl_neg <- all(ctl_vals < 0)
    
    (ad_pos && ctl_pos) ||
      (ad_neg && ctl_neg)
  })
  
  list(
    mat   = mat[pasa, , drop = FALSE],
    mat_p = mat_p[pasa, , drop = FALSE]
  )
}

################################################################################
# FUNCIÓN CENTRAL: GRAFICAR Y GUARDAR EL HEATMAP
# Usa log2FC continuo (gradiente) con asteriscos para padj < 0.05.
# Muestra nombres de genes solo si nrow(mat) < UMBRAL_NOMBRES_GENES.
################################################################################

dibujar_heatmap <- function(mat, mat_p, titulo, nombre_archivo, retornar_objeto = FALSE) {
  
  # Eliminar filas donde TODAS las columnas son NA 
  filas_ok <- apply(mat, 1, function(x) sum(!is.na(x)) >= 2)
  mat   <- mat[filas_ok,   , drop = FALSE]
  mat_p <- mat_p[filas_ok, , drop = FALSE]
  
  n_genes <- nrow(mat)
  
  if (n_genes == 0) {
    message("Sin genes para: ", titulo)
    return(NULL) # Cambiado a NULL para poder controlar si está vacío
  }
  
  # Reemplazar NA por 0 solo para el clustering
  mat_clust           <- mat
  mat_clust[is.na(mat_clust)] <- 0
  
  # Configurar Escala de Color Plano por Signo
  lim <- max(abs(mat), na.rm = TRUE)
  lim <- ifelse(is.finite(lim) && lim > 0, lim, 1)
  
  col_fun <- colorRamp2(
    c(-lim,  -0.00001,  0,       0.00001,  lim),
    c("#2166AC", "#2166AC", "white", "#B2182B", "#B2182B")
  )
  
  bloques <- factor(
    str_extract(colnames(mat), "^[^_]+"),
    levels = c("Ast", "Neu", "Mic", "Oli")
  )
  bloques <- droplevels(bloques)
  
  mostrar_genes <- n_genes < UMBRAL_NOMBRES_GENES
  mat_p_local <- mat_p
  
  ht <- Heatmap(
    mat_clust,
    name            = "log2FC",
    col             = col_fun,
    cluster_rows    = TRUE,
    cluster_columns = FALSE,
    
    show_row_names  = mostrar_genes,
    row_names_gp    = gpar(fontsize = 7),
    
    column_title    = titulo,
    column_title_gp = gpar(fontsize = 11, fontface = "bold"), 
    
    column_split    = bloques,
    top_annotation  = HeatmapAnnotation(
      Cell = anno_block(
        gp        = gpar(fill = "grey95", col = "white"),
        labels    = levels(bloques),
        labels_gp = gpar(fontsize = 9, fontface = "bold")
      )
    ),
    
    rect_gp         = gpar(col = "white", lwd = 0.8),
    column_gap      = unit(2, "mm"),
    column_names_gp = gpar(fontsize = 7),
    border          = TRUE,
    
    heatmap_legend_param = list(
      title          = "log2FC",
      color_bar      = "continuous",
      legend_height  = unit(2, "cm"),
      at             = c(-lim, 0, lim),
      labels         = c("Down", "0 / NA", "Up")
    ),
    
    layer_fun = function(j, i, x, y, width, height, fill) {
      v <- pindex(mat_p_local, i, j)
      sig <- !is.na(v) & v < 0.05
      if (any(sig)) {
        grid.text("*", x[sig], y[sig],
                  gp = gpar(fontsize = 10, fontface = "bold", col = "black"))
      }
    }
  )
  
  # ── AQUÍ CONTROLAMOS SI SE GUARDA O SE RETORNA ─────────────────────────────
  if (retornar_objeto) {
    return(ht)
  }
  
  # Código original de guardado para los gráficos generales
  alto_px <- max(800, n_genes * (if (mostrar_genes) 55 else 22))
  png(nombre_archivo, width = 3200, height = alto_px, res = 300)
  draw(ht)
  dev.off()
  
  cat("  Guardado:", nombre_archivo, "(", n_genes, "genes )\n")
}

################################################################################
# FUNCIÓN UNIFICADORA: GENERA LOS 3 HEATMAPS POR FUNCIÓN
#   1. General   — todos los genes significativos en ≥1 casilla
#   2. Opuestos  — tendencias AD vs Control opuestas en ≥2 tipos celulares
#   3. Iguales   — tendencias AD vs Control iguales  en ≥2 tipos celulares
################################################################################

generar_heatmaps_funcion <- function(funcion_interes, genes_funcion, dea_completo,
                                     orden_columnas_base, output_dir2) {
  
  cat("\n--- Procesando:", funcion_interes,
      "(", length(genes_funcion), "genes en la función ) ---\n")
  
  # ── Construir matrices brutas ────────────────────────────────────────────────
  mats <- construir_matrices(dea_completo, genes_funcion, orden_columnas_base)
  if (is.null(mats)) {
    message("Sin datos DEA para: ", funcion_interes)
    return(invisible(NULL))
  }
  
  # ── Heatmap 1: general (solo filtro significatividad) ───────────────────────
  f1 <- filtro_significativos(mats$mat, mats$mat_p)
  cat("  [General]  genes sig:", nrow(f1$mat), "\n")
  dibujar_heatmap(
    mat            = f1$mat,
    mat_p          = f1$mat_p,
    titulo         = paste0("Función: ", funcion_interes),
    nombre_archivo = file.path(output_dir2,
                               paste0(make.names(funcion_interes), "_general.png"))
  )
  
  # ── Heatmap 2: tendencias opuestas ──────────────────────────────────────────
  f2 <- filtro_tendencias_opuestas(f1$mat, f1$mat_p)
  cat("  [Opuestos] genes:", nrow(f2$mat), "\n")
  dibujar_heatmap(
    mat            = f2$mat,
    mat_p          = f2$mat_p,
    titulo         = paste0("Función: ", funcion_interes,
                            "\n(tendencias AD/Control opuestas, ≥2 tipos celulares)"),
    nombre_archivo = file.path(output_dir2,
                               paste0(make.names(funcion_interes), "_opuestos.png"))
  )
  
  # ── Heatmap 3: tendencias iguales ───────────────────────────────────────────
  f3 <- filtro_tendencias_iguales(f1$mat, f1$mat_p)
  cat("  [Iguales]  genes:", nrow(f3$mat), "\n")
  dibujar_heatmap(
    mat            = f3$mat,
    mat_p          = f3$mat_p,
    titulo         = paste0("Función: ", funcion_interes,
                            "\n(tendencias AD/Control iguales, ≥2 tipos celulares)"),
    nombre_archivo = file.path(output_dir2,
                               paste0(make.names(funcion_interes), "_iguales.png"))
  )
  
  ###############################################################################
  # HEATMAPS POR TIPO CELULAR - CONSOLIDADO DIRECTO EN MEMORIA
  ###############################################################################
  
  tipos <- c("Ast", "Neu", "Mic", "Oli")
  
  # Listas para guardar los objetos en memoria RAM
  lista_ht_opuestos <- list()
  lista_ht_iguales  <- list()
  
  # Variables para calcular dinámicamente el alto del lienzo final
  max_genes_opuestos <- 0
  max_genes_iguales  <- 0
  
  for (tipo in tipos) {
    
    cel <- extraer_tipo_celular(f1$mat, f1$mat_p, tipo)
    
    # --------------------------------------------------------------------------
    # PROCESAR OPUESTOS 
    # --------------------------------------------------------------------------
    op <- filtro_opuestos_tipo(cel$mat, cel$mat_p, tipo)
    n_genes_op <- nrow(op$mat)
    cat("    Tipo:", tipo, "- Opuestos:", n_genes_op, "genes\n")
    
    ht_op <- dibujar_heatmap(
      mat              = op$mat,
      mat_p            = op$mat_p,
      titulo           = paste0(funcion_interes, " - ", tipo, " (OPUESTOS)"),
      nombre_archivo   = NULL,
      retornar_objeto  = TRUE
    )
    
    if (!is.null(ht_op)) {
      lista_ht_opuestos[[tipo]] <- ht_op
      if (n_genes_op > max_genes_opuestos) max_genes_opuestos <- n_genes_op
    }
    
    # --------------------------------------------------------------------------
    # PROCESAR IGUALES 
    # --------------------------------------------------------------------------
    ig <- filtro_iguales_tipo(cel$mat, cel$mat_p, tipo)
    n_genes_ig <- nrow(ig$mat)
    cat("    Tipo:", tipo, "- Iguales:", n_genes_ig, "genes\n")
    
    ht_ig <- dibujar_heatmap(
      mat              = ig$mat,
      mat_p            = ig$mat_p,
      titulo           = paste0(funcion_interes, " - ", tipo, " (IGUALES)"),
      nombre_archivo   = NULL,
      retornar_objeto  = TRUE
    )
    
    if (!is.null(ht_ig)) {
      lista_ht_iguales[[tipo]] <- ht_ig
      if (n_genes_ig > max_genes_iguales) max_genes_iguales <- n_genes_ig
    }
  }
  
  # --------------------------------------------------------------------------
  # GUARDAR DIRECTAMENTE LOS PANELES CONSOLIDADOS 
  # --------------------------------------------------------------------------
  
  # Dibujar Panel de Opuestos si contiene gráficos
  if (length(lista_ht_opuestos) > 0) {
    archivo_opuestos <- file.path(output_dir2, paste0(make.names(funcion_interes), "_opuestos_CONSOLIDADO.png"))
    
    # Calculamos el alto basándonos en el tipo celular que más genes tenga
    mostrar_genes_op <- max_genes_opuestos < UMBRAL_NOMBRES_GENES
    alto_panel_op <- max(1200, max_genes_opuestos * (if (mostrar_genes_op) 55 else 22) * 2)
    
    png(archivo_opuestos, width = 3600, height = alto_panel_op, res = 300)
    
    # Creamos la cuadrícula de 2x2
    grid.newpage()
    lay <- grid.layout(2, 2)
    pushViewport(viewport(layout = lay))
    
    # Dibujar cada cuadrante usando pushViewport + draw + popViewport
    if ("Ast" %in% names(lista_ht_opuestos)) {
      pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
      draw(lista_ht_opuestos[["Ast"]], newpage = FALSE)
      popViewport()
    }
    if ("Neu" %in% names(lista_ht_opuestos)) {
      pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
      draw(lista_ht_opuestos[["Neu"]], newpage = FALSE)
      popViewport()
    }
    if ("Mic" %in% names(lista_ht_opuestos)) {
      pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
      draw(lista_ht_opuestos[["Mic"]], newpage = FALSE)
      popViewport()
    }
    if ("Oli" %in% names(lista_ht_opuestos)) {
      pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
      draw(lista_ht_opuestos[["Oli"]], newpage = FALSE)
      popViewport()
    }
    
    popViewport() # Cierra el layout general
    dev.off()
    cat("  ✓ Guardado Panel Consolidado Opuestos:", basename(archivo_opuestos), "\n")
  }
  
  # Dibujar Panel de Iguales si contiene gráficos
  if (length(lista_ht_iguales) > 0) {
    archivo_iguales <- file.path(output_dir2, paste0(make.names(funcion_interes), "_iguales_CONSOLIDADO.png"))
    
    mostrar_genes_ig <- max_genes_iguales < UMBRAL_NOMBRES_GENES
    alto_panel_ig <- max(1200, max_genes_iguales * (if (mostrar_genes_ig) 55 else 22) * 2)
    
    png(archivo_iguales, width = 3600, height = alto_panel_ig, res = 300)
    
    grid.newpage()
    lay <- grid.layout(2, 2)
    pushViewport(viewport(layout = lay))
    
    if ("Ast" %in% names(lista_ht_iguales)) {
      pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
      draw(lista_ht_iguales[["Ast"]], newpage = FALSE)
      popViewport()
    }
    if ("Neu" %in% names(lista_ht_iguales)) {
      pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
      draw(lista_ht_iguales[["Neu"]], newpage = FALSE)
      popViewport()
    }
    if ("Mic" %in% names(lista_ht_iguales)) {
      pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
      draw(lista_ht_iguales[["Mic"]], newpage = FALSE)
      popViewport()
    }
    if ("Oli" %in% names(lista_ht_iguales)) {
      pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
      draw(lista_ht_iguales[["Oli"]], newpage = FALSE)
      popViewport()
    }
    
    popViewport() # Cierra el layout general
    dev.off()
    cat("  ✓ Guardado Panel Consolidado Iguales:", basename(archivo_iguales), "\n")
  } 
  
} 


################################################################################
# PIPELINE PRINCIPAL
################################################################################

cat("Leyendo todos los archivos DEA...\n")
dea_completo <- leer_dea_completo(rutas_dea, comparaciones, abreviaturas)
cat("Total de filas en DEA completo:", nrow(dea_completo), "\n")

for (bd in names(funciones_por_bd)) {
  
  cat("\n========================================\n")
  cat("Base de datos:", bd, "\n")
  cat("========================================\n")
  
  if (!file.exists(archivos_gsea[[bd]])) {
    warning("Archivo GSEA no encontrado: ", archivos_gsea[[bd]])
    next
  }
  
  tabla_gsea <- read_csv(archivos_gsea[[bd]], show_col_types = FALSE)
  output_dir2 <- file.path("Heatmaps_Individules_de_Funciones", bd)
  dir.create(output_dir2, showWarnings = FALSE, recursive = TRUE)
  
  for (funcion_interes in funciones_por_bd[[bd]]) {
    
    genes_funcion <- tabla_gsea %>%
      filter(Description == funcion_interes) %>%
      pull(core_enrichment) %>%
      str_split("/") %>%
      unlist() %>%
      unique()
    
    if (length(genes_funcion) == 0) {
      message("Función no encontrada en tabla GSEA: ", funcion_interes)
      next
    }
    
    generar_heatmaps_funcion(
      funcion_interes     = funcion_interes,
      genes_funcion       = genes_funcion,
      dea_completo        = dea_completo,
      orden_columnas_base = orden_columnas_base,
      output_dir          = output_dir2
    )
  }
}

cat("\n Análisis terminado.\n")
