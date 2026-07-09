################################################################################
# ANÁLISIS DE LAS RUTAS COMPARTIDAS ENTRE TIPOS CELULARES
# UTILIZANDO LAS INTERSECCIONES DE LAS BASES DE DATOS
# ENTRE GSE157827, SYN18485171 Y SYN52293442
################################################################################

################################################################################
# LIBRERÍAS
################################################################################

library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(tidytext)
library(stringr)

################################################################################
# DIRECTORIOS
################################################################################

base_path <- "/clinicfs/userhomes/ppalero/ALL/Integracion"

output_path <- file.path(base_path,"Rutas_Compartidas")

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

################################################################################
# CARGAR TABLAS 
################################################################################

tabla_GO <- read_csv(file.path(base_path, "GO_Integracion_Total.csv"))

tabla_Reactome <- read_csv(file.path(base_path, "Reactome_Integracion_Total.csv"))

tabla_Aging <- read_csv(file.path(base_path, "Aging_Integracion_Total.csv"))

################################################################################
# CONSTRUIR DATASET DE GENES
################################################################################

procesar_db <- function(df, nombre_db){
  
  df %>%
    select(
      Cell_Type,
      Description,
      core_enrichment
    ) %>%
    mutate(
      Gene = strsplit(core_enrichment, "/")
    ) %>%
    unnest(Gene) %>%
    distinct(
      Cell_Type,
      Description,
      Gene
    ) %>%
    mutate(
      DB = nombre_db
    )
  
}

genes_GO <- procesar_db(tabla_GO, "GO")

genes_Reactome <- procesar_db(tabla_Reactome,"Reactome")

genes_Aging <- procesar_db(tabla_Aging,"Aging")

genes_todos <- bind_rows(
  genes_GO,
  genes_Reactome,
  genes_Aging
)

# Abreviaturas
genes_todos <- genes_todos %>%
  mutate(
    Cell_Type_short = recode(
      Cell_Type,
      "Astrocitos"       = "Ast",
      "Neuronas"         = "Neu",
      "Microglía"        = "Mic",
      "Oligodendrocitos" = "Oli"
    )
  )

################################################################################
# COMBINACIONES CELULARES
################################################################################

tipos_unicos <- sort(unique(genes_todos$Cell_Type_short))

combinaciones_todas <- c(
  combn(tipos_unicos, 2, simplify = FALSE),
  combn(tipos_unicos, 3, simplify = FALSE),
  combn(tipos_unicos, 4, simplify = FALSE)
)

################################################################################
# INTERSECCIÓN DE GENES
################################################################################

rutas_intersecciones <- map_dfr(
  combinaciones_todas,
  function(combo){
    
    nombre_combo <- paste(
      combo,
      collapse = " + "
    )
    
    datos_combo <- genes_todos %>%
      filter(
        Cell_Type_short %in% combo
      )
    
    rutas_validas <- datos_combo %>%
      group_by(DB, Description) %>%
      summarise(
        n_tipos = n_distinct(Cell_Type_short),
        .groups = "drop"
      ) %>%
      filter(
        n_tipos == length(combo)
      )
    
    resultado <- map_dfr(
      seq_len(nrow(rutas_validas)),
      function(i){
        
        ruta_actual <- rutas_validas$Description[i]
        db_actual   <- rutas_validas$DB[i]
        
        df_ruta <- datos_combo %>%
          filter(
            Description == ruta_actual,
            DB == db_actual
          )
        
        genes_por_tipo <- split(
          df_ruta$Gene,
          df_ruta$Cell_Type_short
        )
        
        genes_comunes <- Reduce(
          intersect,
          genes_por_tipo
        )
        
        tibble(
          DB = db_actual,
          Description = ruta_actual,
          tipos = nombre_combo,
          n_genes_compartidos = length(
            genes_comunes
          ),
          genes_compartidos = paste(
            sort(genes_comunes),
            collapse = " | "
          )
        )
        
      }
    )
    
    resultado
    
  }
)

################################################################################
# ELIMINAR RUTAS SIN GENES COMPARTIDOS
################################################################################

rutas_intersecciones <- rutas_intersecciones %>%
  filter(n_genes_compartidos > 0)

################################################################################
# EXPORTAR TABLA
################################################################################

write_csv(rutas_intersecciones, file.path(output_path,
                                          "RutasCompartidas_Celltypes.csv"))

################################################################################
# TOP 15 POR COMBINACIÓN
################################################################################

top_rutas_compartidas <- rutas_intersecciones %>%
  group_by(tipos) %>%
  slice_max(
    n_genes_compartidos,
    n = 15,
    with_ties = FALSE
  ) %>%
  ungroup()

################################################################################
# GRÁFICO 
################################################################################

p <- ggplot(
  top_rutas_compartidas,
  aes(
    x = reorder_within(
      Description,
      n_genes_compartidos,
      tipos
    ),
    y = n_genes_compartidos,
    fill = DB
  )
) +
  geom_col() +
  coord_flip() +
  facet_wrap(
    ~ tipos,
    scales = "free_y",
    ncol = 3
  ) +
  scale_x_reordered() +
  labs(
    title = "Rutas compartidas entre tipos celulares",
    x = NULL,
    y = "Genes compartidos"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(
      face = "bold",
      size = 10
    ),
    legend.position = "bottom",
    axis.text.y = element_text(
      size = 7
    )
  )

ggsave(
  file.path(
    output_path,
    "RutasCompartidas_Celltypes.png"
  ),
  p,
  width = 16,
  height = 10,
  dpi = 300
)

################################################################################
# PDF
################################################################################

pdf(
  file.path(
    output_path,
    "RutasCompartidas_Celltypes.pdf"
  ),
  width = 10,
  height = 8
)

for(combo in unique(top_rutas_compartidas$tipos)){
  
  datos_plot <- top_rutas_compartidas %>%
    filter(
      tipos == combo
    )
  
  p <- ggplot(
    datos_plot,
    aes(
      x = reorder(
        Description,
        n_genes_compartidos
      ),
      y = n_genes_compartidos,
      fill = DB
    )
  ) +
    geom_col() +
    coord_flip() +
    labs(
      title = paste("Rutas compartidas entre los tipos celulares:",combo),
      x = NULL,
      y = "Genes compartidos"
    ) +
    theme_minimal(base_size = 12)
  
  print(p)
  
}

dev.off()



