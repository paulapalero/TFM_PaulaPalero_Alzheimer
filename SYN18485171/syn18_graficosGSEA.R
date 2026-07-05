################################################################################
# SCRIPT AUTOMATIZADO PARA 
# VISUALIZAR LOS RESULTADOS DEL GSEA
# PARA LOS DATOS DE SYN18485171
################################################################################

################################################################################
# LIBRERIAS
################################################################################

library(tidyverse)
library(ggplot2)
library(stringr)

################################################################################
# CONFIGURACIÓN DE DIRECTORIOS
################################################################################

dirs <- list(
  GO = "/clinicfs/userhomes/ppalero/SYN18485171/GSEA/GO",
  Reactome = "/clinicfs/userhomes/ppalero/SYN18485171/GSEA/Reactome",
  Aging = "/clinicfs/userhomes/ppalero/SYN18485171/GSEA/Aging"
)

################################################################################
# FUNCIÓN PARA PROCESAR LOS RESULTADOS GSEA
# Ordenaremos primero por p.adjust (ascendente) y luego por el valor 
# absoluto de NES (descendente). Así, si hay empate en significancia 
# (p.adjust = 1.00), el script elegirá las rutas con mayor impacto biológico.
################################################################################

# Función principal de processamiento
# Y creación de los gráficos

generate_plots <- function(directory, analysis_name) {
  
  output_dir <- file.path(directory, "Plots")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  files <- list.files(directory, pattern = "\\.csv$", full.names = TRUE)
  cell_types <- files %>% 
    basename() %>% 
    str_extract("(?<=_)[^_]+(?=_(AD|Control))") %>% 
    unique() %>% 
    na.omit()
  
  for (cell in cell_types) {
    message(paste("Procesando:", analysis_name, "-", cell))
    
    file_ad <- files[str_detect(files, paste0(cell, "_AD"))]
    file_ctrl <- files[str_detect(files, paste0(cell, "_Control"))]
    
    if (length(file_ad) == 0 | length(file_ctrl) == 0) next
    
    df_ad <- read.csv(file_ad) %>% 
      mutate(Group = paste0(cell, "\n(AD)"))
    
    df_ctrl <- read.csv(file_ctrl) %>% 
      mutate(Group = paste0(cell, "\n(Control)"))
    
    # Ordenamos por p.adjust (menor a mayor) 
    # Y desempatamos por el valor absoluto de NES (mayor a menor)
    top_pathways_ad <- df_ad %>% 
      arrange(p.adjust, desc(abs(NES))) %>% 
      slice_head(n = 15) %>% 
      pull(Description)
    
    top_pathways_ctrl <- df_ctrl %>% 
      arrange(p.adjust, desc(abs(NES))) %>% 
      slice_head(n = 15) %>% 
      pull(Description)
    
    selected_pathways <- unique(c(top_pathways_ad, top_pathways_ctrl))
    
    plot_data <- bind_rows(df_ad, df_ctrl) %>%
      filter(Description %in% selected_pathways) %>%
      mutate(
        log_p = -log10(p.adjust),
        # Si el p.adj es 1, el log_p será 0. Evitamos problemas visuales.
        log_p = ifelse(is.infinite(log_p), 0, log_p),
        Description = str_wrap(Description, width = 60),
        is_sig = ifelse(p.adjust < 0.05, "Significant", "Non-significant")
      )
    
    p <- ggplot(plot_data, aes(x = Group, y = Description)) +
      geom_point(aes(size = log_p, fill = NES, shape = is_sig), 
                 color = "black", stroke = 0.8) +
      # Triángulo con borde (24) y Círculo con borde (21)
      scale_shape_manual(values = c("Significant" = 24, "Non-significant" = 21),
                         name = "Significance") +
      scale_fill_gradient2(low = "dodgerblue4", mid = "white", 
                           high = "firebrick3", midpoint = 0) +
      scale_size_continuous(range = c(2, 9), name = "-log10(p.adj)") +
      theme_bw() +
      labs(title = paste("Análisis Comparativo GSEA", analysis_name),
           subtitle = paste("Late vs Early Aging:", cell),
           x = NULL, y = "Pathways",
           caption = "Triángulo: p.adj < 0.05 | Círculo: No significativo") +
      theme(
        axis.text.x = element_text(face = "bold", color = "black", size = 11),
        axis.text.y = element_text(size = 9, color = "black"),
        panel.grid.major.y = element_line(color = "gray95")
      )
    
    output_file <- file.path(output_dir, paste0("Dotplot_", analysis_name, "_", 
                                                cell, ".png"))
    ggsave(output_file, p, width = 10, height = 10, dpi = 300, bg = "white")
  }
}

# Executar para cada base de datos
generate_plots(dirs$GO, "GO")
generate_plots(dirs$Reactome, "Reactome")
generate_plots(dirs$Aging, "Aging")
