########################################################
# SCRIPT AUTOMATIZADO PARA REALIZAR UN 
# ANÁLISIS DE EXPRESIÓN DIFERENCIAL
# SOBRE LOS DATOS DEL ESTUDIO SYN52293442
########################################################

########################################################
# LIBRERIAS: Cargamos los paquetes necesarios
########################################################

library(Seurat)
library(dplyr)
library(tidyverse)
library(MAST)
library(DESeq2)
library(Matrix)
library(tibble)

########################################################
# CARGAR DATOS 
########################################################

setwd("/clinicfs/userhomes/ppalero/SYN52293442")
seurat <- readRDS("syn52_processed.rds")

########################################################
# HOMOGENEIZAR METADATOS
########################################################

meta <- seurat@meta.data

meta$Sample_ID  <- meta$projid
meta$Diagnosis  <- meta$condition
meta$Sex        <- meta$sex
meta$Age        <- as.integer(as.numeric(gsub("\\+", "", meta$age_death.y)))

meta$Age_Group <- ifelse(meta$Age > 85, "Late", "Early")
meta$Age_Group <- factor(meta$Age_Group,
                         levels=c("Late","Early"))
names <- c(
  "Ast" = "Astrocitos",
  "Mic" = "Microglía",
  "Exc" = "Neuronas",    # Agrupamos Exc -> Neuronas
  "Inh" = "Neuronas",    # Agrupamos Inh -> Neuronas
  "Oli" = "Oligodendrocitos",
  "OPC" = "Oligodendrocitos" # Agrupamos OPC -> Oligo
)

meta$CellType <- names[as.character(meta$major_cell_type)]
seurat@meta.data <- meta

seurat <- subset(seurat, subset = CellType %in% c("Astrocitos", 
                                                  "Microglía", 
                                                  "Neuronas", 
                                                  "Oligodendrocitos"))

########################################################
## SEURAT PROCESS
########################################################

# Cluster por tipos celulares
png("syn52_clusters.png", width = 2000, height = 1000, res = 300)
print(DimPlot(seurat, 
              group.by = c("major_cell_type"), 
              shuffle =  FALSE,
              pt.size = 0.1,
              raster = FALSE))
dev.off()

# Grafico: Número de células por paciente
# Contar número de células por paciente y tipo celular
cell_counts <- seurat@meta.data %>%
  group_by(Sample_ID, Diagnosis, CellType) %>%
  summarise(n_cells = n(), .groups = "drop")

# Calcular proporción de cada tipo celular por paciente
cell_props <- cell_counts %>%
  group_by(Sample_ID) %>%
  mutate(prop = n_cells / sum(n_cells)) %>%
  ungroup()


# Graficar proporciones con barras
png("syn52_cellProp.png", width = 2400, height = 1200, res = 300)

print(
  # El único cambio real está aquí: as.character(Sample_ID)
  ggplot(cell_props, aes(x = as.character(Sample_ID), y = prop, fill = CellType)) +
    geom_col(width = 0.7) +                      
    facet_wrap(~ Diagnosis, scales = "free_x") + 
    theme_minimal(base_size = 12) +              
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.major.x = element_blank(),     
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16)
    ) +
    labs(title = "Proporción de tipos celulares por paciente",
         x = "Paciente (ID)",
         y = "Proporción de células",
         fill = "Tipo Celular")
)

dev.off()

saveRDS(seurat, file = "seurat_syn52.rds")

########################################################
## VARIABLES
# CELLTYPES: Extraer los distintos tipos celulares
# COMPARISONS: lista con las comparaciones
########################################################

celltypes <- unique(seurat$CellType) 

comparisons <- list(
  
  Control_Late_vs_Early = list(
    subset="Control",
    ident1="Late",
    ident2="Early"
  ),
  
  AD_Late_vs_Early = list(
    subset="AD",
    ident1="Late",
    ident2="Early"
  )
  
)

methods <- c("Wilcoxon")

########################################################
# CREAR CARPETAS PARA LOS RESULTADOS
########################################################

dir.create("DEA", showWarnings=FALSE)

for(m in methods){
  dir.create(file.path("DEA",m), showWarnings=FALSE)
}

########################################################
# FUNCION PARA ESTANDARIZAR LOS RESULTADOS
########################################################

format_results <- function(df, method, celltype, comparison){
  
  df <- df %>%
    rownames_to_column("Gene")
  
  if("avg_log2FC" %in% colnames(df)){
    df$log2FC <- df$avg_log2FC
  }
  
  if("log2FoldChange" %in% colnames(df)){
    df$log2FC <- df$log2FoldChange
  }
  
  if("p_val_adj" %in% colnames(df)){
    df$padj <- df$p_val_adj
  }
  
  if("pvalue" %in% colnames(df)){
    df$pval <- df$pvalue
    
  }
  
  df <- df %>%
    dplyr::select(Gene, log2FC, padj)
  
  df$Method <- method
  df$CellType <- celltype
  df$Comparison <- comparison
  
  return(df)
  
}


########################################################
# LOOP PRINCIPAL
########################################################

for(comp_name in names(comparisons)){ # Itera sobre las 3 comparaciones
  
  comp <- comparisons[[comp_name]]
  
  seu <- seurat
  
  if(!is.null(comp$subset)){
    seu <- subset(seu, Diagnosis==comp$subset)
  }
  
  Idents(seu) <- "Age_Group"
  
  for(ct in celltypes){ # Itera sobre cell.type dentro de cada comparación
    
    cells_ct <- WhichCells(seu, expression = CellType==ct)
    
    if(length(cells_ct) < 50) next
    
    seu_ct <- subset(seu, cells=cells_ct)
    
    #######################
    # WILCOXON
    #######################
    
    cat("Running Wilcoxon...\n")
    
    wilcox <- FindMarkers(
      seu_ct,
      ident.1=comp$ident1,
      ident.2=comp$ident2,
      test.use="wilcox",
      logfc.threshold=0
    )
    
    wilcox <- format_results(wilcox,"Wilcoxon",ct,comp_name)
    
    write.csv(
      wilcox,
      file=file.path("DEA","Wilcoxon",
                     paste0(ct,"_",comp_name,".csv")),
      row.names=FALSE
    )
    
    cat("Finished Wilcoxon\n")
    
  }
  
}
