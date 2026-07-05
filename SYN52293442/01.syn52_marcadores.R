########################################################
# ANALISIS DE ANOTACION
# CORRESPONDENCIA DE MARCADORES CON TIPOS CELULARES
########################################################

########################################################
# LIBRERIAS
########################################################

library(Seurat)
library(ggplot2)
library(patchwork)

########################################################
# CARGAR DATOS 
########################################################

setwd("/clinicfs/userhomes/ppalero/SYN52293442")
syn52_raw <- readRDS("syn52_raw.rds")

########################################################
# SEURAT PROCESS
########################################################

seurat <- NormalizeData(syn52_raw)
seurat <- FindVariableFeatures(so1)
seurat <- ScaleData(so1, features = VariableFeatures(object = so1))
seurat <- RunPCA(so1, features = VariableFeatures(object = so1))
seurat <- FindNeighbors(so1) # por defecto
seurat <- FindClusters(so1, resolution = 0.8)
seurat <- RunUMAP(so1, dims = 1:10)


saveRDS(seurat, file = "syn52_processed.rds")

########################################################
# UMAP INDIVIDUALES DE MARCADORES
# ESPECÍFICOS POR TIPO CELULAR
########################################################

# Lista de marcadores
markers_list <- list(
  Neuronas = c("SNAP25","SYT1","RBFOX3","TUBB3"),
  Excitadoras = c("SLC17A7"),
  Inhibidoras = c("GAD1","GAD2"),
  Astrocitos = c("AQP4","GFAP","ALDH1L1","SLC1A2"),
  Microglia = c("C1QA","C1QB","CX3CR1","TMEM119"),
  Oligos = c("MBP","PLP1","MOG","MOBP"),
  OPCs = c("PDGFRA","CSPG4","NG2")
)

pdf("Marcadores.pdf", width = 10, height = 8)

p_umap <- DimPlot(seurat, group.by = "major_cell_type", 
                  label = TRUE) +
  ggtitle("UMAP - Tipos celulares") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold")
  )

print(p_umap)

for (celltype in names(markers_list)) {
  
  genes <- markers_list[[celltype]]
  
  # Obtener plots individuales
  plots <- FeaturePlot(seurat, features = genes, combine = FALSE)
  
  # Añadir título a cada marcador
  plots <- lapply(1:length(plots), function(i) {
    plots[[i]] + ggtitle(genes[i]) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
      )
  })
  
  # Combinar y añadir título (tipo celular)
  p <- patchwork::wrap_plots(plots, ncol = 2) +
    plot_annotation(title = celltype)
  
  print(p)
}

dev.off()

