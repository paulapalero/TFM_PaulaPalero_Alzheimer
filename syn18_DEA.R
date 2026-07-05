################################################################################
# SCRIPT AUTOMATIZADO PARA REALIZAR UN 
# ANÁLISIS DE EXPRESIÓN DIFERENCIAL
# SOBRE LOS DATOS DEL ESTUDIO SYN18485171
################################################################################

################################################################################
# LIBRERIAS
################################################################################

library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(tidyverse)
library(MAST)
library(DESeq2)
library(Matrix)
library(biomaRt)
library(tibble)
library(AnnotationDbi)
library(org.Hs.eg.db)

################################################################################
# CARGAR DATOS
################################################################################

setwd("/clinicfs/userhomes/ppalero/SYN18485171")
sce <- readRDS("sce.annot.clean.rds")

################################################################################
# CONVERTIR A SEURAT
################################################################################

seurat <- as.Seurat(sce, counts = "counts")

################################################################################
# HOMOGENEIZAR METADATOS
################################################################################

meta <- seurat@meta.data

meta$Sample_ID  <- meta$PatientID
meta$Diagnosis  <- meta$condition
meta$Sex        <- meta$sex
meta$Age        <- as.numeric(gsub("\\+", "", seurat$age_death))

meta$CellType   <- meta$cell.type
meta$CellType[meta$CellType %in% c("Excitadoras", "Inhibidoras")] <- "Neuronas"
meta$CellType[meta$CellType %in% c("OPC")] <- "Oligodendrocitos"

meta$Age_Group <- ifelse(meta$Age > 85, "Late", "Early")
meta$Age_Group <- factor(meta$Age_Group,
                         levels=c("Late","Early"))

seurat@meta.data <- meta

################################################################################
# CONVERSION - ENSEMBL TO GENE_SYMBOL
################################################################################

# Extraer la matriz de conteos
counts <- GetAssayData(seurat, layer = "counts")

# Mapear ENSEMBL a Gene Symbol
gene_symbols <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                      keys = rownames(counts),
                                      column = "SYMBOL",
                                      keytype = "ENSEMBL",
                                      multiVals = "first")

# Manejar los NAs (genes que no se pudieron mapear)
gene_symbols[is.na(gene_symbols)] <- names(gene_symbols)[is.na(gene_symbols)]

counts <- rowsum(counts, group = gene_symbols)
counts<- counts[order(rownames(counts)), ]

seurat <- CreateSeuratObject(
  counts = counts,
  meta.data = seurat@meta.data
)

################################################################################
## SEURAT PROCESS
################################################################################

seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat)
seurat <- RunUMAP(seurat, dims = 1:30)
saveRDS(seurat, file = "seurat_syn18.rds")

# --------------------- Gráfico: Clusters --------------------------------------

# Definir la paleta de colores
colores_clusters <- brewer.pal(4, "Set2") 
names(colores_clusters) <- c("Astrocitos", 
                             "Microglía", 
                             "Neuronas", 
                             "Oligodendrocitos")

# Extraer las coordenadas UMAP
# Calcular el centroide (la media) de cada clúster
coordenadas <- as.data.frame(Embeddings(seurat_syn18, reduction = "umap"))
coordenadas$CellType <- seurat_syn18$CellType

etiquetas_pos <- coordenadas %>%
  group_by(CellType) %>%
  summarize(umap_1 = mean(umap_1), umap_2 = mean(umap_2))

png("syn18_clusters.png", width = 2400, height = 1800, res = 300)

DimPlot(seurat_syn18,
        group.by = "CellType",
        shuffle = TRUE,
        pt.size = 0.4,
        label = FALSE, 
        raster = FALSE) + 
  scale_color_manual(values = colores_clusters) +
  
  geom_label(data = etiquetas_pos, 
             aes(x = umap_1, y = umap_2, label = CellType, fill = CellType), 
             color = "black", 
             fontface = "plain", 
             size = 4.5, 
             alpha = 0.75, 
             label.padding = unit(0.3, "lines")) +
  
  scale_fill_manual(values = colores_clusters) +
  
  labs(title = "SYN18485171 CellTypes", 
       x = "UMAP 1", 
       y = "UMAP 2") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "none"
  )

dev.off()

# ----------------- Gráfico: Número de células por paciente --------------------

# Contar número de células por paciente y tipo celular
cell_counts <- seurat@meta.data %>%
  group_by(Sample_ID, Diagnosis, CellType) %>%
  summarise(n_cells = n(), .groups = "drop")

# Calcular proporción de cada tipo celular por paciente
cell_props <- cell_counts %>%
  group_by(Sample_ID) %>%
  mutate(prop = n_cells / sum(n_cells)) %>%
  ungroup()

png("syn18_cellProp.png", width = 2000, height = 1000, res = 300)
print(
ggplot(cell_props, aes(x = Sample_ID, y = prop, fill = CellType)) +
  geom_col() +
  facet_wrap(~ Diagnosis, scales = "free_x") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Proporción de tipos celulares por paciente",
       x = "Paciente (ID)",
       y = "Proporción de células",
       fill = "Tipo Celular")
)
dev.off()


################################################################################
# VARIABLES
# CELLTYPES: Extraer los distintos tipos celulares
# COMPARISONS: lista con las comparaciones
################################################################################

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

################################################################################
# CREAR CARPETAS PARA LOS RESULTADOS
################################################################################

dir.create("DEA", showWarnings=FALSE)

for(m in methods){
  dir.create(file.path("DEA",m), showWarnings=FALSE)
}

################################################################################
# FUNCION PARA ESTANDARIZAR LOS RESULTADOS
################################################################################

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


################################################################################
# LOOP PRINCIPAL
################################################################################

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
    
    cat("Ejecutando Wilcoxon...\n")
    
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
    
    cat("Análisis Terminado\n")
    
  }
  
}
