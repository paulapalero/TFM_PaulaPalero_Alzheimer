################################################################################
# SCRIPT AUTOMATIZADO PARA REALIZAR UN 
# ANÁLISIS DE EXPRESIÓN DIFERENCIAL
# SOBRE LOS DATOS DEL ESTUDIO GSE157827
# EMPLEANDO VARIOS MÉTODOS
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
library(ggplot2)
library(patchwork)

################################################################################
# CARGAR DATOS (SCE PREPROCESADO)
################################################################################

setwd("/clinicfs/userhomes/ppalero/GSE157827")
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
meta$Diagnosis  <- meta$CONDITION
meta$Sex        <- meta$SEX
meta$Age        <- meta$AGE

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

gene_symbols <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                      keys = rownames(seurat),
                                      column = "SYMBOL",
                                      keytype = "ENSEMBL",
                                      multiVals = "first")

gene_symbols[is.na(gene_symbols)] <- names(gene_symbols)[is.na(gene_symbols)]

# Crear nueva columna en Metadata
seurat@meta.data$GeneSymbol <- "all_cells"

# Sumar filas repetidas (Colapsar ENSEMBL -> SYMBOL)
# Usamos una matriz de diseño dispersa
group_factor <- factor(gene_symbols)

# Creamos una matriz que dice qué ENSEMBL pertenece a qué SYMBOL
dot_matrix <- sparse.model.matrix(~ 0 + group_factor)
colnames(dot_matrix) <- levels(group_factor)

# Sumar las filas que comparten el mismo símbolo
counts_collapsed <- t(dot_matrix) %*% GetAssayData(seurat,
                                                   assay = "originalexp",
                                                   layer = "counts")
rownames(counts_collapsed) <- levels(group_factor)

meta_data_backup <- seurat@meta.data

seurat <- CreateSeuratObject(
  counts = counts_collapsed,
  meta.data = meta_data_backup
)


################################################################################
## SEURAT PROCESS
################################################################################

seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat)
seurat <- RunUMAP(seurat, dims = 1:30)
saveRDS(seurat, file = "seurat_gse15.rds")


################################################################################
# VARIABLES
# CELLTYPES: Extraer los distintos tipos celulares
# COMPARISONS: lista con las comparaciones
################################################################################

celltypes <- unique(seurat$CellType) 

comparisons <- list(
  
  Control_Late_vs_Early = list(
    subset="NC",
    ident1="Late",
    ident2="Early"
  ),
  
  AD_Late_vs_Early = list(
    subset="AD",
    ident1="Late",
    ident2="Early"
  )
  
)

methods <- c("MAST","Wilcoxon","Pseudobulk")

################################################################################
# CREAR CARPETAS PARA LOS RESULTADOS
################################################################################

dir.create("DEA2", showWarnings=FALSE)

for(m in methods){
  dir.create(file.path("DEA2",m), showWarnings=FALSE)
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

########################################################
# FUNCION PSEUDOBULK
# Hace el Análisis de Expresión Diferencial con pb
########################################################

run_pseudobulk <- function(seurat_obj){
  
  counts <- GetAssayData(seurat_obj, layer="counts")
  
  meta <- seurat_obj@meta.data
  
  group <- paste(meta$Sample_ID, meta$Age_Group, sep="_")
  
  pb <- rowsum(t(counts), group)
  pb <- t(pb)
  
  meta_pb <- data.frame(sample=colnames(pb))
  
  meta_pb$Age_Group <- sapply(strsplit(meta_pb$sample,"_"), `[`, 2)
  
  dds <- DESeqDataSetFromMatrix(
    countData=pb,
    colData=meta_pb,
    design=~Age_Group
  )
  counts <- GetAssayData(seurat, layer = "counts")
  rownames(counts) <- sub("\\..*", "", rownames(counts))
  dds <- DESeq(dds)
  
  res <- results(dds, independentFiltering = FALSE, contrast=c("Age_Group","Late","Early")) 
  
  return(as.data.frame(res))
  
}

################################################################################
# DATAFRAMES PARA GUARDAR TIEMPOS Y GENES SIGNIFICATIVOS
################################################################################

tiempos_ejecucion <- data.frame(
  Comparison = character(),
  CellType   = character(),
  Method     = character(),
  Time_sec   = numeric(),
  stringsAsFactors = FALSE
)

genes_significativos <- data.frame(
  Comparison = character(),
  CellType   = character(),
  Method     = character(),
  Num_Genes  = numeric(),
  stringsAsFactors = FALSE
)

################################################################################
# LOOP PRINCIPAL
################################################################################

for(comp_name in names(comparisons)){
  
  cat("\n=============================\n")
  cat("COMPARACIÓN:", comp_name,"\n")
  cat("=============================\n")
  
  comp <- comparisons[[comp_name]]
  seu <- seurat
  
  if(!is.null(comp$subset)){
    seu <- subset(seu, Diagnosis==comp$subset)
  }
  
  Idents(seu) <- "Age_Group"
  
  for(ct in celltypes){
    
    cat("\nCellType:", ct,"\n")
    cells_ct <- WhichCells(seu, expression=CellType==ct)
    
    if(length(cells_ct) < 50){
      cat("Saltando",ct,"(muy pocas células)\n")
      next
    }
    
    seu_ct <- subset(seu, cells=cells_ct)
    
    #######################
    # WILCOXON
    #######################
    
    cat("Ejecutando Wilcoxon...\n")
    inicio <- Sys.time()
    
    wilcox <- FindMarkers(
      seu_ct,
      ident.1=comp$ident1,
      ident.2=comp$ident2,
      test.use="wilcox",
      logfc.threshold=0
    )
    
    wilcox_fmt <- format_results(wilcox,"Wilcoxon",ct,comp_name)
    write.csv(wilcox_fmt, file=file.path("DEA2","Wilcoxon", paste0(ct,"_",comp_name,".csv")), 
              row.names=FALSE)
    
    fin <- Sys.time()
    
    # Tiempos
    tiempos_ejecucion <- rbind(tiempos_ejecucion, data.frame(
      Comparison = comp_name, CellType = ct, Method = "Wilcoxon",
      Time_sec = as.numeric(difftime(fin, inicio, units = "secs"))
    ))
    # Genes (padj < 0.05)
    genes_significativos <- rbind(genes_significativos, data.frame(
      Comparison = comp_name, CellType = ct, Method = "Wilcoxon",
      Num_Genes = sum(wilcox_fmt$padj < 0.05, na.rm = TRUE)
    ))
    
    #######################
    # MAST
    #######################
    
    cat("Ejecutando MAST...\n")
    inicio <- Sys.time()
    
    mast <- FindMarkers(
      seu_ct,
      ident.1=comp$ident1,
      ident.2=comp$ident2,
      test.use="MAST",
      logfc.threshold=0
    )
    
    mast_fmt <- format_results(mast,"MAST",ct,comp_name) 
    write.csv(mast_fmt, file=file.path("DEA2","MAST", paste0(ct,"_",comp_name,".csv")), 
              row.names=FALSE)
    
    fin <- Sys.time()
    
    # Tiempos
    tiempos_ejecucion <- rbind(tiempos_ejecucion, data.frame(
      Comparison = comp_name, CellType = ct, Method = "MAST",
      Time_sec = as.numeric(difftime(fin, inicio, units = "secs"))
    ))
    # Genes
    genes_significativos <- rbind(genes_significativos, data.frame(
      Comparison = comp_name, CellType = ct, Method = "MAST",
      Num_Genes = sum(mast_fmt$padj < 0.05, na.rm = TRUE)
    ))
    
    #######################
    # PSEUDOBULK
    #######################
    
    cat("Ejecutando Pseudobulk...\n")
    inicio <- Sys.time()
    
    pb <- run_pseudobulk(seu_ct)
    pb_fmt <- format_results(pb,"Pseudobulk",ct,comp_name)
    write.csv(pb_fmt, file=file.path("DEA2","Pseudobulk", paste0(ct,"_",comp_name,".csv")), row.names=FALSE)
    
    fin <- Sys.time()
    
    # Tiempos
    tiempos_ejecucion <- rbind(tiempos_ejecucion, data.frame(
      Comparison = comp_name, CellType = ct, Method = "Pseudobulk",
      Time_sec = as.numeric(difftime(fin, inicio, units = "secs"))
    ))
    # Genes
    genes_significativos <- rbind(genes_significativos, data.frame(
      Comparison = comp_name, CellType = ct, Method = "Pseudobulk",
      Num_Genes = sum(pb_fmt$padj < 0.05, na.rm = TRUE)
    ))
  }
}

# Guardar los csv resultantes
write.csv(tiempos_ejecucion, "DEA2/Tiempos_Ejecucion.csv", row.names = FALSE)
write.csv(genes_significativos, "DEA2/Genes_Significativos.csv", row.names = FALSE)

################################################################################
# GUARDAR RESULTADOS Y GRAFICAR
################################################################################

# Nombres limpios (AD y Control)
tiempos_ejecucion$Grupo <- ifelse(grepl("^AD", tiempos_ejecucion$Comparison), "AD", "Control")
genes_significativos$Grupo <- ifelse(grepl("^AD", genes_significativos$Comparison), "AD", "Control")

# Paleta de colores 
colores_metodos <- c("Pseudobulk" = "#4ec3a3", "MAST" = "#ff8c5a", "Wilcoxon" = "#8fa1d4")

################################################################################
# GRÁFICO A: TIEMPOS DE EJECUCIÓN
################################################################################

plot_A <- ggplot(tiempos_ejecucion, aes(x = CellType, y = Time_sec + 1, fill = Method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~Grupo) +
  scale_y_log10(labels = scales::comma) + # Escala logarítmica 
  scale_fill_manual(values = colores_metodos) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"), 
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    legend.position = "none" 
  ) +
  labs(
    title = "A  Tiempo de Ejecución por Método",
    x = "Tipo Celular",
    y = "Tiempo (Segundos, escala: log(x+1))"
  )

################################################################################
# GRÁFICO B: GENES SIGNIFICATIVOS DETECTADOS
################################################################################

plot_B <- ggplot(genes_significativos, aes(x = CellType, y = Num_Genes + 1, fill = Method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~Grupo) +
  scale_y_log10(labels = scales::comma) + 
  scale_fill_manual(values = colores_metodos) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"), # 
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  ) +
  labs(
    title = "B  Genes Significativos Detectados",
    x = "Tipo Celular",
    y = "Nº de Genes (escala: log(x+1))",
    fill = "Método"
  )

################################################################################
# JUNTAR AMBOS GRÁFICOS 
################################################################################

grafico_final <- plot_A / plot_B

# Guardar en alta calidad
ggsave("DEA2/Grafico_Tiempos_y_Genes.png", plot = grafico_final, width = 11, 
       height = 7, dpi = 300)
