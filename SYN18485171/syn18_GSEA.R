################################################################################
# SCRIPT AUTOMATIZADO PARA 
# REALIZAR UN ANÁLISIS FUNCIONAL
# SOBRE LOS DATOS DE SYN18485171
################################################################################

################################################################################
# LIBRERIAS
################################################################################

library(clusterProfiler)
library(biomaRt)
library(org.Hs.eg.db)
library(GO.db)
library(AnnotationDbi)
library(reactome.db)
library(dplyr)
library(tibble)

################################################################################
# CONFIGURAR DIRECTORIOS
################################################################################

input_dir <- "/clinicfs/userhomes/ppalero/SYN18485171/DEA/Wilcoxon"
base_out  <- "/clinicfs/userhomes/ppalero/SYN18485171"
genesets_path <- "/clinicfs/userhomes/ppalero/genesets"

# Definir y crear directorios para GO
annot_dir_go  <- file.path(base_out, "annotations/GO")
output_gsea_go <- file.path(base_out, "GSEA/GO")

# Definir y crear directorios para Reactome
annot_dir_reac  <- file.path(base_out, "annotations/Reactome")
output_gsea_reac <- file.path(base_out, "GSEA/Reactome")

# Definir y crear directorios para Aging
annot_dir_aging   <- file.path(base_out, "annotations/Aging")
output_gsea_aging <- file.path(base_out, "GSEA/Aging")

# Crear todos los directorios
dirs <- c(annot_dir_go, output_gsea_go, 
          annot_dir_reac, output_gsea_reac,
          annot_dir_aging, output_gsea_aging)
lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

################################################################################
# PREPARACIÓN DE LAS ANOTACIONES 
################################################################################

message(">> Preparando bases de datos de anotación...")

# ------------------ Preparar GO (Biological Process) --------------------------

message("   -> Procesando Gene Ontology...")

keys_go <- keys(org.Hs.eg.db, keytype="SYMBOL")
human_go <- AnnotationDbi::select(org.Hs.eg.db, 
                                  keys = keys_go,
                                  columns = c("SYMBOL", "GO", "ONTOLOGY"),
                                  keytype = "SYMBOL")
go_terms <- AnnotationDbi::select(GO.db, 
                                  keys = unique(human_go$GO), 
                                  columns = "TERM", 
                                  keytype = "GOID")
h_bp <- merge(human_go, go_terms, by.x = "GO", by.y = "GOID") %>%
  filter(ONTOLOGY == "BP") %>%
  na.omit() %>%
  distinct()

T2G_GO <- h_bp[, c("GO", "SYMBOL")]
colnames(T2G_GO) <- c("Term", "Gene")
T2N_GO <- h_bp[, c("GO", "TERM")]
colnames(T2N_GO) <- c("Term", "Name")

# ---------------------- Preparar Reactome -------------------------------------

message("   -> Procesando Reactome...")

mart <- useMart(biomart = "ensembl", 
                dataset = "hsapiens_gene_ensembl", 
                host = "https://www.ensembl.org")
h_reac <- getBM(attributes=c('external_gene_name', 'reactome'), 
                values=TRUE, mart=mart)

h_reac <- h_reac %>%
  filter(reactome != "" & external_gene_name != "") %>%
  distinct()

T2G_REAC <- h_reac[, c("reactome", "external_gene_name")]
colnames(T2G_REAC) <- c("Term", "Gene")

react_ids <- unique(h_reac$reactome)
react_terms <- mapIds(reactome.db, 
                      keys = react_ids, 
                      column = "PATHNAME", 
                      keytype = "PATHID", 
                      multiVals = "first")

T2N_REAC <- data.frame(
  Term = names(react_terms),
  Name = gsub("^Homo sapiens: ", "", as.character(react_terms)),
  stringsAsFactors = FALSE
)

# --------------------------- Preparar AGING -----------------------------------

message("   -> Procesando Gene Sets de Aging...")

# 1. MSigDB 
DEMA_up <- read.gmt(file.path(genesets_path, "DEMAGALHAES_AGING_UP.v2026.1.Hs.gmt"))
DEMA_dn <- read.gmt(file.path(genesets_path, "DEMAGALHAES_AGING_DN.v2026.1.Hs.gmt"))
colnames(DEMA_up) <- colnames(DEMA_dn) <- c("Term", "Gene")

ALLOGRAFT_up <- read.gmt(file.path(genesets_path,"HALLMARK_ALLOGRAFT_REJECTION.v2026.1.Hs.gmt")) %>%
  mutate(term = "HALLMARK_ALLOGRAFT_REJECTION_UP") %>%
  rename(Term = term, Gene = gene)

IL6_up <- read.gmt(file.path(genesets_path,"HALLMARK_IL6_JAK_STAT3_SIGNALING.v2026.1.Hs.gmt")) %>%
  mutate(term = "HALLMARK_IL6_JAK_STAT3_SIGNALING_UP") %>%
  rename(Term = term, Gene = gene)

INFLAMM <- read.gmt(file.path(genesets_path, "HALLMARK_INFLAMMATORY_RESPONSE.v2026.1.Hs.gmt"))
colnames(INFLAMM) <- c("Term", "Gene")

interAlpha_up <- read.gmt(file.path(genesets_path,"HALLMARK_INTERFERON_ALPHA_RESPONSE.v2026.1.Hs.gmt")) %>%
  mutate(term = "HALLMARK_INTERFERON_ALPHA_RESPONSE_UP") %>%
  rename(Term = term, Gene = gene)

interGamma_up <- read.gmt(file.path(genesets_path,"HALLMARK_INTERFERON_GAMMA_RESPONSE.v2026.1.Hs.gmt")) %>%
  mutate(term = "HALLMARK_INTERFERON_GAMMA_RESPONSE_UP") %>%
  rename(Term = term, Gene = gene)

TNFA <- read.gmt(file.path(genesets_path, "HALLMARK_TNFA_SIGNALING_VIA_NFKB.v2026.1.Hs.gmt"))
colnames(TNFA) <- c("Term", "Gene")

GOBP <- read.gmt(file.path(genesets_path, "GOBP_CELLULAR_SENESCENCE.v2026.1.Hs.gmt"))
colnames(GOBP) <- c("Term", "Gene")

HP <- read.gmt(file.path(genesets_path, "HP_NEUROINFLAMMATION.v2026.1.Hs.gmt"))
colnames(HP) <- c("Term", "Gene")


# 2. CellAge 
cellage <- read.csv(file.path(genesets_path, "cellAge_signatures2020.csv"), sep = ";")
ca_up <- data.frame(Term = "CellAge_UP", 
                    Gene = cellage$gene_symbol[cellage$ovevrexp == 1], 
                    stringsAsFactors = FALSE)
ca_dn <- data.frame(Term = "CellAge_DOWN", 
                    Gene = cellage$gene_symbol[cellage$underexp == 1], 
                    stringsAsFactors = FALSE)

# 3. GenAge 
ga_up_raw <- read.csv(file.path(genesets_path, "geneage_up.csv"), sep = ";")
ga_up <- data.frame(Term = "GeneAge_UP", 
                    Gene = ga_up_raw$Gene, 
                    stringsAsFactors = FALSE)

ga_dn_raw <- read.csv(file.path(genesets_path, "geneage_down.csv"), sep = ";")
ga_dn <- data.frame(Term = "GeneAge_DOWN", 
                    Gene = ga_dn_raw$Gene, 
                    stringsAsFactors = FALSE)

# 4. EnrichR
enrich_up_raw <- read.csv(file.path(genesets_path, "enrichR_up.csv"), sep = ",")
enrich_up <- data.frame(Term = "Enrich_UP", 
                        Gene = enrich_up_raw$Gene, 
                        stringsAsFactors = FALSE)

enrich_dn_raw <- read.csv(file.path(genesets_path, "enrichR_down.csv"), sep = ",")
enrich_dn <- data.frame(Term = "Enrich_DOWN", 
                        Gene = enrich_dn_raw$Gene, 
                        stringsAsFactors = FALSE)


# -------------- Consolidación AGING -------------------
T2G_AGING <- rbind(DEMA_up, DEMA_dn,
                   ALLOGRAFT_up,
                   IL6_up,
                   INFLAMM,
                   TNFA,
                   GOBP,
                   HP,
                   interAlpha_up,
                   interGamma_up,
                   ca_up, ca_dn, 
                   ga_up, ga_dn, 
                   enrich_up, enrich_dn)

T2G_AGING <- T2G_AGING %>% 
  filter(!is.na(Gene) & Gene != "") %>% 
  distinct()

T2N_AGING <- data.frame(Term = unique(T2G_AGING$Term), 
                        Name = unique(T2G_AGING$Term), 
                        stringsAsFactors = FALSE)

################################################################################
# BUCLE PRINCIPAL
################################################################################

files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

for (file_path in files) {
  
  file_name <- basename(file_path)
  base_id <- gsub(".csv", "", file_name)
  message(paste("\n********** PROCESANDO:", base_id, "**********"))
  
  diffexp <- read.csv(file_path, sep = ",") 
  diffexp <- diffexp[!is.na(diffexp$log2FC) & !is.na(diffexp$Gene), ]
  
  geneList <- diffexp$log2FC
  names(geneList) <- diffexp$Gene
  geneList <- sort(geneList, decreasing = TRUE)
  
  # ---------- EJECUTAR GSEA GO -------------
  message(">> Ejecutando GSEA GO (BP)...")
  res_go <- GSEA(geneList = geneList, 
                 TERM2GENE = T2G_GO, 
                 TERM2NAME = T2N_GO,
                 nPermSimple = 10000,
                 minGSSize = 10, 
                 maxGSSize = 500, 
                 pvalueCutoff = 1, 
                 verbose = FALSE)
  
  # --------- EJECUTAR GSEA REACTOME --------
  message(">> Ejecutando GSEA Reactome...")
  res_reac <- GSEA(geneList = geneList, 
                   TERM2GENE = T2G_REAC, 
                   TERM2NAME = T2N_REAC,
                   nPermSimple = 10000, 
                   minGSSize = 10, 
                   maxGSSize = 500, 
                   pvalueCutoff = 1, 
                   verbose = FALSE)
  
  # ---------- EJECUTAR GSEA AGING ---------
  message(">> Ejecutando GSEA Aging...")
  res_aging <- GSEA(geneList = geneList, 
                    TERM2GENE = T2G_AGING, 
                    TERM2NAME = T2N_AGING,
                    nPermSimple = 10000, 
                    minGSSize = 5,  
                    maxGSSize = 500, 
                    pvalueCutoff = 1, 
                    verbose = FALSE)
  
  # --------------------- GUARDAR RESULTADOS GO --------------------------------
  write.table(diffexp, file = file.path(annot_dir_go, paste0("annot_", 
                                                             base_id, ".txt")), 
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.csv(as.data.frame(res_go), file = file.path(output_gsea_go, 
                                                    paste0("gseaGO_", 
                                                           base_id, ".csv")), 
            row.names = FALSE)
  saveRDS(res_go, file = file.path(output_gsea_go, paste0("gseaGO_", 
                                                          base_id, ".rds")))
  
  # -------------------- GUARDAR RESULTADOS REACTOME ---------------------------
  write.table(diffexp, file = file.path(annot_dir_reac, paste0("annot_", 
                                                               base_id, ".txt")), 
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.csv(as.data.frame(res_reac), file = file.path(output_gsea_reac, 
                                                      paste0("gseaReac_", 
                                                             base_id, ".csv")), 
            row.names = FALSE)
  saveRDS(res_reac, file = file.path(output_gsea_reac, paste0("gseaReac_", 
                                                              base_id, ".rds")))
  
  # --------------------- GUARDAR RESULTADOS AGING -----------------------------
  write.table(diffexp, file = file.path(annot_dir_aging, paste0("annot_", 
                                                                base_id, ".txt")), 
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.csv(as.data.frame(res_aging), 
            file = file.path(output_gsea_aging, paste0("gseaAging_", base_id, ".csv")), 
            row.names = FALSE)
  
  saveRDS(res_aging, 
          file = file.path(output_gsea_aging, paste0("gseaAging_", base_id, ".rds")))
  
  message(paste(">> Finalizado con éxito:", base_id))
}

message("\n>> Proceso completado.")