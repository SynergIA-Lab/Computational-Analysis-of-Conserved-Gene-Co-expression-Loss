# ============================================================
# TCGA RNA-seq (STAR - Counts) - ADQUISICIÓN Y CURACIÓN
# Versión: Limpia (Sin Pureza Tumoral / ESTIMATE)
# ============================================================

# --- Paquetes ---
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

for (pkg in c("TCGAbiolinks","SummarizedExperiment","DESeq2", "dplyr", "readr", "tibble", "matrixStats", "stringr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(DESeq2)
  library(dplyr)
  library(readr)
  library(tibble)
  library(matrixStats)
  library(stringr)
})

# ============ CONFIGURACIÓN ============
setwd("/Users/mriosc/Documents/NETWORKS")

projects <- c("TCGA-BRCA")  
out_dir  <- "TCGA_RNAseq"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sample_tumor   <- "Primary Tumor"
sample_normal  <- "Solid Tissue Normal"
data_category  <- "Transcriptome Profiling"
experimental_strategy <- "RNA-Seq"
data_type      <- "Gene Expression Quantification"
workflow_used  <- "STAR - Counts"

# Requisitos mínimos
require_both <- TRUE
min_tumors   <- 20
min_normals  <- 10
min_primary_samples_per_tissue <- 20 

# Deduplicación
do_dedup <- TRUE
dedup_priority <- "library_size" 

set.seed(1234)

# ============================================================
# FUNCIONES: DEDUPLICACIÓN
# ============================================================

.patient_id_from_tcga <- function(x) substr(as.character(x), 1, 12)

dedup_one_per_patient <- function(df, counts_mat = NULL, sample_type_col = "sample_type") {
  id_col <- if ("barcode_se" %in% names(df)) "barcode_se" else "barcode"
  
  df <- df %>%
    mutate(patient_id = .patient_id_from_tcga(.data[[id_col]])) %>%
    filter(!is.na(patient_id))
  
  if (is.null(counts_mat)) {
    return(df %>% group_by(patient_id, .data[[sample_type_col]]) %>% slice(1) %>% ungroup())
  }
  
  lib_sizes <- colSums(counts_mat[, intersect(df[[id_col]], colnames(counts_mat)), drop = FALSE])
  lib_size_df <- tibble(sample_id = names(lib_sizes), library_size = as.numeric(lib_sizes))
  
  df_dedup <- df %>%
    left_join(lib_size_df, by = setNames("sample_id", id_col)) %>%
    arrange(desc(library_size)) %>%
    group_by(patient_id, .data[[sample_type_col]]) %>%
    slice(1) %>%
    ungroup()
  
  return(df_dedup)
}

# ============================================================
# FUNCIÓN PRINCIPAL: DESCARGA
# ============================================================

.assay_counts <- function(se) {
  an <- assayNames(se)
  if ("unstranded" %in% an) assay(se, "unstranded") else assay(se)
}

download_project <- function(project_id) {
  prj_dir <- file.path(out_dir, project_id)
  dir.create(prj_dir, showWarnings = FALSE, recursive = TRUE)
  
  message(">> [", project_id, "] Iniciando descarga...")
  
  q <- GDCquery(
    project = project_id,
    data.category = data_category,
    data.type = data_type,
    workflow.type = workflow_used,
    experimental.strategy = experimental_strategy
  )
  
  GDCdownload(q, method="api", files.per.chunk=20)
  se <- GDCprepare(q)
  
  counts_mat <- .assay_counts(se)
  counts_mat <- counts_mat[!grepl("^N_", rownames(counts_mat)), ]
  
  meta_df <- as.data.frame(colData(se))
  meta_df$barcode_se <- colnames(se)
  
  if (do_dedup) {
    meta_df <- dedup_one_per_patient(meta_df, counts_mat)
    counts_mat <- counts_mat[, intersect(meta_df$barcode_se, colnames(counts_mat))]
  }
  
  # Exportación base
  write_tsv(as.data.frame(counts_mat) %>% rownames_to_column("gene_id"), 
            file.path(prj_dir, paste0(project_id, "_counts.tsv")))
  write_tsv(meta_df, file.path(prj_dir, paste0(project_id, "_metadata.tsv")))
  write_tsv(as.data.frame(rowData(se)), file.path(prj_dir, paste0(project_id, "_gene_info.tsv")))
  
  return(list(prj_dir = prj_dir))
}

# ============ EJECUCIÓN DESCARGA ============
for (p in projects) { info <- download_project(p) }

# ============================================================
# FASE: ORGANIZACIÓN POR SUBTIPOS (PAM50 PARA BRCA)
# ============================================================
message("\n>> Iniciando partición por subtipos PAM50 para BRCA...")

p <- "TCGA-BRCA"
prj_dir <- file.path(out_dir, p)

metadata  <- read_tsv(file.path(prj_dir, paste0(p, "_metadata.tsv")), show_col_types = FALSE)
counts    <- read_tsv(file.path(prj_dir, paste0(p, "_counts.tsv")), show_col_types = FALSE) %>% column_to_rownames("gene_id")
gene_info <- read_tsv(file.path(prj_dir, paste0(p, "_gene_info.tsv")), show_col_types = FALSE)

subtype_col <- "paper_BRCA_Subtype_PAM50"

metadata_filt <- metadata %>%
  filter(
    sample_type %in% c(sample_tumor, sample_normal),
    (prior_malignancy != "YES" | is.na(prior_malignancy)),
    (prior_treatment != "YES" | is.na(prior_treatment))
  )

subsets_base_dir <- file.path(prj_dir, "by_Subtype")
dir.create(subsets_base_dir, showWarnings = FALSE, recursive = TRUE)

controls_meta <- metadata_filt %>% filter(sample_type == sample_normal)
tumors_meta   <- metadata_filt %>% filter(sample_type == sample_tumor, !is.na(.data[[subtype_col]]))

valid_subtypes <- tumors_meta %>%
  group_by(.data[[subtype_col]]) %>%
  summarise(n = n()) %>%
  filter(n >= min_primary_samples_per_tissue) %>%
  pull(.data[[subtype_col]])

for (sub in valid_subtypes) {
  message("--- Procesando subtipo: ", sub, " ---")
  
  sub_meta <- bind_rows(controls_meta, tumors_meta %>% filter(.data[[subtype_col]] == sub))
  common_ids <- intersect(sub_meta$barcode_se, colnames(counts))
  sub_meta <- sub_meta %>% 
    filter(barcode_se %in% common_ids) %>%
    mutate(base_label = ifelse(sample_type == sample_tumor, "Primary", "Control")) %>%
    group_by(base_label) %>%
    mutate(new_name = paste0(base_label, row_number())) %>%
    ungroup()
  
  # Mapeo y Colapso a Gene Symbol
  sub_counts <- counts[, sub_meta$barcode_se, drop = FALSE]
  mat_mapped <- as.data.frame(sub_counts) %>%
    rownames_to_column("gene_id") %>%
    inner_join(gene_info %>% select(gene_id, gene_name), by = "gene_id") %>%
    select(-gene_id) %>%
    group_by(gene_name) %>%
    summarise(across(everything(), sum)) %>%
    filter(!is.na(gene_name), gene_name != "")
  
  colnames(mat_mapped) <- c("gene", sub_meta$new_name[match(colnames(mat_mapped)[-1], sub_meta$barcode_se)])
  
  sub_dir <- file.path(subsets_base_dir, gsub("[ /,]", "_", sub))
  dir.create(sub_dir, showWarnings = FALSE, recursive = TRUE)
  
  write_csv(mat_mapped, file.path(sub_dir, "counts_renombrado.csv"))
  write_csv(mat_mapped %>% select(gene, starts_with("Primary")), file.path(sub_dir, "counts_primary.csv"))
  write_csv(mat_mapped %>% select(gene, starts_with("Control")), file.path(sub_dir, "counts_control.csv"))
  write_csv(sub_meta %>% select(barcode_se, sample_type, new_name, !!sym(subtype_col)), file.path(sub_dir, "mapeo_muestras.csv"))
}

message("\n✅ PROCESO FINALIZADO.")
