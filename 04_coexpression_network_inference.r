# Load the WGCNA library
library(WGCNA)
library(data.table)   # para fread
library(utils)        # txtProgressBar, setTxtProgressBar (normalmente ya cargado)


################################################################################
# CONFIGURACIÓN DE RUTAS
################################################################################

# Directorio base
BASE_DIR <- "/Users/mriosc/Documents/NETWORKS/TCGA-BRCA/by_Subtype/LumB"

# Directorio de entrada (DEGs)
INPUT_DIR <- file.path(BASE_DIR, "1DEGs_analysis/FC3")

# Directorio de salida
OUTPUT_DIR <- file.path(BASE_DIR, "2resultsWGCNA/07_03/FC3")

# Crear directorio de salida si no existe
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
  cat(sprintf("✓ Directorio de salida creado: %s\n", OUTPUT_DIR))
}

################################################################################
# PARÁMETROS DE ANÁLISIS
################################################################################

CONFIG <- list(
  # Archivos específicos (ya están correctos)
  control_file = "normal_DEGs_FC3.csv",
  tumor_file = "tumor_DEGs_FC3.csv",
  
  # NO necesitas degs_file porque tus archivos YA están filtrados
  degs_file = NULL,
  
  # IMPORTANTE: Aumenta este límite para tus datos
  max_expected_genes = 3000,  # Cambia de 2000 a 3000
  min_expected_genes = 100,
  
  # Umbrales de correlación
  cor_threshold_control = 0.7,   # Mínima correlación en control
  cor_threshold_tumor = 0.3,     # Máxima correlación en tumor
  min_delta_correlation = 0.3,   # Mínimo cambio en correlación
  
  # Método de correlación
  cor_method = "pearson",  # "pearson" o "spearman"
  
  # Archivos de salida
  output_interactions = "lost_interactions.csv",
  output_summary = "analysis_summary.txt",
  output_cor_control = "correlation_matrix_control.rds",
  output_cor_tumor = "correlation_matrix_tumor.rds",
  output_genes_used = "genes_used.txt"  # NUEVO: lista de genes incluidos
)

################################################################################
# FUNCIONES
################################################################################

# NUEVA FUNCIÓN: Cargar lista de DEGs
load_degs_list <- function(degs_file_path, gene_column = NULL) {
  if (is.null(degs_file_path)) {
    cat("\n⚠️  ADVERTENCIA: No se especificó archivo de DEGs\n")
    cat("   Se usarán TODOS los genes del archivo de expresión\n")
    cat("   Esto puede resultar en redes muy grandes y poco específicas\n")
    cat("   Recomendación: Especifica CONFIG$degs_file\n\n")
    return(NULL)
  }
  
  if (!file.exists(degs_file_path)) {
    stop(sprintf("❌ Archivo de DEGs no encontrado: %s", degs_file_path))
  }
  
  cat(sprintf("\n[DEGs] Cargando lista de genes: %s\n", basename(degs_file_path)))
  
  degs_data <- fread(degs_file_path, header = TRUE)
  
  cat(sprintf("  Dimensiones del archivo: %d filas × %d columnas\n", 
              nrow(degs_data), ncol(degs_data)))
  
  # Autodetectar columna de genes si no se especificó
  if (is.null(gene_column)) {
    # Nombres comunes de columnas de genes
    possible_names <- c("gene", "gene_id", "gene_name", "symbol", "Gene", 
                        "GENE", "gene_symbol", "GeneSymbol", "ensembl_gene_id")
    
    found_col <- NULL
    for (col_name in possible_names) {
      if (col_name %in% colnames(degs_data)) {
        found_col <- col_name
        break
      }
    }
    
    if (is.null(found_col)) {
      # Usar primera columna como fallback
      found_col <- colnames(degs_data)[1]
      cat(sprintf("  ⚠️  No se detectó columna estándar de genes\n"))
      cat(sprintf("  Usando primera columna: '%s'\n", found_col))
    } else {
      cat(sprintf("  ✓ Columna de genes detectada: '%s'\n", found_col))
    }
    
    gene_column <- found_col
  }
  
  if (!gene_column %in% colnames(degs_data)) {
    stop(sprintf("❌ La columna '%s' no existe en el archivo de DEGs", gene_column))
  }
  
  degs_list <- as.character(degs_data[[gene_column]])
  degs_list <- unique(degs_list[!is.na(degs_list) & degs_list != ""])
  
  cat(sprintf("  ✓ DEGs únicos cargados: %d\n", length(degs_list)))
  
  if (length(degs_list) < CONFIG$min_expected_genes) {
    cat(sprintf("  ⚠️  ADVERTENCIA: Número de DEGs muy bajo (%d)\n", length(degs_list)))
  } else if (length(degs_list) > CONFIG$max_expected_genes) {
    cat(sprintf("  ⚠️  ADVERTENCIA: Número de DEGs muy alto (%d)\n", length(degs_list)))
    cat("     ¿Estás seguro de que estos son solo DEGs significativos?\n")
  }
  
  return(degs_list)
}

# Buscar archivos en el directorio de entrada
find_input_files <- function() {
  cat("\n[BÚSQUEDA] Buscando archivos CSV en el directorio de entrada...\n")
  cat(sprintf("  Directorio: %s\n", INPUT_DIR))
  
  all_files <- list.files(INPUT_DIR, pattern = "\\.csv$", full.names = TRUE)
  
  if (length(all_files) == 0) {
    stop("❌ No se encontraron archivos CSV en el directorio de entrada")
  }
  
  cat(sprintf("\n  Archivos CSV encontrados: %d\n", length(all_files)))
  for (f in all_files) {
    cat(sprintf("    - %s\n", basename(f)))
  }
  
  # Si no se especificaron archivos, intentar encontrarlos por patrón
  control_file <- CONFIG$control_file
  tumor_file <- CONFIG$tumor_file
  degs_file <- CONFIG$degs_file
  
  if (is.null(control_file)) {
    # Buscar archivo de control
    control_candidates <- all_files[grepl("control|normal", basename(all_files), ignore.case = TRUE)]
    
    if (length(control_candidates) == 1) {
      control_file <- control_candidates[1]
      cat(sprintf("\n  ✓ Archivo CONTROL detectado: %s\n", basename(control_file)))
    } else if (length(control_candidates) > 1) {
      cat("\n  ⚠️  Múltiples candidatos para CONTROL encontrados:\n")
      for (i in seq_along(control_candidates)) {
        cat(sprintf("    %d. %s\n", i, basename(control_candidates[i])))
      }
      stop("Por favor, especifica CONFIG$control_file en el script")
    } else {
      cat("\n  📋 Archivos disponibles:\n")
      for (i in seq_along(all_files)) {
        cat(sprintf("    %d. %s\n", i, basename(all_files[i])))
      }
      stop("No se pudo detectar automáticamente el archivo CONTROL. Por favor, especifica CONFIG$control_file")
    }
  } else {
    control_file <- file.path(INPUT_DIR, control_file)
  }
  
  if (is.null(tumor_file)) {
    # Buscar archivo de tumor
    tumor_candidates <- all_files[grepl("tumor|tumour", basename(all_files), ignore.case = TRUE)]
    
    if (length(tumor_candidates) == 1) {
      tumor_file <- tumor_candidates[1]
      cat(sprintf("  ✓ Archivo TUMOR detectado: %s\n", basename(tumor_file)))
    } else if (length(tumor_candidates) > 1) {
      cat("\n  ⚠️  Múltiples candidatos para TUMOR encontrados:\n")
      for (i in seq_along(tumor_candidates)) {
        cat(sprintf("    %d. %s\n", i, basename(tumor_candidates[i])))
      }
      stop("Por favor, especifica CONFIG$tumor_file en el script")
    } else {
      stop("No se pudo detectar automáticamente el archivo TUMOR. Por favor, especifica CONFIG$tumor_file")
    }
  } else {
    tumor_file <- file.path(INPUT_DIR, tumor_file)
  }
  
  # Buscar archivo de DEGs (nuevo)
  if (is.null(degs_file)) {
    # Intentar encontrar archivo de DEGs
    degs_candidates <- all_files[grepl("deg|significant|filtered", basename(all_files), ignore.case = TRUE)]
    degs_candidates <- degs_candidates[!grepl("control|tumor|normal|tumour", 
                                              basename(degs_candidates), ignore.case = TRUE)]
    
    if (length(degs_candidates) == 1) {
      degs_file <- degs_candidates[1]
      cat(sprintf("  ✓ Archivo DEGs detectado: %s\n", basename(degs_file)))
    } else if (length(degs_candidates) > 1) {
      cat("\n  ℹ️  Múltiples candidatos para archivo DEGs encontrados:\n")
      for (i in seq_along(degs_candidates)) {
        cat(sprintf("    %d. %s\n", i, basename(degs_candidates[i])))
      }
      cat("  Especifica CONFIG$degs_file si quieres usar uno específico\n")
      degs_file <- NULL
    }
  } else {
    degs_file <- file.path(INPUT_DIR, degs_file)
  }
  
  # Verificar que existen
  if (!file.exists(control_file)) {
    stop(sprintf("❌ Archivo no encontrado: %s", control_file))
  }
  if (!file.exists(tumor_file)) {
    stop(sprintf("❌ Archivo no encontrado: %s", tumor_file))
  }
  
  return(list(control = control_file, tumor = tumor_file, degs = degs_file))
}

# MODIFICADA: Cargar datos de expresión con filtrado opcional por DEGs
load_expression_data <- function(file_path, degs_list = NULL) {
  cat(sprintf("\n[CARGA] Leyendo: %s\n", basename(file_path)))
  
  data <- fread(file_path, header = TRUE)
  
  # Detectar formato
  cat(sprintf("  Dimensiones originales: %d filas × %d columnas\n", 
              nrow(data), ncol(data)))
  
  # Asumir que primera columna son genes
  gene_col <- 1
  gene_names <- data[[gene_col]]
  
  # Resto son muestras
  expr_matrix <- as.matrix(data[, -gene_col, with = FALSE])
  rownames(expr_matrix) <- gene_names
  
  # NUEVO: Filtrar por DEGs si se proporcionó la lista
  if (!is.null(degs_list)) {
    genes_before <- nrow(expr_matrix)
    genes_to_keep <- intersect(gene_names, degs_list)
    
    if (length(genes_to_keep) == 0) {
      stop("❌ CRÍTICO: Ningún gen de la lista de DEGs se encontró en los datos de expresión")
    }
    
    expr_matrix <- expr_matrix[genes_to_keep, , drop = FALSE]
    
    cat(sprintf("  🔍 Filtrado por DEGs:\n"))
    cat(sprintf("     - Genes antes: %d\n", genes_before))
    cat(sprintf("     - Genes después: %d\n", nrow(expr_matrix)))
    cat(sprintf("     - Genes removidos: %d\n", genes_before - nrow(expr_matrix)))
    
    missing_degs <- length(degs_list) - length(genes_to_keep)
    if (missing_degs > 0) {
      cat(sprintf("     - DEGs no encontrados en datos: %d\n", missing_degs))
    }
  } else {
    cat("  ⚠️  Sin filtrado por DEGs - usando todos los genes del archivo\n")
  }
  
  # Transponer: muestras en filas, genes en columnas (formato WGCNA)
  expr_matrix <- t(expr_matrix)
  
  # Validar número de genes
  n_genes <- ncol(expr_matrix)
  if (n_genes > CONFIG$max_expected_genes) {
    cat(sprintf("\n  ⚠️  ADVERTENCIA: Número de genes muy alto (%d)\n", n_genes))
    cat("     Esto puede indicar que no se están usando solo DEGs\n")
    cat("     Considera especificar CONFIG$degs_file para filtrar\n\n")
  }
  
  # Filtrar genes con varianza cero o NA
  gene_vars <- apply(expr_matrix, 2, var, na.rm = TRUE)
  bad_genes <- which(gene_vars == 0 | is.na(gene_vars))
  
  if (length(bad_genes) > 0) {
    cat(sprintf("  Removiendo %d genes con varianza cero/NA\n", length(bad_genes)))
    expr_matrix <- expr_matrix[, -bad_genes, drop = FALSE]
  }
  
  # Verificar valores NA
  na_count <- sum(is.na(expr_matrix))
  if (na_count > 0) {
    cat(sprintf("  Advertencia: %d valores NA detectados (%.2f%%)\n", 
                na_count, 100 * na_count / length(expr_matrix)))
  }
  
  cat(sprintf("  ✓ Datos cargados: %d muestras × %d genes\n", 
              nrow(expr_matrix), ncol(expr_matrix)))
  
  return(expr_matrix)
}

# Calcular correlaciones
calculate_correlations <- function(expr_data, condition_name) {
  cat(sprintf("\n[CORRELACIÓN] Calculando para %s...\n", condition_name))
  
  n_genes <- ncol(expr_data)
  n_pairs <- choose(n_genes, 2)
  
  cat(sprintf("  Genes: %d\n", n_genes))
  cat(sprintf("  Pares posibles: %s\n", format(n_pairs, big.mark = ",")))
  
  # Advertencia si hay demasiados pares
  if (n_pairs > 50000000) {  # 50 millones
    cat(sprintf("  ⚠️  ADVERTENCIA: Número muy alto de pares (%s)\n", 
                format(n_pairs, big.mark = ",")))
    cat("     El cálculo puede tomar mucho tiempo\n")
  }
  
  cor_matrix <- cor(expr_data, 
                    method = CONFIG$cor_method, 
                    use = "pairwise.complete.obs")
  
  cat("  ✓ Correlación completada\n")
  
  return(cor_matrix)
}

# Extraer interacciones perdidas
extract_lost_interactions <- function(cor_control, cor_tumor) {
  cat("\n[EXTRACCIÓN] Identificando interacciones perdidas...\n")
  cat(sprintf("  Umbral control: >= %.2f\n", CONFIG$cor_threshold_control))
  cat(sprintf("  Umbral tumor: <= %.2f\n", CONFIG$cor_threshold_tumor))
  cat(sprintf("  Delta mínimo: >= %.2f\n", CONFIG$min_delta_correlation))
  
  # Genes comunes
  common_genes <- intersect(colnames(cor_control), colnames(cor_tumor))
  n_genes <- length(common_genes)
  
  cat(sprintf("  Genes comunes: %d\n", n_genes))
  
  if (n_genes == 0) {
    stop("❌ No hay genes comunes entre control y tumor")
  }
  
  # Subset matrices
  cor_control <- cor_control[common_genes, common_genes]
  cor_tumor <- cor_tumor[common_genes, common_genes]
  
  # Pre-asignar vectores
  max_interactions <- choose(n_genes, 2)
  gene1_vec <- character(max_interactions)
  gene2_vec <- character(max_interactions)
  cor_ctrl_vec <- numeric(max_interactions)
  cor_tmr_vec <- numeric(max_interactions)
  
  counter <- 0
  
  cat("  Procesando pares...\n")
  pb <- txtProgressBar(min = 0, max = n_genes - 1, style = 3, width = 50)
  
  for (i in 1:(n_genes - 1)) {
    j_indices <- (i + 1):n_genes
    
    cor_ctrl_vals <- cor_control[i, j_indices]
    cor_tmr_vals <- cor_tumor[i, j_indices]
    
    # Aplicar filtros
    valid_pairs <- which(
      !is.na(cor_ctrl_vals) & 
        !is.na(cor_tmr_vals) &
        cor_ctrl_vals >= CONFIG$cor_threshold_control &
        cor_tmr_vals <= CONFIG$cor_threshold_tumor &
        (cor_ctrl_vals - cor_tmr_vals) >= CONFIG$min_delta_correlation
    )
    
    if (length(valid_pairs) > 0) {
      actual_j <- j_indices[valid_pairs]
      n_valid <- length(valid_pairs)
      
      idx_start <- counter + 1
      idx_end <- counter + n_valid
      
      gene1_vec[idx_start:idx_end] <- common_genes[i]
      gene2_vec[idx_start:idx_end] <- common_genes[actual_j]
      cor_ctrl_vec[idx_start:idx_end] <- cor_ctrl_vals[valid_pairs]
      cor_tmr_vec[idx_start:idx_end] <- cor_tmr_vals[valid_pairs]
      
      counter <- counter + n_valid
    }
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  cat(sprintf("\n\n  ✓ Interacciones encontradas: %s\n", format(counter, big.mark = ",")))
  
  if (counter == 0) {
    cat("\n  ⚠️  ADVERTENCIA: No se encontraron interacciones que cumplan los criterios\n")
    cat("     Considera:\n")
    cat(sprintf("       - Reducir cor_threshold_control (actualmente %.2f)\n", CONFIG$cor_threshold_control))
    cat(sprintf("       - Aumentar cor_threshold_tumor (actualmente %.2f)\n", CONFIG$cor_threshold_tumor))
    cat(sprintf("       - Reducir min_delta_correlation (actualmente %.2f)\n", CONFIG$min_delta_correlation))
    return(NULL)
  }
  
  # Crear data frame
  interactions_df <- data.frame(
    Gene1 = gene1_vec[1:counter],
    Gene2 = gene2_vec[1:counter],
    Correlation_Control = cor_ctrl_vec[1:counter],
    Correlation_Tumor = cor_tmr_vec[1:counter],
    Delta_Correlation = cor_ctrl_vec[1:counter] - cor_tmr_vec[1:counter],
    stringsAsFactors = FALSE
  )
  
  # Ordenar por delta
  interactions_df <- interactions_df[order(-interactions_df$Delta_Correlation), ]
  
  return(interactions_df)
}

# MODIFICADA: Generar resumen (sin hardcoded dataset name)
generate_summary <- function(lost_interactions, input_files, genes_used) {
  summary_file <- file.path(OUTPUT_DIR, CONFIG$output_summary)
  
  sink(summary_file)
  
  cat(rep("=", 80), "\n", sep = "")
  cat("ANÁLISIS WGCNA - INTERACCIONES PERDIDAS\n")
  cat(rep("=", 80), "\n", sep = "")
  cat(sprintf("\nFecha: %s\n", Sys.time()))
  cat(sprintf("Directorio de análisis: %s\n", BASE_DIR))
  
  cat("\nARCHIVOS DE ENTRADA:\n")
  cat(sprintf("  Control: %s\n", basename(input_files$control)))
  cat(sprintf("  Tumor:   %s\n", basename(input_files$tumor)))
  if (!is.null(input_files$degs)) {
    cat(sprintf("  DEGs:    %s\n", basename(input_files$degs)))
  } else {
    cat("  DEGs:    [No especificado - se usaron todos los genes]\n")
  }
  
  cat("\nGENES UTILIZADOS:\n")
  cat(sprintf("  Total de genes en el análisis: %d\n", length(genes_used)))
  if (!is.null(input_files$degs)) {
    cat("  Filtrado: Sí (solo DEGs)\n")
  } else {
    cat("  Filtrado: No (todos los genes del archivo)\n")
  }
  
  cat("\nPARÁMETROS:\n")
  cat(sprintf("  Umbral correlación control: >= %.2f\n", CONFIG$cor_threshold_control))
  cat(sprintf("  Umbral correlación tumor:   <= %.2f\n", CONFIG$cor_threshold_tumor))
  cat(sprintf("  Delta mínimo requerido:     >= %.2f\n", CONFIG$min_delta_correlation))
  cat(sprintf("  Método de correlación:      %s\n", CONFIG$cor_method))
  
  if (!is.null(lost_interactions)) {
    cat("\nRESULTADOS:\n")
    cat(sprintf("  Total interacciones perdidas: %s\n", 
                format(nrow(lost_interactions), big.mark = ",")))
    
    cat("\nESTADÍSTICAS DE CORRELACIÓN:\n")
    cat(sprintf("  Control - Media: %.3f | Rango: [%.3f, %.3f]\n",
                mean(lost_interactions$Correlation_Control),
                min(lost_interactions$Correlation_Control),
                max(lost_interactions$Correlation_Control)))
    cat(sprintf("  Tumor   - Media: %.3f | Rango: [%.3f, %.3f]\n",
                mean(lost_interactions$Correlation_Tumor),
                min(lost_interactions$Correlation_Tumor),
                max(lost_interactions$Correlation_Tumor)))
    cat(sprintf("  Delta   - Media: %.3f | Rango: [%.3f, %.3f]\n",
                mean(lost_interactions$Delta_Correlation),
                min(lost_interactions$Delta_Correlation),
                max(lost_interactions$Delta_Correlation)))
    
    cat("\nTOP 20 INTERACCIONES:\n")
    cat(rep("-", 80), "\n", sep = "")
    top20 <- head(lost_interactions, 20)
    for (i in 1:nrow(top20)) {
      cat(sprintf("%2d. %-15s -- %-15s | Control: %.3f | Tumor: %.3f | Δ: %.3f\n",
                  i, top20$Gene1[i], top20$Gene2[i],
                  top20$Correlation_Control[i],
                  top20$Correlation_Tumor[i],
                  top20$Delta_Correlation[i]))
    }
  } else {
    cat("\nRESULTADOS:\n")
    cat("  No se encontraron interacciones que cumplan los criterios\n")
  }
  
  cat("\nARCHIVOS DE SALIDA:\n")
  cat(sprintf("  %s\n", CONFIG$output_interactions))
  cat(sprintf("  %s\n", CONFIG$output_summary))
  cat(sprintf("  %s\n", CONFIG$output_cor_control))
  cat(sprintf("  %s\n", CONFIG$output_cor_tumor))
  cat(sprintf("  %s\n", CONFIG$output_genes_used))
  
  cat("\n", rep("=", 80), "\n", sep = "")
  
  sink()
  
  cat(sprintf("\n✓ Resumen guardado en: %s\n", summary_file))
}

################################################################################
# MAIN
################################################################################

main <- function() {
  cat("\n")
  cat(rep("=", 80), "\n", sep = "")
  cat("ANÁLISIS WGCNA - INTERACCIONES PERDIDAS\n")
  cat("Detección de correlaciones perdidas en tumor vs control\n")
  cat(rep("=", 80), "\n", sep = "")
  
  # Verificar directorios
  if (!dir.exists(INPUT_DIR)) {
    stop(sprintf("❌ Directorio de entrada no existe: %s", INPUT_DIR))
  }
  
  cat(sprintf("\n📁 Directorio de entrada: %s\n", INPUT_DIR))
  cat(sprintf("📁 Directorio de salida:  %s\n", OUTPUT_DIR))
  
  # 1. Buscar archivos
  cat("\n[1/6] BÚSQUEDA DE ARCHIVOS\n")
  cat(rep("-", 80), "\n", sep = "")
  input_files <- find_input_files()
  
  # 1.5. Cargar lista de DEGs (NUEVO)
  cat("\n[2/6] CARGA DE LISTA DE DEGs\n")
  cat(rep("-", 80), "\n", sep = "")
  degs_list <- load_degs_list(input_files$degs, CONFIG$degs_gene_column)
  
  # 2. Cargar datos
  cat("\n[3/6] CARGA DE DATOS DE EXPRESIÓN\n")
  cat(rep("-", 80), "\n", sep = "")
  expr_control <- load_expression_data(input_files$control, degs_list)
  expr_tumor <- load_expression_data(input_files$tumor, degs_list)
  
  # 3. Control de calidad
  cat("\n[4/6] CONTROL DE CALIDAD\n")
  cat(rep("-", 80), "\n", sep = "")
  
  gsg_control <- goodSamplesGenes(expr_control, verbose = 0)
  gsg_tumor <- goodSamplesGenes(expr_tumor, verbose = 0)
  
  if (!gsg_control$allOK) {
    cat("  Limpiando datos control...\n")
    expr_control <- expr_control[gsg_control$goodSamples, gsg_control$goodGenes]
  } else {
    cat("  ✓ Control OK\n")
  }
  
  if (!gsg_tumor$allOK) {
    cat("  Limpiando datos tumor...\n")
    expr_tumor <- expr_tumor[gsg_tumor$goodSamples, gsg_tumor$goodGenes]
  } else {
    cat("  ✓ Tumor OK\n")
  }
  
  # Obtener lista final de genes usados
  common_genes <- intersect(colnames(expr_control), colnames(expr_tumor))
  cat(sprintf("\n  Genes comunes finales: %d\n", length(common_genes)))
  
  # Guardar lista de genes utilizados
  genes_used_file <- file.path(OUTPUT_DIR, CONFIG$output_genes_used)
  writeLines(sort(common_genes), genes_used_file)
  cat(sprintf("  ✓ Lista de genes guardada en: %s\n", basename(genes_used_file)))
  
  # 4. Calcular correlaciones
  cat("\n[5/6] CÁLCULO DE CORRELACIONES\n")
  cat(rep("-", 80), "\n", sep = "")
  
  cor_control <- calculate_correlations(expr_control, "CONTROL")
  cor_tumor <- calculate_correlations(expr_tumor, "TUMOR")
  
  # Guardar matrices de correlación
  saveRDS(cor_control, file.path(OUTPUT_DIR, CONFIG$output_cor_control))
  saveRDS(cor_tumor, file.path(OUTPUT_DIR, CONFIG$output_cor_tumor))
  cat(sprintf("\n  ✓ Matrices guardadas en: %s\n", OUTPUT_DIR))
  
  # 5. Extraer interacciones
  cat("\n[6/6] EXTRACCIÓN DE INTERACCIONES PERDIDAS\n")
  cat(rep("-", 80), "\n", sep = "")
  
  lost_interactions <- extract_lost_interactions(cor_control, cor_tumor)
  
  # Guardar resultados
  if (!is.null(lost_interactions)) {
    output_file <- file.path(OUTPUT_DIR, CONFIG$output_interactions)
    write.csv(lost_interactions, output_file, row.names = FALSE, quote = FALSE)
    cat(sprintf("\n✓ Interacciones guardadas en: %s\n", output_file))
  }
  
  # Generar resumen
  generate_summary(lost_interactions, input_files, common_genes)
  
  cat("\n")
  cat(rep("=", 80), "\n", sep = "")
  cat("✨ ANÁLISIS COMPLETADO\n")
  cat(rep("=", 80), "\n", sep = "")
  
  # Resumen final en consola
  cat("\n📊 RESUMEN FINAL:\n")
  cat(sprintf("   - Genes analizados: %d\n", length(common_genes)))
  if (!is.null(lost_interactions)) {
    cat(sprintf("   - Interacciones perdidas: %s\n", 
                format(nrow(lost_interactions), big.mark = ",")))
  } else {
    cat("   - Interacciones perdidas: 0\n")
  }
  cat(sprintf("   - Resultados en: %s\n\n", OUTPUT_DIR))
}

# Ejecutar
tryCatch({
  main()
}, error = function(e) {
  cat("\n❌ ERROR:\n")
  cat(conditionMessage(e), "\n\n")
  traceback()
  quit(status = 1)
})

