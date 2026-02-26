################################################################################
# EVALUACIÓN DE RED DE INTERACCIONES PERDIDAS vs GENEMANIA
################################################################################
# Este script compara las interacciones perdidas detectadas por WGCNA
# con una red GoldStandard de GeneMania y calcula métricas de evaluación
################################################################################

library(data.table)

################################################################################
# CONFIGURACIÓN GLOBAL (no cambia entre ejecuciones)
################################################################################

CONFIG <- list(
  # Directorio base
  base_dir = "/Users/mriosc/Documents/NETWORKS/TCGA-STAD/by_Molecular_Subtype/MSI",
  
  # ► AQUÍ defines todos los FC que quieres evaluar
  fc_list = c("FC2", "FC3"),
  
  # Carpeta base de resultados WGCNA (el FC se añade automáticamente)
  wgcna_results_base = "2resultsWGCNA/06_03",
  
  # Archivos de entrada
  lost_interactions_file = "lost_interactions.csv",
  genemania_file = "/Users/mriosc/Documents/NETWORKS/red_coexpresion_final.txt",
  
  # Columnas en el archivo de GeneMania
  genemania_gene1_col = "Gene1",
  genemania_gene2_col = "Gene2",
  genemania_weight_col = "Weight",
  
  # Filtro de peso (NULL = sin filtro)
  genemania_min_weight = NULL,
  
  # Nombres de archivos de salida (se añadirá el sufijo _FC2, _FC3, etc.)
  output_report  = "network_evaluation_report.txt",
  output_metrics = "evaluation_metrics.csv",
  output_confusion = "confusion_matrix.csv"
)

################################################################################
# FUNCIONES (sin cambios respecto al script original)
################################################################################

create_gene_pair <- function(gene1, gene2) {
  pair <- ifelse(gene1 < gene2,
                 paste(gene1, gene2, sep = "~"),
                 paste(gene2, gene1, sep = "~"))
  return(pair)
}

load_lost_interactions <- function(filepath) {
  cat("\n[1/4] CARGANDO INTERACCIONES PERDIDAS\n")
  cat(rep("-", 80), "\n", sep = "")
  
  if (!file.exists(filepath)) stop(sprintf("❌ Archivo no encontrado: %s", filepath))
  
  cat(sprintf("  Archivo: %s\n", basename(filepath)))
  data <- fread(filepath, header = TRUE)
  cat(sprintf("  ✓ Interacciones cargadas: %s\n", format(nrow(data), big.mark = ",")))
  
  data$pair_id <- create_gene_pair(data$Gene1, data$Gene2)
  genes <- unique(c(data$Gene1, data$Gene2))
  cat(sprintf("  ✓ Genes únicos: %d\n", length(genes)))
  
  return(list(interactions = data, genes = genes, pairs = unique(data$pair_id)))
}

load_genemania_network <- function(filepath, gene1_col, gene2_col) {
  cat("\n[2/4] CARGANDO RED GOLDSTANDARD (GENEMANIA)\n")
  cat(rep("-", 80), "\n", sep = "")
  
  if (!file.exists(filepath)) stop(sprintf("❌ Archivo no encontrado: %s", filepath))
  
  cat(sprintf("  Archivo: %s\n", basename(filepath)))
  
  first_line <- readLines(filepath, n = 1)
  if (grepl("\t", first_line)) {
    cat("  Separador detectado: TAB\n")
    data <- fread(filepath, header = TRUE, sep = "\t")
  } else if (grepl(",", first_line)) {
    cat("  Separador detectado: COMA\n")
    data <- fread(filepath, header = TRUE, sep = ",")
  } else {
    cat("  Usando separador por defecto\n")
    data <- fread(filepath, header = TRUE)
  }
  
  cat(sprintf("  Dimensiones originales: %d filas × %d columnas\n", nrow(data), ncol(data)))
  cat(sprintf("  Columnas: %s\n", paste(colnames(data), collapse = ", ")))
  
  if (!gene1_col %in% colnames(data))
    stop(sprintf("❌ Columna '%s' no encontrada. Disponibles: %s", gene1_col, paste(colnames(data), collapse = ", ")))
  if (!gene2_col %in% colnames(data))
    stop(sprintf("❌ Columna '%s' no encontrada. Disponibles: %s", gene2_col, paste(colnames(data), collapse = ", ")))
  
  genemania_pairs <- data.frame(
    Gene1 = as.character(data[[gene1_col]]),
    Gene2 = as.character(data[[gene2_col]]),
    stringsAsFactors = FALSE
  )
  
  if (!is.null(CONFIG$genemania_weight_col) && CONFIG$genemania_weight_col %in% colnames(data)) {
    genemania_pairs$Weight <- as.numeric(data[[CONFIG$genemania_weight_col]])
    cat(sprintf("  ✓ Columna de pesos: %s | Rango: [%.4f, %.4f]\n",
                CONFIG$genemania_weight_col,
                min(genemania_pairs$Weight, na.rm = TRUE),
                max(genemania_pairs$Weight, na.rm = TRUE)))
    
    if (!is.null(CONFIG$genemania_min_weight)) {
      n_before <- nrow(genemania_pairs)
      genemania_pairs <- genemania_pairs[genemania_pairs$Weight >= CONFIG$genemania_min_weight, ]
      cat(sprintf("  ✓ Filtrado por peso >= %.4f: %s → %s interacciones\n",
                  CONFIG$genemania_min_weight,
                  format(n_before, big.mark = ","),
                  format(nrow(genemania_pairs), big.mark = ",")))
    }
  }
  
  genemania_pairs$pair_id <- create_gene_pair(genemania_pairs$Gene1, genemania_pairs$Gene2)
  
  cat(sprintf("  ✓ Interacciones totales: %s\n", format(length(unique(genemania_pairs$pair_id)), big.mark = ",")))
  
  genes <- unique(c(genemania_pairs$Gene1, genemania_pairs$Gene2))
  cat(sprintf("  ✓ Genes únicos: %s\n", format(length(genes), big.mark = ",")))
  
  return(list(interactions = genemania_pairs, genes = genes, pairs = unique(genemania_pairs$pair_id)))
}

filter_genemania_by_genes <- function(genemania_data, target_genes) {
  cat("\n[3/4] FILTRANDO GENEMANIA POR GENES DE ANÁLISIS\n")
  cat(rep("-", 80), "\n", sep = "")
  
  cat(sprintf("  Genes objetivo: %d\n", length(target_genes)))
  
  filtered <- genemania_data$interactions[
    genemania_data$interactions$Gene1 %in% target_genes &
      genemania_data$interactions$Gene2 %in% target_genes, ]
  
  cat(sprintf("  ✓ Interacciones filtradas: %s\n", format(length(unique(filtered$pair_id)), big.mark = ",")))
  
  genes_in_both <- intersect(genemania_data$genes, target_genes)
  cat(sprintf("  ✓ Genes comunes: %d\n", length(genes_in_both)))
  
  n_genes <- length(genes_in_both)
  max_possible_interactions <- (n_genes * (n_genes - 1)) / 2
  cat(sprintf("  ℹ️  Máximo interacciones posibles: %s\n", format(max_possible_interactions, big.mark = ",")))
  
  return(list(interactions = filtered, pairs = unique(filtered$pair_id),
              genes = genes_in_both, max_possible = max_possible_interactions))
}

calculate_metrics <- function(predicted_pairs, goldstandard_pairs, all_genes, verbose = TRUE) {
  cat("\n[4/4] CALCULANDO MÉTRICAS DE EVALUACIÓN\n")
  cat(rep("-", 80), "\n", sep = "")
  
  n_genes <- length(all_genes)
  total_possible <- (n_genes * (n_genes - 1)) / 2
  
  TP_pairs <- intersect(predicted_pairs, goldstandard_pairs); TP <- length(TP_pairs)
  FP_pairs <- setdiff(predicted_pairs, goldstandard_pairs);   FP <- length(FP_pairs)
  FN_pairs <- setdiff(goldstandard_pairs, predicted_pairs);   FN <- length(FN_pairs)
  TN <- total_possible - TP - FP - FN
  
  if (verbose) {
    cat("\n  MATRIZ DE CONFUSIÓN:\n")
    cat(sprintf("    True Positives (TP):   %s\n", format(TP, big.mark = ",")))
    cat(sprintf("    False Positives (FP):  %s\n", format(FP, big.mark = ",")))
    cat(sprintf("    False Negatives (FN):  %s\n", format(FN, big.mark = ",")))
    cat(sprintf("    True Negatives (TN):   %s\n", format(TN, big.mark = ",")))
    cat(sprintf("    Total posible:         %s\n", format(total_possible, big.mark = ",")))
  }
  
  precision    <- ifelse(TP + FP > 0, TP / (TP + FP), 0)
  recall       <- ifelse(TP + FN > 0, TP / (TP + FN), 0)
  f_score      <- ifelse(precision + recall > 0, 2 * (precision * recall) / (precision + recall), 0)
  accuracy     <- (TP + TN) / total_possible
  specificity  <- ifelse(TN + FP > 0, TN / (TN + FP), 0)
  
  if (verbose) {
    cat("\n  MÉTRICAS DE EVALUACIÓN:\n")
    cat(sprintf("    Precision:    %.4f (%.2f%%)\n", precision, precision * 100))
    cat(sprintf("    Recall:       %.4f (%.2f%%)\n", recall, recall * 100))
    cat(sprintf("    F-score:      %.4f\n", f_score))
    cat(sprintf("    Accuracy:     %.4f (%.2f%%)\n", accuracy, accuracy * 100))
    cat(sprintf("    Specificity:  %.4f (%.2f%%)\n", specificity, specificity * 100))
  }
  
  return(list(
    confusion_matrix = data.frame(
      Metric = c("TP", "FP", "FN", "TN", "Total_Possible"),
      Value  = c(TP, FP, FN, TN, total_possible),
      stringsAsFactors = FALSE
    ),
    metrics = data.frame(
      Metric     = c("Precision", "Recall", "F_score", "Accuracy", "Specificity"),
      Value      = c(precision, recall, f_score, accuracy, specificity),
      Percentage = c(precision * 100, recall * 100, f_score * 100, accuracy * 100, specificity * 100),
      stringsAsFactors = FALSE
    ),
    TP_pairs = TP_pairs,
    FP_pairs = FP_pairs,
    FN_pairs = FN_pairs
  ))
}

generate_report <- function(lost_data, genemania_data, filtered_gm, metrics, output_dir, fc_label) {
  # Añadir sufijo del FC al nombre del archivo para no sobreescribir
  report_filename <- sub("\\.txt$", sprintf("_%s.txt", fc_label), CONFIG$output_report)
  report_file <- file.path(output_dir, report_filename)
  
  sink(report_file)
  
  cat(rep("=", 80), "\n", sep = "")
  cat(sprintf("EVALUACIÓN DE RED DE INTERACCIONES PERDIDAS vs GENEMANIA — %s\n", fc_label))
  cat(rep("=", 80), "\n", sep = "")
  cat(sprintf("\nFecha: %s\n", Sys.time()))
  
  cat("\n", rep("-", 80), "\n", sep = "")
  cat("DATOS DE ENTRADA\n")
  cat(rep("-", 80), "\n", sep = "")
  
  cat("\nRED PREDICHA (Interacciones Perdidas):\n")
  cat(sprintf("  Genes únicos:           %d\n", length(lost_data$genes)))
  cat(sprintf("  Interacciones totales:  %s\n", format(length(lost_data$pairs), big.mark = ",")))
  
  cat("\nRED GOLDSTANDARD (GeneMania - Original):\n")
  cat(sprintf("  Genes únicos:           %s\n", format(length(genemania_data$genes), big.mark = ",")))
  cat(sprintf("  Interacciones totales:  %s\n", format(length(genemania_data$pairs), big.mark = ",")))
  
  cat("\nRED GOLDSTANDARD (GeneMania - Filtrada):\n")
  cat(sprintf("  Genes comunes:          %d\n", length(filtered_gm$genes)))
  cat(sprintf("  Interacciones filtradas: %s\n", format(length(filtered_gm$pairs), big.mark = ",")))
  cat(sprintf("  Máximo posible:         %s\n", format(filtered_gm$max_possible, big.mark = ",")))
  
  cat("\n", rep("-", 80), "\n", sep = "")
  cat("MATRIZ DE CONFUSIÓN\n")
  cat(rep("-", 80), "\n", sep = "")
  cm <- metrics$confusion_matrix
  for (i in 1:nrow(cm)) cat(sprintf("  %-20s %s\n", paste0(cm$Metric[i], ":"), format(cm$Value[i], big.mark = ",")))
  
  cat("\n", rep("-", 80), "\n", sep = "")
  cat("MÉTRICAS DE EVALUACIÓN\n")
  cat(rep("-", 80), "\n", sep = "")
  m <- metrics$metrics
  for (i in 1:nrow(m)) cat(sprintf("  %-20s %.6f  (%.2f%%)\n", paste0(m$Metric[i], ":"), m$Value[i], m$Percentage[i]))
  
  cat("\n", rep("-", 80), "\n", sep = "")
  cat("INTERPRETACIÓN\n")
  cat(rep("-", 80), "\n", sep = "")
  precision <- m$Value[m$Metric == "Precision"]
  recall    <- m$Value[m$Metric == "Recall"]
  f_score   <- m$Value[m$Metric == "F_score"]
  
  cat("\nPRECISION (", sprintf("%.2f%%", precision * 100), "):\n", sep = "")
  if      (precision >= 0.7) cat("  ✓ EXCELENTE - La mayoría de las predicciones están validadas\n")
  else if (precision >= 0.5) cat("  ✓ BUENA - Más de la mitad de las predicciones son correctas\n")
  else if (precision >= 0.3) cat("  ⚠️  MODERADA - Hay bastantes falsos positivos\n")
  else                       cat("  ⚠️  BAJA - Muchas predicciones no están validadas en GeneMania\n")
  
  cat("\nRECALL (", sprintf("%.2f%%", recall * 100), "):\n", sep = "")
  if      (recall >= 0.7) cat("  ✓ EXCELENTE - Se detecta la mayoría de las interacciones conocidas\n")
  else if (recall >= 0.5) cat("  ✓ BUENA - Se detecta más de la mitad de las interacciones conocidas\n")
  else if (recall >= 0.3) cat("  ⚠️  MODERADA - Se pierden bastantes interacciones conocidas\n")
  else                    cat("  ⚠️  BAJA - Se detectan pocas de las interacciones conocidas\n")
  
  cat("\nF-SCORE (", sprintf("%.4f", f_score), "):\n", sep = "")
  if      (f_score >= 0.7) cat("  ✓ EXCELENTE balance entre precisión y recall\n")
  else if (f_score >= 0.5) cat("  ✓ BUEN balance entre precisión y recall\n")
  else if (f_score >= 0.3) cat("  ⚠️  Balance MODERADO entre precisión y recall\n")
  else                     cat("  ⚠️  Balance BAJO - considerar ajustar umbrales\n")
  
  cat("\n", rep("-", 80), "\n", sep = "")
  cat("EJEMPLOS DE INTERACCIONES\n")
  cat(rep("-", 80), "\n", sep = "")
  
  if (length(metrics$TP_pairs) > 0) {
    cat("\nTRUE POSITIVES (Primeras 10):\n")
    for (pair in head(metrics$TP_pairs, 10)) {
      g <- strsplit(pair, "~")[[1]]
      cat(sprintf("  ✓ %s -- %s\n", g[1], g[2]))
    }
  }
  if (length(metrics$FP_pairs) > 0) {
    cat("\nFALSE POSITIVES (Primeras 10):\n")
    for (pair in head(metrics$FP_pairs, 10)) {
      g <- strsplit(pair, "~")[[1]]
      cat(sprintf("  ✗ %s -- %s\n", g[1], g[2]))
    }
  }
  if (length(metrics$FN_pairs) > 0) {
    cat("\nFALSE NEGATIVES (Primeras 10):\n")
    for (pair in head(metrics$FN_pairs, 10)) {
      g <- strsplit(pair, "~")[[1]]
      cat(sprintf("  ⊗ %s -- %s\n", g[1], g[2]))
    }
  }
  
  cat("\n", rep("=", 80), "\n", sep = "")
  sink()
  
  cat(sprintf("\n  ✓ Reporte guardado: %s\n", report_filename))
}

################################################################################
# FUNCIÓN PRINCIPAL PARA UN SOLO FC
################################################################################

run_for_fc <- function(fc_label, genemania_data) {
  
  cat("\n")
  cat(rep("█", 80), "\n", sep = "")
  cat(sprintf("  PROCESANDO: %s\n", fc_label))
  cat(rep("█", 80), "\n", sep = "")
  
  # Construir ruta específica de este FC
  wgcna_dir    <- file.path(CONFIG$base_dir, CONFIG$wgcna_results_base, fc_label)
  lost_int_file <- file.path(wgcna_dir, CONFIG$lost_interactions_file)
  
  if (!file.exists(lost_int_file)) {
    cat(sprintf("⚠️  SKIPPING %s — Archivo no encontrado: %s\n", fc_label, lost_int_file))
    return(NULL)
  }
  
  # 1. Cargar interacciones perdidas de este FC
  lost_data <- load_lost_interactions(lost_int_file)
  
  # 2. GeneMania ya está cargada (se reutiliza)
  
  # 3. Filtrar GeneMania
  filtered_gm <- filter_genemania_by_genes(genemania_data, lost_data$genes)
  
  # 4. Calcular métricas
  metrics <- calculate_metrics(
    predicted_pairs    = lost_data$pairs,
    goldstandard_pairs = filtered_gm$pairs,
    all_genes          = filtered_gm$genes,
    verbose            = TRUE
  )
  
  # 5. Guardar resultados (con sufijo del FC para no sobreescribir)
  cat("\n[GUARDANDO RESULTADOS]\n")
  cat(rep("-", 80), "\n", sep = "")
  
  metrics_filename  <- sub("\\.csv$", sprintf("_%s.csv", fc_label), CONFIG$output_metrics)
  confusion_filename <- sub("\\.csv$", sprintf("_%s.csv", fc_label), CONFIG$output_confusion)
  
  write.csv(metrics$metrics,          file.path(wgcna_dir, metrics_filename),   row.names = FALSE, quote = FALSE)
  write.csv(metrics$confusion_matrix, file.path(wgcna_dir, confusion_filename), row.names = FALSE, quote = FALSE)
  
  cat(sprintf("  ✓ Métricas:          %s\n", metrics_filename))
  cat(sprintf("  ✓ Matriz confusión:  %s\n", confusion_filename))
  
  # 6. Reporte
  generate_report(lost_data, genemania_data, filtered_gm, metrics, wgcna_dir, fc_label)
  
  # Devolver resumen para la tabla comparativa final
  return(data.frame(
    FC         = fc_label,
    Precision  = metrics$metrics$Value[1],
    Recall     = metrics$metrics$Value[2],
    F_score    = metrics$metrics$Value[3],
    Accuracy   = metrics$metrics$Value[4],
    stringsAsFactors = FALSE
  ))
}

################################################################################
# MAIN — itera sobre todos los FC definidos en CONFIG$fc_list
################################################################################

main <- function() {
  
  cat("\n")
  cat(rep("=", 80), "\n", sep = "")
  cat("EVALUACIÓN AUTOMATIZADA DE REDES vs GENEMANIA\n")
  cat(sprintf("FC a procesar: %s\n", paste(CONFIG$fc_list, collapse = ", ")))
  cat(rep("=", 80), "\n", sep = "")
  
  # Verificar archivo GeneMania
  if (!file.exists(CONFIG$genemania_file)) {
    stop(sprintf("❌ Archivo GeneMania no encontrado: %s", CONFIG$genemania_file))
  }
  
  # Cargar GeneMania UNA SOLA VEZ (ahorra tiempo si son muchos FC)
  genemania_data <- load_genemania_network(
    CONFIG$genemania_file,
    CONFIG$genemania_gene1_col,
    CONFIG$genemania_gene2_col
  )
  
  # Iterar sobre cada FC
  resultados <- list()
  for (fc in CONFIG$fc_list) {
    resultados[[fc]] <- tryCatch(
      run_for_fc(fc, genemania_data),
      error = function(e) {
        cat(sprintf("\n❌ ERROR en %s: %s\n", fc, conditionMessage(e)))
        NULL
      }
    )
  }
  
  # ── TABLA COMPARATIVA FINAL ──────────────────────────────────────────────────
  cat("\n")
  cat(rep("=", 80), "\n", sep = "")
  cat("📊 RESUMEN COMPARATIVO\n")
  cat(rep("=", 80), "\n", sep = "")
  
  resumen <- do.call(rbind, Filter(Negate(is.null), resultados))
  
  if (!is.null(resumen) && nrow(resumen) > 0) {
    cat(sprintf("\n  %-6s  %-10s  %-10s  %-10s  %-10s\n",
                "FC", "Precision", "Recall", "F-score", "Accuracy"))
    cat(rep("-", 52), "\n", sep = "")
    for (i in seq_len(nrow(resumen))) {
      cat(sprintf("  %-6s  %-10s  %-10s  %-10s  %-10s\n",
                  resumen$FC[i],
                  sprintf("%.4f", resumen$Precision[i]),
                  sprintf("%.4f", resumen$Recall[i]),
                  sprintf("%.4f", resumen$F_score[i]),
                  sprintf("%.4f", resumen$Accuracy[i])))
    }
    
    # Guardar tabla comparativa en el directorio base
    summary_file <- file.path(CONFIG$base_dir,
                              CONFIG$wgcna_results_base,
                              "comparison_summary.csv")
    write.csv(resumen, summary_file, row.names = FALSE, quote = FALSE)
    cat(sprintf("\n  ✓ Tabla comparativa guardada: comparison_summary.csv\n"))
  }
  
  cat("\n")
  cat(rep("=", 80), "\n", sep = "")
  cat("✨ EVALUACIÓN COMPLETADA\n")
  cat(rep("=", 80), "\n", sep = "")
}

# Ejecutar
tryCatch({
  main()
}, error = function(e) {
  cat("\n❌ ERROR FATAL:\n")
  cat(conditionMessage(e), "\n\n")
  traceback()
  quit(status = 1)
})

