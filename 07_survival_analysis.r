################################################################################
# ANÁLISIS DE SUPERVIVENCIA - INTERACCIONES PERDIDAS COMUNES
################################################################################
# Este script analiza si la pérdida de correlación entre genes se asocia
# con supervivencia en pacientes con cáncer de TCGA
################################################################################

library(survival)
library(survminer)
library(dplyr)
library(readr)
library(ggplot2)
library(data.table)
library(gridExtra)

################################################################################
# CONFIGURACIÓN
################################################################################

CONFIG <- list(
  # Directorio base
  base_dir = "/Users/mriosc/Documents/NETWORKS",
  
  # Directorio de interacciones comunes
  common_interactions_dir = "_Common_Interactions/07_03/FC2",
  
  # Archivos de interacciones comunes (pares)
  interaction_files = c(
    "02_PAR_TCGA-LUAD_vs_TCGA-BRCA.csv",
    "02_PAR_TCGA-LUAD_vs_TCGA-STAD.csv",
    "02_PAR_TCGA-HNSC_vs_TCGA-STAD.csv"
  ),
  
  # Información de datasets por proyecto
  # Formato: Proyecto -> lista de subdatasets
  datasets = list(
    "TCGA-BRCA" = list(
      base_path = "TCGA-BRCA/by_Subtype",
      subdatasets = c("Basal", "Her2", "LumA", "LumB"),
      metadata_file = "TCGA-BRCA/TCGA-BRCA_metadata.tsv"
    ),
    "TCGA-LUAD" = list(
      base_path = "TCGA-LUAD/by_Tissue",
      subdatasets = c("Upper_lobe__lung", "Lower_lobe__lung"),
      metadata_file = "TCGA-LUAD/TCGA-LUAD_metadata.tsv"
    ),
    "TCGA-HNSC" = list(
      base_path = "TCGA-HNSC/by_Anatomic_Region",
      subdatasets = c("Larynx", "Oral_Cavity", "Oropharynx"),
      metadata_file = "TCGA-HNSC/TCGA-HNSC_metadata.tsv"
    ),
    "TCGA-STAD" = list(
      base_path = "TCGA-STAD/by_Molecular_Subtype",
      subdatasets = c("MSI", "CIN"),
      metadata_file = "TCGA-STAD/TCGA-STAD_metadata.tsv"
    )
  ),
  
  # Nombre del archivo de expresión en cada subdataset
  expression_filename = "counts_renombrado.csv",
  
  # Parámetros de análisis
  stratification_method = "tertiles",  # "tertiles", "median", "quartiles"
  min_patients_per_group = 10,  # Mínimo de pacientes por grupo
  
  # Directorio de salida
  output_dir = "Survival_Analysis_Results",
  
  # Nombres de columnas en metadata
  survival_cols = list(
    barcode = "barcode",
    vital_status = "vital_status",
    days_to_death = "days_to_death",
    days_to_last_follow_up = "days_to_last_follow_up",
    age = "age_at_diagnosis",
    stage = "ajcc_pathologic_stage",
    gender = "gender"
  )
)

################################################################################
# FUNCIONES AUXILIARES
################################################################################

# Crear directorio de salida
create_output_dirs <- function() {
  main_dir <- file.path(CONFIG$base_dir, CONFIG$output_dir)
  dir.create(main_dir, showWarnings = FALSE, recursive = TRUE)
  
  subdirs <- c("Figures", "Tables", "Reports")
  for (subdir in subdirs) {
    dir.create(file.path(main_dir, subdir), showWarnings = FALSE)
  }
  
  return(main_dir)
}

# Limpiar nombres de pacientes (normalizar separadores y eliminar sufijos)
clean_patient_id <- function(barcode) {
  # Normalizar separadores (puntos a guiones) por si acaso
  barcode <- gsub("\\.", "-", barcode)
  
  # Cortar: TCGA-BH-A0AY-01A -> TCGA-BH-A0AY
  sapply(strsplit(barcode, "-"), function(x) {
    if (length(x) >= 3) {
      paste(x[1:3], collapse = "-")
    } else {
      paste(x, collapse = "-")
    }
  })
}

################################################################################
# FUNCIONES DE CARGA DE DATOS
################################################################################

# Cargar interacciones comunes
load_common_interactions <- function() {
  cat("\n[1/5] CARGANDO INTERACCIONES COMUNES\n")
  cat(rep("=", 80), "\n", sep = "")
  
  interactions_dir <- file.path(CONFIG$base_dir, CONFIG$common_interactions_dir)
  
  all_interactions <- list()
  
  for (file in CONFIG$interaction_files) {
    filepath <- file.path(interactions_dir, file)
    
    if (!file.exists(filepath)) {
      warning(sprintf("Archivo no encontrado: %s", filepath))
      next
    }
    
    cat(sprintf("\n  Cargando: %s\n", file))
    
    data <- read_csv(filepath, show_col_types = FALSE)
    
    # Extraer proyectos del nombre del archivo
    # 02_PAR_TCGA-BRCA_vs_TCGA-LUAD.csv -> TCGA-BRCA, TCGA-LUAD
    projects <- gsub("02_PAR_|.csv", "", file)
    projects <- strsplit(projects, "_vs_")[[1]]
    
    data$Project1 <- projects[1]
    data$Project2 <- projects[2]
    data$Comparison <- paste(projects[1], "vs", projects[2])
    
    cat(sprintf("    Interacciones: %d\n", nrow(data)))
    cat(sprintf("    Genes únicos: %d\n", 
                length(unique(c(data$Gene1, data$Gene2)))))
    
    all_interactions[[file]] <- data
  }
  
  combined <- bind_rows(all_interactions)
  
  cat(sprintf("\n  ✓ Total interacciones cargadas: %d\n", nrow(combined)))
  cat(sprintf("  ✓ Comparaciones: %d\n", length(unique(combined$Comparison))))
  
  return(combined)
}

# Cargar expresión de un subdataset (MODIFICADA para usar mapeo)
load_expression_data <- function(project, subdataset) {
  project_info <- CONFIG$datasets[[project]]
  
  # Rutas de archivos
  base_path <- file.path(CONFIG$base_dir, project_info$base_path, subdataset)
  expr_file <- file.path(base_path, CONFIG$expression_filename) # counts_renombrado.csv
  map_file  <- file.path(base_path, "mapeo_muestras.csv")       # Archivo de mapeo
  
  # Verificaciones
  if (!file.exists(expr_file)) {
    warning(sprintf("Archivo de expresión no encontrado: %s", expr_file))
    return(NULL)
  }
  
  if (!file.exists(map_file)) {
    warning(sprintf("Archivo de mapeo no encontrado: %s. No se pueden recuperar los IDs originales.", map_file))
    return(NULL)
  }
  
  # 1. Cargar conteos (Generic names: Control1, Primary1...)
  expr_data <- fread(expr_file, header = TRUE)
  
  # 2. Cargar mapeo
  mapping <- fread(map_file)
  
  # 3. Crear vector de traducción: Nombre -> ID Real
  # Aseguramos que 'new_name' coincida con las columnas de expr_data
  # mapping$new_name -> mapping$barcode_se
  
  # Extraer nombres actuales de las columnas (quitando la primera que es 'gene')
  current_cols <- colnames(expr_data)[-1]
  
  # Filtrar el mapeo para quedarse solo con las muestras presentes en el archivo de counts
  # (Por seguridad, aunque deberían ser las mismas)
  map_subset <- mapping[match(current_cols, mapping$new_name), ]
  
  # Verificar si el orden se mantuvo (match lo garantiza, pero verificamos NAs)
  if(any(is.na(map_subset$barcode_se))) {
    warning("Algunas columnas de expresión no tienen correspondencia en el archivo de mapeo.")
  }
  
  # 4. Renombrar columnas en el dataframe de expresión
  # Sustituimos (Control1, Primary1...) por (TCGA-XX-XXXX...)
  colnames(expr_data)[-1] <- map_subset$barcode_se
  
  # 5. Convertir a matriz y devolver
  genes <- expr_data[[1]]
  expr_matrix <- as.matrix(expr_data[, -1])
  rownames(expr_matrix) <- genes
  
  return(expr_matrix)
}

# Cargar metadata de un proyecto
load_metadata <- function(project) {
  project_info <- CONFIG$datasets[[project]]
  
  metadata_file <- file.path(CONFIG$base_dir, project_info$metadata_file)
  
  if (!file.exists(metadata_file)) {
    stop(sprintf("Archivo de metadata no encontrado: %s", metadata_file))
  }
  
  metadata <- read_delim(metadata_file, delim = "\t", show_col_types = FALSE)
  
  # Limpiar patient IDs
  metadata$patient_clean <- clean_patient_id(metadata[[CONFIG$survival_cols$barcode]])
  
  return(metadata)
}

# Combinar expresión de todos los subdatasets de un proyecto
combine_project_expression <- function(project) {
  cat(sprintf("\n  Combinando expresión para %s...\n", project))
  
  project_info <- CONFIG$datasets[[project]]
  all_expr <- list()
  
  for (subdataset in project_info$subdatasets) {
    expr <- load_expression_data(project, subdataset)
    if (!is.null(expr)) {
      all_expr[[subdataset]] <- expr
      cat(sprintf("    ✓ %s: %d genes × %d pacientes\n", 
                  subdataset, nrow(expr), ncol(expr)))
    }
  }
  
  if (length(all_expr) == 0) {
    return(NULL)
  }
  
  # Combinar por genes comunes
  common_genes <- Reduce(intersect, lapply(all_expr, rownames))
  
  combined <- do.call(cbind, lapply(all_expr, function(x) x[common_genes, ]))
  
  cat(sprintf("    → Combinado: %d genes × %d pacientes\n", 
              nrow(combined), ncol(combined)))
  
  return(combined)
}

################################################################################
# FUNCIONES DE ANÁLISIS
################################################################################

# Preparar datos de supervivencia
prepare_survival_data <- function(metadata, expr_data, gene1, gene2, verbose = TRUE) {
  
  # Verificar que los genes existan
  if (!gene1 %in% rownames(expr_data) || !gene2 %in% rownames(expr_data)) {
    if (verbose) {
      cat(sprintf("        ⚠️  Genes no encontrados en expresión: %s, %s\n", gene1, gene2))
    }
    return(NULL)
  }
  
  # Extraer expresión
  expr1 <- expr_data[gene1, ]
  expr2 <- expr_data[gene2, ]
  
  # Limpiar IDs de pacientes de expresión
  patient_ids_expr <- clean_patient_id(colnames(expr_data))
  
  # Crear dataframe de expresión
  expr_df <- data.frame(
    patient_clean = patient_ids_expr,
    expr_gene1 = as.numeric(expr1),
    expr_gene2 = as.numeric(expr2),
    stringsAsFactors = FALSE
  )
  
  # Calcular score de co-expresión (producto de Z-scores)
  expr_df$z_gene1 <- scale(expr_df$expr_gene1)[,1]
  expr_df$z_gene2 <- scale(expr_df$expr_gene2)[,1]
  expr_df$coexpr_score <- expr_df$z_gene1 * expr_df$z_gene2
  
  # Merge con metadata
  surv_cols <- CONFIG$survival_cols
  
  metadata_subset <- metadata %>%
    select(
      patient_clean,
      vital_status = all_of(surv_cols$vital_status),
      days_to_death = all_of(surv_cols$days_to_death),
      days_to_last_follow_up = all_of(surv_cols$days_to_last_follow_up),
      any_of(c(surv_cols$age, surv_cols$stage, surv_cols$gender))
    )
  
  # DIAGNÓSTICO: Ver IDs antes del merge
  if (verbose) {
    cat(sprintf("        IDs en expresión (primeros 3): %s\n", 
                paste(head(expr_df$patient_clean, 3), collapse = ", ")))
    cat(sprintf("        IDs en metadata (primeros 3): %s\n", 
                paste(head(metadata_subset$patient_clean, 3), collapse = ", ")))
    cat(sprintf("        IDs únicos en expresión: %d\n", 
                length(unique(expr_df$patient_clean))))
    cat(sprintf("        IDs únicos en metadata: %d\n", 
                length(unique(metadata_subset$patient_clean))))
  }
  
  # Unir datos
  merged <- expr_df %>%
    inner_join(metadata_subset, by = "patient_clean")
  
  if (verbose) {
    cat(sprintf("        Pacientes después del merge: %d\n", nrow(merged)))
  }
  
  if (nrow(merged) == 0) {
    if (verbose) {
      # Buscar coincidencias parciales para diagnóstico
      common_ids <- intersect(expr_df$patient_clean, metadata_subset$patient_clean)
      cat(sprintf("        ⚠️  IDs en común: %d\n", length(common_ids)))
      
      if (length(common_ids) == 0) {
        cat("        ❌ Ningún ID coincide - Verificar formato de IDs\n")
        cat(sprintf("        Ejemplo expr: '%s'\n", expr_df$patient_clean[1]))
        cat(sprintf("        Ejemplo meta: '%s'\n", metadata_subset$patient_clean[1]))
      }
    }
    return(NULL)
  }
  
  # Calcular tiempo de supervivencia y evento
  merged <- merged %>%
    mutate(
      event = case_when(
        vital_status == "Dead" ~ 1,
        vital_status == "Alive" ~ 0,
        TRUE ~ NA_real_
      ),
      time = case_when(
        vital_status == "Dead" & !is.na(days_to_death) ~ as.numeric(days_to_death),
        vital_status == "Alive" & !is.na(days_to_last_follow_up) ~ as.numeric(days_to_last_follow_up),
        TRUE ~ NA_real_
      )
    ) %>%
    filter(!is.na(event) & !is.na(time) & time > 0)
  
  if (nrow(merged) < CONFIG$min_patients_per_group * 2) {
    if (verbose) {
      cat(sprintf("        ⚠️  Pacientes insuficientes después de filtros: %d (mínimo: %d)\n", 
                  nrow(merged), CONFIG$min_patients_per_group * 2))
    }
    return(NULL)
  }
  
  # Estratificar pacientes por score de co-expresión
  if (CONFIG$stratification_method == "tertiles") {
    merged$risk_group <- cut(
      merged$coexpr_score,
      breaks = quantile(merged$coexpr_score, probs = c(0, 1/3, 2/3, 1)),
      labels = c("Low_CoExpr", "Medium_CoExpr", "High_CoExpr"),
      include.lowest = TRUE
    )
  } else if (CONFIG$stratification_method == "median") {
    merged$risk_group <- ifelse(
      merged$coexpr_score >= median(merged$coexpr_score),
      "High_CoExpr",
      "Low_CoExpr"
    )
  }
  
  return(merged)
}

# Realizar análisis de supervivencia para una interacción
analyze_interaction_survival <- function(surv_data, gene1, gene2, project, output_dir) {
  
  if (is.null(surv_data) || nrow(surv_data) == 0) {
    return(NULL)
  }
  
  interaction_name <- paste(gene1, gene2, sep = "_")
  
  # Kaplan-Meier
  fit <- survfit(Surv(time, event) ~ risk_group, data = surv_data)
  
  # Log-rank test
  logrank <- survdiff(Surv(time, event) ~ risk_group, data = surv_data)
  logrank_pval <- 1 - pchisq(logrank$chisq, length(logrank$n) - 1)
  
  # Cox proportional hazards
  cox_model <- coxph(Surv(time, event) ~ risk_group, data = surv_data)
  cox_summary <- summary(cox_model)
  
  # Extraer hazard ratios
  if (CONFIG$stratification_method == "tertiles") {
    # Medium vs Low
    hr_medium <- cox_summary$conf.int[1, c(1, 3, 4)]  # HR, lower, upper
    p_medium <- cox_summary$coefficients[1, 5]
    
    # High vs Low
    hr_high <- cox_summary$conf.int[2, c(1, 3, 4)]
    p_high <- cox_summary$coefficients[2, 5]
    
    results <- data.frame(
      Gene1 = gene1,
      Gene2 = gene2,
      Project = project,
      N_patients = nrow(surv_data),
      N_events = sum(surv_data$event),
      LogRank_pval = logrank_pval,
      HR_Medium_vs_Low = hr_medium[1],
      HR_Medium_CI_lower = hr_medium[2],
      HR_Medium_CI_upper = hr_medium[3],
      P_Medium = p_medium,
      HR_High_vs_Low = hr_high[1],
      HR_High_CI_lower = hr_high[2],
      HR_High_CI_upper = hr_high[3],
      P_High = p_high,
      stringsAsFactors = FALSE
    )
  } else {
    # Median split
    hr <- cox_summary$conf.int[1, c(1, 3, 4)]
    p_val <- cox_summary$coefficients[1, 5]
    
    results <- data.frame(
      Gene1 = gene1,
      Gene2 = gene2,
      Project = project,
      N_patients = nrow(surv_data),
      N_events = sum(surv_data$event),
      LogRank_pval = logrank_pval,
      HR = hr[1],
      HR_CI_lower = hr[2],
      HR_CI_upper = hr[3],
      P_value = p_val,
      stringsAsFactors = FALSE
    )
  }
  
  # Generar curva de Kaplan-Meier
  plot_file <- file.path(output_dir, "Figures", 
                         sprintf("KM_%s_%s_%s.pdf", project, gene1, gene2))
  
  pdf(plot_file, width = 10, height = 8)
  
  km_plot <- ggsurvplot(
    fit,
    data = surv_data,
    pval = TRUE,
    pval.method = TRUE,
    conf.int = TRUE,
    risk.table = TRUE,
    risk.table.height = 0.25,
    ggtheme = theme_bw(),
    palette = c("#E41A1C", "#377EB8", "#4DAF4A"),
    title = sprintf("%s: %s - %s\nCo-expression Score and Survival", 
                    project, gene1, gene2),
    xlab = "Time (days)",
    ylab = "Overall Survival Probability",
    legend.title = "Co-expression Group",
    legend.labs = levels(surv_data$risk_group)
  )
  
  print(km_plot)
  dev.off()
  
  cat(sprintf("    ✓ KM plot guardado: %s\n", basename(plot_file)))
  
  return(results)
}

################################################################################
# PIPELINE PRINCIPAL
################################################################################

main <- function() {
  cat("\n")
  cat(rep("=", 80), "\n", sep = "")
  cat("ANÁLISIS DE SUPERVIVENCIA - INTERACCIONES PERDIDAS\n")
  cat(rep("=", 80), "\n", sep = "")
  
  # Crear directorios de salida
  output_dir <- create_output_dirs()
  cat(sprintf("\n📁 Directorio de salida: %s\n", output_dir))
  
  # 1. Cargar interacciones comunes
  interactions <- load_common_interactions()
  
  # 2. Cargar metadata de todos los proyectos
  cat("\n[2/5] CARGANDO METADATA\n")
  cat(rep("=", 80), "\n", sep = "")
  
  metadata_list <- list()
  for (project in names(CONFIG$datasets)) {
    cat(sprintf("\n  Cargando metadata de %s...\n", project))
    metadata_list[[project]] <- load_metadata(project)
    cat(sprintf("    ✓ Pacientes: %d\n", nrow(metadata_list[[project]])))
  }
  
  # 3. Cargar expresión de todos los proyectos
  cat("\n[3/5] CARGANDO EXPRESIÓN GÉNICA\n")
  cat(rep("=", 80), "\n", sep = "")
  
  expr_list <- list()
  for (project in names(CONFIG$datasets)) {
    expr_list[[project]] <- combine_project_expression(project)
  }
  
  # 4. Análisis de supervivencia
  cat("\n[4/5] ANÁLISIS DE SUPERVIVENCIA\n")
  cat(rep("=", 80), "\n", sep = "")
  
  all_results <- list()
  counter <- 1
  
  # Agrupar interacciones por comparación
  comparisons <- unique(interactions$Comparison)
  
  for (comp in comparisons) {
    cat(sprintf("\n  Analizando: %s\n", comp))
    cat(rep("-", 80), "\n", sep = "")
    
    comp_interactions <- interactions %>% filter(Comparison == comp)
    projects <- c(comp_interactions$Project1[1], comp_interactions$Project2[1])
    
    for (i in 1:nrow(comp_interactions)) {
      gene1 <- comp_interactions$Gene1[i]
      gene2 <- comp_interactions$Gene2[i]
      
      cat(sprintf("\n    [%d/%d] %s -- %s\n", i, nrow(comp_interactions), gene1, gene2))
      
      # Analizar en ambos proyectos
      for (project in projects) {
        cat(sprintf("      Proyecto: %s\n", project))
        
        # Preparar datos
        surv_data <- prepare_survival_data(
          metadata_list[[project]],
          expr_list[[project]],
          gene1,
          gene2
        )
        
        if (is.null(surv_data)) {
          cat(sprintf("        ⚠️  Datos insuficientes\n"))
          next
        }
        
        cat(sprintf("        Pacientes: %d | Eventos: %d\n", 
                    nrow(surv_data), sum(surv_data$event)))
        
        # Análisis de supervivencia
        result <- analyze_interaction_survival(
          surv_data, gene1, gene2, project, output_dir
        )
        
        if (!is.null(result)) {
          result$Comparison <- comp
          all_results[[counter]] <- result
          counter <- counter + 1
          
          cat(sprintf("        ✓ Log-rank p = %.4f\n", result$LogRank_pval))
        }
      }
    }
  }
  
  # 5. Consolidar y guardar resultados
  cat("\n[5/5] GUARDANDO RESULTADOS\n")
  cat(rep("=", 80), "\n", sep = "")
  
  if (length(all_results) > 0) {
    results_df <- bind_rows(all_results)
    
    # Ordenar por p-value
    results_df <- results_df %>% arrange(LogRank_pval)
    
    # Guardar tabla
    results_file <- file.path(output_dir, "Tables", "survival_results_all.csv")
    write_csv(results_df, results_file)
    cat(sprintf("\n  ✓ Resultados guardados: %s\n", results_file))
    
    # Crear resumen
    summary_file <- file.path(output_dir, "Reports", "survival_analysis_summary.txt")
    
    sink(summary_file)
    cat(rep("=", 80), "\n", sep = "")
    cat("RESUMEN DE ANÁLISIS DE SUPERVIVENCIA\n")
    cat(rep("=", 80), "\n", sep = "")
    cat(sprintf("\nFecha: %s\n", Sys.time()))
    cat(sprintf("Total de análisis: %d\n", nrow(results_df)))
    cat(sprintf("Interacciones únicas: %d\n", 
                length(unique(paste(results_df$Gene1, results_df$Gene2)))))
    cat(sprintf("Proyectos analizados: %s\n", 
                paste(unique(results_df$Project), collapse = ", ")))
    
    # Resultados significativos
    sig_results <- results_df %>% filter(LogRank_pval < 0.05)
    cat(sprintf("\nResultados significativos (p < 0.05): %d\n", nrow(sig_results)))
    
    if (nrow(sig_results) > 0) {
      cat("\nTop 10 resultados más significativos:\n")
      cat(rep("-", 80), "\n", sep = "")
      
      top10 <- head(sig_results, 10)
      for (i in 1:nrow(top10)) {
        cat(sprintf("%2d. %s -- %s [%s] | p = %.4f\n",
                    i, top10$Gene1[i], top10$Gene2[i], 
                    top10$Project[i], top10$LogRank_pval[i]))
      }
    }
    
    cat("\n", rep("=", 80), "\n", sep = "")
    sink()
    
    cat(sprintf("  ✓ Resumen guardado: %s\n", basename(summary_file)))
    
    # Mostrar resumen en consola
    cat("\n")
    cat(rep("=", 80), "\n", sep = "")
    cat("✨ ANÁLISIS COMPLETADO\n")
    cat(rep("=", 80), "\n", sep = "")
    cat(sprintf("\n  Total análisis: %d\n", nrow(results_df)))
    cat(sprintf("  Significativos (p<0.05): %d\n", nrow(sig_results)))
    cat(sprintf("\n  Resultados en: %s\n\n", output_dir))
    
  } else {
    cat("\n  ⚠️  No se generaron resultados\n\n")
  }
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
