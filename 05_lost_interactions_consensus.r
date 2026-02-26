# Script para encontrar interacciones comunes entre redes WGCNA de diferentes datasets
# Fecha: 2026-02-04

# Cargar librerías necesarias
library(dplyr)
library(readr)

# ============================================================================
# CONFIGURACIÓN
# ============================================================================

# Ruta base
base_path <- "/Users/mriosc/Documents/NETWORKS"

# Definir los 7 datasets
datasets <- list(
  "TCGA-LUAD/by_Tissue/Upper_lobe__lung",
  "TCGA-LUAD/by_Tissue/Lower_lobe__lung",
  "TCGA-HNSC/by_Anatomic_Region/Larynx",
  "TCGA-HNSC/by_Anatomic_Region/Oral_Cavity",
  "TCGA-HNSC/by_Anatomic_Region/Oropharynx",
  "TCGA-STAD/by_Molecular_Subtype/MSI",
  "TCGA-STAD/by_Molecular_Subtype/CIN", 
  "TCGA-BRCA/by_Subtype/Basal",
  "TCGA-BRCA/by_Subtype/Her2",
  "TCGA-BRCA/by_Subtype/LumA",
  "TCGA-BRCA/by_Subtype/LumB"
)

# Fold changes a analizar
fold_changes <- c("FC2", "FC3")

# Directorio de salida
output_dir <- file.path(base_path, "_Common_Interactions/08_03")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# FUNCIONES
# ============================================================================

# Función para crear identificador único de interacción (no direccional)
create_interaction_id <- function(gene1, gene2) {
  # Ordenar los genes alfabéticamente para que A-B sea lo mismo que B-A
  sorted_genes <- t(apply(cbind(gene1, gene2), 1, sort))
  interaction_id <- paste(sorted_genes[, 1], sorted_genes[, 2], sep = "_")
  return(interaction_id)
}

# Función para leer y procesar un archivo de interacciones
read_interactions <- function(file_path, dataset_name) {
  if (!file.exists(file_path)) {
    warning(paste("Archivo no encontrado:", file_path))
    return(NULL)
  }
  
  cat(paste("Leyendo:", file_path, "\n"))
  
  # Leer el archivo
  interactions <- read_csv(file_path, show_col_types = FALSE)
  
  # Crear ID de interacción
  interactions$Interaction_ID <- create_interaction_id(
    interactions$Gene1, 
    interactions$Gene2
  )
  
  # Añadir columna con el nombre del dataset
  interactions$Dataset <- dataset_name
  
  cat(paste("  - Interacciones encontradas:", nrow(interactions), "\n"))
  
  return(interactions)
}

# Función para encontrar interacciones comunes por tipo de cáncer
find_project_specific_interactions <- function(all_interactions) {
  
  # 1. Contar cuántos datasets tenemos realmente por cada proyecto (LUAD=2, HNSC=3, STAD=2)
  project_counts <- all_interactions %>%
    group_by(Project) %>%
    summarise(total_datasets_in_group = n_distinct(Dataset), .groups = 'drop')
  
  # 2. Agrupar interacciones por ID y por Proyecto
  project_results <- all_interactions %>%
    group_by(Interaction_ID, Project) %>%
    summarise(
      n_present = n_distinct(Dataset),
      Gene1 = first(Gene1),
      Gene2 = first(Gene2),
      Delta_Corr_Avg = mean(Delta_Correlation, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    # Unir con los conteos totales para filtrar
    left_join(project_counts, by = "Project") %>%
    # FILTRO CRÍTICO: Que esté en TODOS los datasets de SU proyecto
    filter(n_present == total_datasets_in_group)
  
  return(project_results)
}

# Función para encontrar interacciones comunes entre pares de proyectos
find_pairwise_intersections <- function(project_results) {
  projects <- unique(project_results$Project)
  pairwise_list <- list()
  
  if (length(projects) < 2) return(pairwise_list)
  
  # Generar todas las combinaciones posibles de 2 proyectos
  pares <- combn(projects, 2, simplify = FALSE)
  
  for (par in pares) {
    p1 <- par[1]
    p2 <- par[2]
    
    # Filtrar IDs que aparecen en ambos proyectos
    ids_p1 <- project_results %>% filter(Project == p1) %>% pull(Interaction_ID)
    ids_p2 <- project_results %>% filter(Project == p2) %>% pull(Interaction_ID)
    common_ids <- intersect(ids_p1, ids_p2)
    
    if (length(common_ids) > 0) {
      pair_data <- project_results %>%
        filter(Interaction_ID %in% common_ids & (Project == p1 | Project == p2)) %>%
        group_by(Interaction_ID) %>%
        summarise(
          Gene1 = first(Gene1),
          Gene2 = first(Gene2),
          Delta_Corr_P1 = Delta_Corr_Avg[Project == p1][1],  # Añadir [1]
          Delta_Corr_P2 = Delta_Corr_Avg[Project == p2][1],  # Añadir [1]
          Delta_Corr_Global = (Delta_Corr_P1 + Delta_Corr_P2) / 2,
          .groups = 'drop'
        )
      
      label <- paste0(p1, "_y_", p2)
      pairwise_list[[label]] <- pair_data
    }
  }
  return(pairwise_list)
}

# Función para encontrar el "Core Universal" (común a los 3 proyectos)
find_universal_core <- function(project_results) {
  n_proyectos <- n_distinct(project_results$Project)
  
  universal <- project_results %>%
    group_by(Interaction_ID) %>%
    filter(n_distinct(Project) == n_proyectos) %>%
    summarise(
      Gene1 = first(Gene1),
      Gene2 = first(Gene2),
      Delta_Corr_Global = mean(Delta_Corr_Avg),
      Proyectos_Incluidos = paste(sort(unique(Project)), collapse = ", "),
      .groups = 'drop'
    )
  return(universal)
}


# Función para guardar resultados y crear resumen automático
save_with_summary <- function(data, folder, filename_base, title, info_extra = "") {
  # 1. Guardar CSV
  write_csv(data, file.path(folder, paste0(filename_base, ".csv")))
  
  # 2. Crear Resumen TXT
  summary_file <- file.path(folder, paste0(filename_base, "_Resumen.txt"))
  sink(summary_file)
  cat("============================================================\n")
  cat(title, "\n")
  cat("============================================================\n\n")
  cat("Fecha:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Total de interacciones comunes:", nrow(data), "\n")
  cat(info_extra, "\n\n")
  
  if(nrow(data) > 0) {
    cat("Top 15 interacciones (por magnitud de Delta_Corr):\n")
    # Intentar usar la columna de correlación disponible
    col_corr <- intersect(names(data), c("Delta_Corr_Avg", "Delta_Correlation", "Delta_Corr_Global"))[1]
    data_top <- data %>% arrange(desc(abs(get(col_corr)))) %>% head(15)
    print(as.data.frame(data_top[, c("Gene1", "Gene2", col_corr)]))
  } else {
    cat("No se encontraron interacciones que cumplan este criterio.\n")
  }
  sink()
}



# ============================================================================
# ANÁLISIS PRINCIPAL CORREGIDO
# ============================================================================

for (fc in fold_changes) {
  cat("\nIniciando proceso para:", fc, "\n")
  
  fc_dir <- file.path(output_dir, fc)
  dir.create(fc_dir, showWarnings = FALSE)
  
  all_data <- list()
  
  # --- 1. CARGA DE DATOS ---
  for (i in seq_along(datasets)) {
    path <- datasets[[i]]
    
    # MEJORA: Extraer el nombre que empieza por "TCGA-" sin importar la profundidad
    path_parts <- strsplit(path, "/")[[1]]
    proj <- path_parts[grep("TCGA-", path_parts)][1] 
    
    df <- read_interactions(file.path(base_path, path, "2resultsWGCNA/08_03", fc, "lost_interactions.csv"), 
                            gsub("/", "_", path))
    
    if (!is.null(df)) {
      df$Project <- proj
      all_data[[i]] <- df
    }
  }
  combined <- bind_rows(all_data)
  
  # --- 2. NIVEL PROYECTOS ---
  cat("  > Analizando Proyectos...\n")
  project_results <- list()
  for (proj in unique(combined$Project)) {
    n_datasets_proj <- n_distinct(combined$Dataset[combined$Project == proj])
    
    res_proj <- combined %>%
      filter(Project == proj) %>%
      group_by(Interaction_ID) %>%
      filter(n_distinct(Dataset) == n_datasets_proj) %>%
      summarise(
        Gene1 = first(Gene1), 
        Gene2 = first(Gene2), 
        Delta_Corr_Avg = mean(Delta_Correlation, na.rm = TRUE), 
        .groups = 'drop'
      )
    
    project_results[[proj]] <- res_proj
    save_with_summary(res_proj, fc_dir, paste0("01_PROYECTO_", proj), 
                      paste("RESUMEN PROYECTO:", proj), 
                      paste("Criterio: Presente en los", n_datasets_proj, "datasets de este proyecto."))
  }
  
  # --- 3. NIVEL PARES (SOLUCIÓN AL ERROR 'Delta_Corr_Avg') ---
  cat("  > Analizando Pares...\n")
  proyectos <- names(project_results)
  if (length(proyectos) >= 2) {
    pares <- combn(proyectos, 2, simplify = FALSE)
    for (par in pares) {
      p1 <- par[1]; p2 <- par[2]
      
      common_ids <- intersect(project_results[[p1]]$Interaction_ID, 
                              project_results[[p2]]$Interaction_ID)
      
      if (length(common_ids) > 0) {
        res_par <- project_results[[p1]] %>%
          filter(Interaction_ID %in% common_ids) %>%
          inner_join(
            project_results[[p2]] %>% 
              select(Interaction_ID, Delta_Corr_P2 = Delta_Corr_Avg), 
            by = "Interaction_ID"
          ) %>%
          # Evitamos el error manteniendo los nombres claros
          mutate(Delta_Corr_P1 = Delta_Corr_Avg) %>%
          mutate(Delta_Corr_Global = (Delta_Corr_P1 + Delta_Corr_P2) / 2) %>%
          select(Interaction_ID, Gene1, Gene2, Delta_Corr_P1, Delta_Corr_P2, Delta_Corr_Global)
        
        save_with_summary(res_par, fc_dir, paste0("02_PAR_", p1, "_vs_", p2), 
                          paste("COMPARATIVA:", p1, "vs", p2),
                          "Criterio: Intersección de los consensos de ambos proyectos.")
      }
    }
  }
  
  # --- 4. NIVEL UNIVERSAL (7 datasets de 4 proyectos) ---
  cat("  > Analizando Core Universal...\n")
  n_total_ds <- n_distinct(combined$Dataset)
  res_univ <- combined %>%
    group_by(Interaction_ID) %>%
    filter(n_distinct(Dataset) == n_total_ds) %>% 
    summarise(
      Gene1 = first(Gene1), 
      Gene2 = first(Gene2), 
      Delta_Corr_Global = mean(Delta_Correlation, na.rm = TRUE), 
      .groups = 'drop'
    )
  
  save_with_summary(res_univ, fc_dir, "03_UNIVERSAL_Core_Total", 
                    "CORE UNIVERSAL: INTERSECCIÓN DE TODOS LOS DATASETS",
                    paste("Criterio: Interacción presente en los", n_total_ds, "datasets analizados."))
}


cat("\n==================================================\n")
cat("ANÁLISIS COMPLETADO\n")
cat("==================================================\n\n")
cat(paste("Los resultados se han guardado en:", output_dir, "\n\n"))

