#!/usr/bin/env Rscript
# ============================================================================
# correlation_distribution_and_robustness.R
#
# Responde al Major Comment 3 del Reviewer 2 (Network Modeling Analysis in
# Health Informatics and Bioinformatics):
#   (i)  Reporta la distribucion de valores de correlacion (r_N, r_T) por
#        dataset a partir de las matrices ya calculadas por el pipeline
#        (correlation_matrix_control.rds / correlation_matrix_tumor.rds,
#        FC2), sin recalcular nada.
#   (ii) Repite la logica de consenso cross-cancer (misma que
#        05_lost_interactions_consensus.r) pero usando criterios
#        alternativos top-1% / top-5% / top-10% basados en rango de
#        Delta_r = r_N - r_T en vez del doble umbral fijo (0.7 / 0.3), y
#        compara cada resultado con las 18 interacciones actuales via
#        indice de Jaccard (mismo formato que la Tabla S2 del Online
#        Resource 1).
#
# USO:
#   1. Revisa CONFIGURACION abajo (BASE_DIR, SUBCOHORT_FOLDER).
#   2. Ejecuta primero solo el chequeo de rutas:
#        Rscript correlation_distribution_and_robustness.R --check
#   3. Cuando todas las rutas esten OK, ejecuta el analisis completo:
#        Rscript correlation_distribution_and_robustness.R
#   4. Resultados en ./correlation_analysis_output/
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
})

# ----------------------------------------------------------------------------
# CONFIGURACION
# ----------------------------------------------------------------------------

BASE_DIR <- "/Users/mriosc/Documents/Doctorado/papers/01_Ríos-Cadenas et al./Computational Analysis of Conserved Gene Co-expression Loss Reveals Prognostic Regulatory Network Disruption Across Solid Tumors/1_Data_and_Results"

OUTPUT_DIR <- "./correlation_analysis_output"

COHORT_TCGA_FOLDER <- c(BRCA = "TCGA-BRCA", LUAD = "TCGA-LUAD",
                         HNSC = "TCGA-HNSC", STAD = "TCGA-STAD")

COHORT_BY_FOLDER <- c(BRCA = "by_Subtype", HNSC = "by_Anatomic_Region",
                       LUAD = "by_Tissue", STAD = "by_Molecular_Subtype")

# label -> c(cohorte, subcarpeta exacta) -- mismos valores confirmados que en
# el script de Major Comment 2 (05_lost_interactions_consensus.r / 07_survival_analysis.r)
SUBCOHORT_FOLDER <- list(
  Basal      = c("BRCA", "Basal"),
  Her2       = c("BRCA", "Her2"),
  LumA       = c("BRCA", "LumA"),
  LumB       = c("BRCA", "LumB"),
  LUAD_Upper = c("LUAD", "Upper_lobe__lung"),
  LUAD_Lower = c("LUAD", "Lower_lobe__lung"),
  Larynx     = c("HNSC", "Larynx"),
  OralCavity = c("HNSC", "Oral_Cavity"),
  Oropharynx = c("HNSC", "Oropharynx"),
  CIN        = c("STAD", "CIN"),
  MSI        = c("STAD", "MSI")
)

FC <- "FC2"
TOP_PERCENTS <- c(0.01, 0.05, 0.10)   # 1%, 5%, 10% -- tres puntos de robustez
PLOT_SAMPLE_PER_DATASET <- 50000  # submuestreo solo para el grafico de densidad
SEED <- 1234

# Confirmado: las 11 subcohortes tienen las mismas 5 carpetas de run
# (06_02, 06_03, 07_02, 07_03, 08_03) bajo 2resultsWGCNA, y para BRCA Basal
# se verifico que las 5 son identicas byte a byte (misma dimension, mismos
# genes, mismos valores de correlacion). Se usa la mas reciente como
# preferida; si para algun dataset esa carpeta no existiera, el script cae
# automaticamente al modo de busqueda/deteccion de ambiguedad de antes.
PREFERRED_RUN_ID <- "08_03"

# Las 18 interacciones conservadas de la Tabla 5 del manuscrito (referencia).
# Interaction_ID = genes ordenados alfabeticamente y unidos por "_", igual
# que create_interaction_id() en 05_lost_interactions_consensus.r.
make_interaction_id <- function(g1, g2) {
  ord <- ifelse(g1 < g2, paste(g1, g2, sep = "_"), paste(g2, g1, sep = "_"))
  ord
}

REFERENCE_18 <- tibble::tribble(
  ~Gene1, ~Gene2,
  "CST1", "MMP11",
  "COL10A1", "CST1",
  "IGHV5-78", "TNFSF11",
  "IGKV3-11", "JSRP1",
  "FCRL4", "TNFSF11",
  "IGHV1OR15-2", "IGHV4-4",
  "B4GALNT4", "TMEM59L",
  "IGHV5-51", "JSRP1",
  "IGKV1-9", "JSRP1",
  "IGHV1OR15-2", "IGHV3-20",
  "IGHV4-55", "JSRP1",
  "HOXC11", "ZIC2",
  "TNXB", "WNT2",
  "GRIN2D", "ZIC2",
  "AC010789.1", "AL162413.1",
  "AL162413.1", "GRIN2D",
  "GRIN2D", "TRPM2-AS",
  "ADAMTS2", "TNXB"
) %>% mutate(Interaction_ID = make_interaction_id(Gene1, Gene2))

# ----------------------------------------------------------------------------
# DESCUBRIMIENTO DE RUTAS
# ----------------------------------------------------------------------------

subdataset_dir <- function(label) {
  info <- SUBCOHORT_FOLDER[[label]]
  cohort <- info[1]; subfolder <- info[2]
  file.path(BASE_DIR, COHORT_TCGA_FOLDER[[cohort]], COHORT_BY_FOLDER[[cohort]], subfolder)
}

# Busca correlation_matrix_control.rds / correlation_matrix_tumor.rds.
# Primero intenta la ruta directa con PREFERRED_RUN_ID (rapido y
# deterministico); si esa carpeta concreta no existe para este dataset, cae
# al modo de busqueda recursiva bajo cualquier */FC2/* y avisa si hay
# ambiguedad entre varios runs.
find_fc2_matrices <- function(label) {
  base <- file.path(subdataset_dir(label), "2resultsWGCNA")
  if (!dir.exists(base)) {
    return(list(ok = FALSE, error = paste("carpeta no encontrada:", base)))
  }

  # --- Intento 1: ruta directa con el run preferido ---
  direct_control <- file.path(base, PREFERRED_RUN_ID, FC, "correlation_matrix_control.rds")
  direct_tumor <- file.path(base, PREFERRED_RUN_ID, FC, "correlation_matrix_tumor.rds")
  if (file.exists(direct_control) && file.exists(direct_tumor)) {
    return(list(ok = TRUE, control = direct_control, tumor = direct_tumor))
  }

  # --- Intento 2: busqueda recursiva (fallback si falta ese run concreto) ---
  all_control <- list.files(base, pattern = "^correlation_matrix_control\\.rds$",
                             recursive = TRUE, full.names = TRUE)
  fc2_control <- all_control[grepl(paste0("/", FC, "/"), all_control, fixed = TRUE)]

  if (length(fc2_control) == 0) {
    return(list(ok = FALSE, error = paste0("no se encontro correlation_matrix_control.rds bajo ningun */",
                                            FC, "/* dentro de: ", base)))
  }
  if (length(fc2_control) > 1) {
    return(list(ok = FALSE, error = paste0(
      "no existe el run preferido ('", PREFERRED_RUN_ID, "') para este dataset, y hay VARIOS ",
      "correlation_matrix_control.rds bajo */", FC, "/* -- ambiguo:\n            ",
      paste(fc2_control, collapse = "\n            "),
      "\n          -> edita SUBCOHORT_FOLDER o PREFERRED_RUN_ID para '", label, "'")))
  }

  control_path <- fc2_control[1]
  tumor_path <- file.path(dirname(control_path), "correlation_matrix_tumor.rds")
  if (!file.exists(tumor_path)) {
    return(list(ok = FALSE, error = paste("existe el control pero falta el tumor esperado en:", tumor_path)))
  }

  list(ok = TRUE, control = control_path, tumor = tumor_path)
}

check_all_paths <- function() {
  cat("=== Chequeo de rutas (matrices de correlacion FC2) ===\n")
  resolved <- list()
  any_missing <- FALSE
  for (label in names(SUBCOHORT_FOLDER)) {
    res <- find_fc2_matrices(label)
    if (!res$ok) {
      cat(sprintf("  [FALTA] %-12s -> %s\n", label, res$error))
      any_missing <- TRUE
      resolved[[label]] <- NULL
    } else {
      cat(sprintf("  [OK]    %-12s -> %s\n", label, res$control))
      resolved[[label]] <- res
    }
  }
  if (any_missing) {
    cat("\nHay rutas sin resolver. Ajusta SUBCOHORT_FOLDER / BASE_DIR y repite --check.\n")
  } else {
    cat("\nTodas las rutas se han resuelto correctamente.\n")
  }
  list(resolved = resolved, any_missing = any_missing)
}

# ----------------------------------------------------------------------------
# EXTRACCION DE PARES A PARTIR DE LAS MATRICES
# ----------------------------------------------------------------------------

# Devuelve un data.frame con Gene1, Gene2 (orden alfabetico), r_N, r_T,
# Delta_Correlation, para todos los pares del triangulo superior.
extract_pairs <- function(control_path, tumor_path) {
  cor_control <- readRDS(control_path)
  cor_tumor <- readRDS(tumor_path)

  common_genes <- intersect(rownames(cor_control), rownames(cor_tumor))
  if (length(common_genes) < ncol(cor_control) || length(common_genes) < ncol(cor_tumor)) {
    cat(sprintf("    [aviso] genes no coincidentes entre control (%d) y tumor (%d); usando interseccion (%d)\n",
                ncol(cor_control), ncol(cor_tumor), length(common_genes)))
  }
  cor_control <- cor_control[common_genes, common_genes]
  cor_tumor   <- cor_tumor[common_genes, common_genes]

  idx <- which(upper.tri(cor_control), arr.ind = TRUE)
  g1_raw <- common_genes[idx[, 1]]
  g2_raw <- common_genes[idx[, 2]]
  r_N <- cor_control[idx]
  r_T <- cor_tumor[idx]

  swap <- g1_raw > g2_raw
  Gene1 <- ifelse(swap, g2_raw, g1_raw)
  Gene2 <- ifelse(swap, g1_raw, g2_raw)

  out <- tibble(Gene1 = Gene1, Gene2 = Gene2, r_N = r_N, r_T = r_T,
                Delta_Correlation = r_N - r_T)
  out <- out %>% filter(!is.na(r_N), !is.na(r_T))
  out$Interaction_ID <- make_interaction_id(out$Gene1, out$Gene2)
  out
}

# ----------------------------------------------------------------------------
# CONSENSO CROSS-CANCER (misma logica que 05_lost_interactions_consensus.r,
# pasos 2 y 3: consenso a nivel de proyecto + interseccion por pares)
# ----------------------------------------------------------------------------

run_consensus <- function(per_dataset_top1pct) {
  combined <- bind_rows(per_dataset_top1pct, .id = "Dataset")
  combined$Project <- sapply(combined$Dataset, function(d) SUBCOHORT_FOLDER[[d]][1])

  # Nivel proyecto: interaccion debe estar en TODOS los subdatasets del proyecto
  project_results <- list()
  for (proj in unique(combined$Project)) {
    n_datasets_proj <- n_distinct(combined$Dataset[combined$Project == proj])
    res_proj <- combined %>%
      filter(Project == proj) %>%
      group_by(Interaction_ID) %>%
      filter(n_distinct(Dataset) == n_datasets_proj) %>%
      summarise(Gene1 = first(Gene1), Gene2 = first(Gene2),
                Delta_Corr_Avg = mean(Delta_Correlation, na.rm = TRUE),
                .groups = "drop")
    project_results[[proj]] <- res_proj
  }

  # Nivel pares: interseccion de los consensos de cada dos proyectos
  proyectos <- names(project_results)
  pairwise_all <- list()
  if (length(proyectos) >= 2) {
    pares <- combn(proyectos, 2, simplify = FALSE)
    for (par in pares) {
      p1 <- par[1]; p2 <- par[2]
      common_ids <- intersect(project_results[[p1]]$Interaction_ID,
                               project_results[[p2]]$Interaction_ID)
      if (length(common_ids) > 0) {
        res_par <- project_results[[p1]] %>%
          filter(Interaction_ID %in% common_ids) %>%
          select(Interaction_ID, Gene1, Gene2)
        pairwise_all[[paste(p1, p2, sep = "_vs_")]] <- res_par
      }
    }
  }

  bind_rows(pairwise_all, .id = "Cancer_Pair")
}

# ----------------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------------

main <- function() {
  check <- check_all_paths()
  if ("--check" %in% commandArgs(trailingOnly = TRUE)) return(invisible())
  if (check$any_missing) stop("Hay rutas sin resolver (ver arriba). Corrige la configuracion primero.")

  dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
  set.seed(SEED)

  summary_rows <- list()
  plot_sample_rows <- list()
  top_by_pct <- setNames(vector("list", length(TOP_PERCENTS)), as.character(TOP_PERCENTS))
  for (pct in TOP_PERCENTS) top_by_pct[[as.character(pct)]] <- list()

  for (label in names(SUBCOHORT_FOLDER)) {
    cat(sprintf("\n=== Procesando %s ===\n", label))
    paths <- check$resolved[[label]]
    pairs <- extract_pairs(paths$control, paths$tumor)
    n_pairs <- nrow(pairs)
    cat(sprintf("  Pares extraidos: %s\n", format(n_pairs, big.mark = ",")))

    probs <- c(0.25, 0.5, 0.75, 0.90, 0.95, 0.99, 0.999)
    qN <- quantile(pairs$r_N, probs = probs, na.rm = TRUE)
    qT <- quantile(pairs$r_T, probs = probs, na.rm = TRUE)

    summary_rows[[label]] <- tibble(
      Dataset = label, N_pairs = n_pairs,
      Mean_rN = mean(pairs$r_N), SD_rN = sd(pairs$r_N),
      Mean_rT = mean(pairs$r_T), SD_rT = sd(pairs$r_T),
      P25_rN = qN[1], P50_rN = qN[2], P75_rN = qN[3], P90_rN = qN[4], P95_rN = qN[5], P99_rN = qN[6], P999_rN = qN[7],
      P25_rT = qT[1], P50_rT = qT[2], P75_rT = qT[3], P90_rT = qT[4], P95_rT = qT[5], P99_rT = qT[6], P999_rT = qT[7]
    )

    # submuestreo solo para graficar (no afecta a los percentiles de arriba)
    n_sample <- min(PLOT_SAMPLE_PER_DATASET, n_pairs)
    idx_sample <- sample.int(n_pairs, n_sample)
    plot_sample_rows[[label]] <- bind_rows(
      tibble(Dataset = label, Condition = "Normal", r = pairs$r_N[idx_sample]),
      tibble(Dataset = label, Condition = "Tumor",  r = pairs$r_T[idx_sample])
    )

    # Ordenar una sola vez por Delta_r descendente; los top-X% son prefijos
    # crecientes de la misma lista ordenada (mas eficiente que reordenar
    # para cada porcentaje).
    sorted_pairs <- pairs %>% arrange(desc(Delta_Correlation))
    for (pct in TOP_PERCENTS) {
      key <- as.character(pct)
      n_top <- ceiling(pct * n_pairs)
      top_by_pct[[key]][[label]] <- sorted_pairs %>%
        slice_head(n = n_top) %>%
        select(Gene1, Gene2, Interaction_ID, Delta_Correlation)
      cat(sprintf("  Top %.0f%%: %s pares (umbral Delta_r >= %.3f)\n",
                  pct * 100, format(n_top, big.mark = ","),
                  min(top_by_pct[[key]][[label]]$Delta_Correlation)))
    }

    rm(pairs, sorted_pairs); gc(verbose = FALSE)
  }

  # ---- (i) Tabla de distribuciones ----
  summary_table <- bind_rows(summary_rows)
  write_csv(summary_table, file.path(OUTPUT_DIR, "correlation_distribution_summary.csv"))
  cat("\n✓ Tabla de distribuciones guardada: correlation_distribution_summary.csv\n")

  # ---- (i) Grafico de densidad ----
  plot_data <- bind_rows(plot_sample_rows)
  p <- ggplot(plot_data, aes(x = r, fill = Condition, color = Condition)) +
    geom_density(alpha = 0.35) +
    facet_wrap(~Dataset, ncol = 3) +
    theme_minimal(base_size = 10) +
    labs(x = "Pearson correlation (r)", y = "Density",
         title = "Distribution of pairwise correlation values per dataset (FC2)",
         subtitle = sprintf("Normal vs tumor tissue (subsampled to %s pairs/condition for plotting)",
                             format(PLOT_SAMPLE_PER_DATASET, big.mark = ",")))
  ggsave(file.path(OUTPUT_DIR, "correlation_distributions.pdf"), p, width = 11, height = 9)
  ggsave(file.path(OUTPUT_DIR, "correlation_distributions.png"), p, width = 11, height = 9, dpi = 200)
  cat("✓ Grafico guardado: correlation_distributions.pdf / .png\n")

  # ---- (ii) Robustness check: top-1% / top-5% / top-10% ----
  cat("\n=== Robustness check: consenso cross-cancer bajo criterios rank-based (Delta_r) ===\n")
  ref_ids <- unique(REFERENCE_18$Interaction_ID)
  robustness_rows <- list()

  for (pct in TOP_PERCENTS) {
    key <- as.character(pct)
    cat(sprintf("\n--- Top %.0f%% ---\n", pct * 100))
    new_conserved <- run_consensus(top_by_pct[[key]])
    write_csv(new_conserved, file.path(OUTPUT_DIR, sprintf("top%.0fpct_conserved_interactions.csv", pct * 100)))

    new_ids <- unique(new_conserved$Interaction_ID)
    shared <- intersect(new_ids, ref_ids)
    union_ids <- union(new_ids, ref_ids)
    jaccard <- length(shared) / length(union_ids)

    cat(sprintf("  Interacciones conservadas: %d\n", length(new_ids)))
    cat(sprintf("  Compartidas con las 18 (0.7/0.3): %d\n", length(shared)))
    cat(sprintf("  Union: %d\n", length(union_ids)))
    cat(sprintf("  Jaccard: %.3f\n", jaccard))
    if (length(shared) > 0) {
      cat("  Interacciones de referencia recuperadas:\n")
      for (id in shared) {
        row <- REFERENCE_18[REFERENCE_18$Interaction_ID == id, ]
        cat(sprintf("    - %s -- %s\n", row$Gene1[1], row$Gene2[1]))
      }
    }

    robustness_rows[[key]] <- tibble(
      Configuration = sprintf("Top %.0f%% (rank-based, Delta_r)", pct * 100),
      n_total = length(new_ids),
      Shared_with_0.7_0.3 = length(shared),
      Union = length(union_ids),
      Jaccard = round(jaccard, 3)
    )
  }

  robustness_summary <- bind_rows(robustness_rows)
  write_csv(robustness_summary, file.path(OUTPUT_DIR, "toppct_vs_reference_summary.csv"))
  cat("\n✓ Resumen robustness check guardado: toppct_vs_reference_summary.csv\n")
  print(as.data.frame(robustness_summary))

  cat(sprintf("\nTodos los resultados se han guardado en: %s\n", normalizePath(OUTPUT_DIR)))
}

main()