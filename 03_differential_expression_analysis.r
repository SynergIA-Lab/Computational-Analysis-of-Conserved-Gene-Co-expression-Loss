##########################################################################################################################################
# ANÁLISIS DE GENES DIFERENCIALMENTE EXPRESADOS (DEGs)
# Control vs Tumor - Después de eliminar muestras conflictivas del QC
# Dataset: Upper lobe lung
##########################################################################################################################################

# Cargar librerías
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(DESeq2)
  library(RColorBrewer)
  library(pheatmap)
  library(ggplot2)
  library(ggrepel)
  library(edgeR)
})

# Cambiar ruta directorio de trabajo
setwd("/Users/mriosc/Documents/NETWORKS/TCGA-BRCA")



##########################################################################################################################################
# FUNCIÓN AUXILIAR: LEER COUNTS
##########################################################################################################################################

read_counts_csv <- function(f) {
  df <- readr::read_csv(f, show_col_types = FALSE)
  gene_ids <- df[[1]]
  mat <- as.matrix(df[, -1, drop = FALSE])
  rownames(mat) <- gene_ids
  colnames(mat) <- make.names(colnames(mat), unique = TRUE)
  mode(mat) <- "numeric"
  mat <- round(mat)
  storage.mode(mat) <- "integer"
  mat
}

##########################################################################################################################################
# PASO 1: CARGAR DATOS Y ELIMINAR MUESTRAS CONFLICTIVAS
##########################################################################################################################################

# Cargar lista de muestras a remover (generada en el script de QC)
samples_to_remove_file <- "./by_Subtype/Basal/FINAL_samples_to_remove.txt"

if (file.exists(samples_to_remove_file)) {
  samples_to_remove <- readLines(samples_to_remove_file)
  cat(sprintf("✓ Muestras conflictivas cargadas desde QC: %d\n", length(samples_to_remove)))
} else {
  warning("⚠️ No se encuentra archivo de muestras a remover del QC. Continuando sin filtrado.")
  samples_to_remove <- character(0)
}


# --- CARGAR COUNTS PRIMARY ---
count_data_primary <- read_counts_csv("./by_Subtype/Basal/counts_primary.csv")
cat(sprintf("Counts Primary - Dimensiones originales: %d genes x %d muestras\n", 
            nrow(count_data_primary), ncol(count_data_primary)))

# Eliminar muestras conflictivas
if (length(samples_to_remove) > 0) {
  samples_to_keep_primary <- !(colnames(count_data_primary) %in% samples_to_remove)
  count_data_primary_clean <- count_data_primary[, samples_to_keep_primary]
  cat(sprintf("Counts Primary - Dimensiones limpias: %d genes x %d muestras\n", 
              nrow(count_data_primary_clean), ncol(count_data_primary_clean)))
  cat(sprintf("  → Muestras eliminadas: %d\n\n", sum(!samples_to_keep_primary)))
} else {
  count_data_primary_clean <- count_data_primary
  cat("Sin muestras para eliminar\n\n")
}

# --- CARGAR COUNTS SOLID ---
count_data_solid <- read_counts_csv("./by_Subtype/Basal/counts_control.csv")
cat(sprintf("Counts Solid - Dimensiones: %d genes x %d muestras\n", 
            nrow(count_data_solid), ncol(count_data_solid)))

# Remover muestras conflictivas de solid también
if (length(samples_to_remove) > 0) {
  samples_to_keep_solid <- !(colnames(count_data_solid) %in% samples_to_remove)
  count_data_solid <- count_data_solid[, samples_to_keep_solid]
  cat(sprintf("Counts Solid - Después de limpieza: %d genes x %d muestras\n", 
              nrow(count_data_solid), ncol(count_data_solid)))
}

##########################################################################################################################################
# PASO 2: VERIFICAR GENES COMUNES Y COMBINAR DATASETS
##########################################################################################################################################

# Verificar genes comunes
common_genes <- intersect(rownames(count_data_primary_clean), rownames(count_data_solid))

if (length(common_genes) == 0) {
  stop("ERROR: No hay genes comunes entre Primary y Solid")
}

cat(sprintf("Genes comunes entre Primary y Solid: %d\n", length(common_genes)))

if (length(common_genes) < nrow(count_data_primary_clean) || 
    length(common_genes) < nrow(count_data_solid)) {
  cat(sprintf("  ⚠ Usando solo genes en intersección\n"))
}

# Filtrar por genes comunes
count_data_primary_clean <- count_data_primary_clean[common_genes, ]
count_data_solid <- count_data_solid[common_genes, ]

# Combinar ambas matrices
counts_combined <- cbind(count_data_primary_clean, count_data_solid)

cat(sprintf("\nMatriz combinada: %d genes x %d muestras\n", 
            nrow(counts_combined), ncol(counts_combined)))

##########################################################################################################################################
# PASO 3: CREAR METADATA
##########################################################################################################################################

# Crear metadata
col_data <- data.frame(
  sample = colnames(counts_combined),
  condition = c(rep("Tumor", ncol(count_data_primary_clean)),
                rep("Control", ncol(count_data_solid))),
  stringsAsFactors = FALSE
)
rownames(col_data) <- col_data$sample

# Convertir a factor con Control como referencia
col_data$condition <- factor(col_data$condition, levels = c("Control", "Tumor"))

cat("\nMetadata creada:\n")
cat(sprintf("  - Control: %d muestras\n", sum(col_data$condition == "Control")))
cat(sprintf("  - Tumor: %d muestras\n", sum(col_data$condition == "Tumor")))

# Verificación
stopifnot(all(colnames(counts_combined) == col_data$sample))
stopifnot(all(col_data$condition %in% c("Control", "Tumor")))

##########################################################################################################################################
# PASO 4: FILTRADO POR EXPRESIÓN
##########################################################################################################################################

cat("\nFiltrado de genes con baja expresión...\n")

group <- factor(col_data$condition)
keep <- filterByExpr(counts_combined, group = group)
counts_filtered <- counts_combined[keep, ]

cat(sprintf("Genes retenidos: %d de %d (%.1f%%)\n", 
            sum(keep), length(keep), 100 * sum(keep) / length(keep)))

##########################################################################################################################################
# PASO 5: ANÁLISIS DESeq2
##########################################################################################################################################

cat("\n", rep("=", 80), "\n", sep = "")
cat("ANÁLISIS DIFERENCIAL CON DESeq2\n")
cat(rep("=", 80), "\n", sep = "")

# Crear directorio de salida
output_dir <- "./by_Subtype/Basal/2DEGs_analysis"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Crear objeto DESeq2
dds <- DESeqDataSetFromMatrix(countData = counts_filtered, 
                               colData = col_data, 
                               design = ~ condition)

cat("\nEjecutando DESeq2...\n")
dds <- DESeq(dds)


# --- 1.1 TRANSFORMACIÓN DE ESTABILIZACIÓN DE VARIANZA (VST) ---
# Justificación: Desacopla la varianza de la media para evitar sesgos en la correlación.
cat("Aplicando VST (Variance Stabilizing Transformation)...\n")
vsd <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)

# --- 1.2 FILTRADO POR VARIANZA / CV (Top 25%) ---
# Justificación: Seleccionar genes dinámicos (rewiring) y evitar sesgo de correlación por DEGs.
cat("Calculando Coeficiente de Variación (CV) para selección de genes de red...\n")

# Calcular media y desviación estándar fila por fila
means <- rowMeans(vst_mat)
sds   <- rowSds(vst_mat)
cv    <- sds / means

# Seleccionar el umbral (Top 25% CV)
quantile_threshold <- 0.75 # Top 25%
cv_cutoff <- quantile(cv, probs = quantile_threshold, na.rm = TRUE)

keep_network <- cv >= cv_cutoff & !is.na(cv)
network_matrix <- vst_mat[keep_network, ]

cat(sprintf("\n[RED] Umbral de CV (Top 25%%): %.4f\n", cv_cutoff))
cat(sprintf("[RED] Genes seleccionados para EnGNet: %d de %d originales\n", 
            nrow(network_matrix), nrow(vst_mat)))

# EXPORTAR MATRIZ PARA ENGNET
# Esta es la matriz que debes meter en tu algoritmo de correlación (Spearman/Kendall/NMI)
write.csv(network_matrix, 
          file = file.path(output_dir, "VST_Top25CV.csv"), 
          quote = FALSE)

cat("✓ Matriz para EnGNet exportada: VST_Top25CV.csv\n")


# Conteos normalizados
normalized_counts <- counts(dds, normalized = TRUE)

cat("✓ Normalización completada\n")

##########################################################################################################################################
# PASO 6: PCA
##########################################################################################################################################

cat("\nGenerando PCA...\n")

vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = name), size = 2.5, max.overlaps = 20) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA: Control vs Tumor (Lower lobe lung - CLEAN)") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")

pdf(file.path(output_dir, "PCA_Control_vs_Tumor.pdf"), width = 10, height = 7)
print(pca_plot)
dev.off()

cat("✓ PCA generado\n")


##########################################################################################################################################
# PASO 7: IDENTIFICACIÓN DE BIOMARCADORES (MÉTODO TREAT)
##########################################################################################################################################

cat("\n", rep("=", 80), "\n", sep = "")
cat("FASE 1.3: ESTADÍSTICA TREAT PARA BIOMARCADORES\n")
cat(rep("=", 80), "\n", sep = "")

# Definir umbral biológico tau (FC = 1.2 => log2(1.2) ≈ 0.263)
# Justificación: Penaliza genes ruidosos y detecta cambios modestos pero consistentes.
tau_fc <- 1.2
tau_log2 <- log2(tau_fc) 

cat(sprintf("Aplicando test TREAT con umbral log2FC > %.3f (Fold-Change > %.1f)...\n", tau_log2, tau_fc))

# Ejecutar resultados con lfcThreshold (Equivalente a TREAT en DESeq2)
# altHypothesis = "greaterAbs" comprueba si |beta| > tau
res_treat <- results(dds, 
                     contrast = c("condition", "Tumor", "Control"),
                     lfcThreshold = tau_log2,
                     altHypothesis = "greaterAbs",
                     alpha = 0.05) # FDR threshold

cat("\nResumen de resultados TREAT:\n")
summary(res_treat)

# Filtrar significativos (padj < 0.05)
# Nota: En TREAT, el p-value ya considera el umbral de FC, por lo que es más robusto.
res_treat_df <- as.data.frame(res_treat)
sig_treat <- res_treat_df %>% 
  filter(padj < 0.05) %>%
  arrange(padj)

# Separar Upregulated y Downregulated
up_genes   <- rownames(sig_treat)[sig_treat$log2FoldChange > 0]
down_genes <- rownames(sig_treat)[sig_treat$log2FoldChange < 0]

cat("\n[BIOMARCADORES] Genes identificados robustamente:\n")
cat(sprintf("  - Upregulated (Tumor > Control + umbral): %d\n", length(up_genes)))
cat(sprintf("  - Downregulated (Tumor < Control - umbral): %d\n", length(down_genes)))
cat(sprintf("  - Total Candidatos Validados: %d\n", nrow(sig_treat)))

# EXPORTAR RESULTADOS TREAT
write.csv(sig_treat, 
          file = file.path(output_dir, "biomarkers_TREAT_significant.csv"))

write.table(rownames(sig_treat),
            file = file.path(output_dir, "list_biomarkers_TREAT.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("✓ Lista de biomarcadores robustos exportada.\n")






##########################################################################################################################################
# PASO 7: RESULTADOS DESeq2 Y DEGS
##########################################################################################################################################

cat("\nExtrayendo resultados DESeq2...\n")

res <- results(dds, contrast = c("condition", "Tumor", "Control"))

cat("\nResumen de resultados:\n")
print(summary(res))

# Extraer información
fold_change <- res$log2FoldChange
adj_pval <- res$padj
gene_ids <- rownames(res)

# Identificar DEGs con criterios: |log2FC| > 1, 2 y 3 y padj < 0.05
activated_genes <- gene_ids[fold_change > 2 & adj_pval < 0.05 & !is.na(adj_pval)]
repressed_genes <- gene_ids[fold_change < -2 & adj_pval < 0.05 & !is.na(adj_pval)]
all_DEGs <- gene_ids[abs(fold_change) > 2 & adj_pval < 0.05 & !is.na(adj_pval)]

cat("\n", rep("-", 80), "\n", sep = "")
cat("GENES DIFERENCIALMENTE EXPRESADOS (|log2FC| > 3, padj < 0.05):\n")
cat(rep("-", 80), "\n", sep = "")
cat(sprintf("  - Genes activados (upregulated): %d\n", length(activated_genes)))
cat(sprintf("  - Genes reprimidos (downregulated): %d\n", length(repressed_genes)))
cat(sprintf("  - Total DEGs: %d\n", length(all_DEGs)))

##########################################################################################################################################
# PASO 8: VOLCANO PLOT
##########################################################################################################################################

cat("\nGenerando Volcano Plots...\n")

log_padj <- -log10(adj_pval)

# Volcano plot con ggplot2
volcano_data <- data.frame(
  log2FC = fold_change,
  log10padj = log_padj,
  gene = gene_ids,
  significant = ifelse(gene_ids %in% activated_genes, "Upregulated",
                      ifelse(gene_ids %in% repressed_genes, "Downregulated", "Not significant"))
)

volcano_plot <- ggplot(volcano_data, aes(x = log2FC, y = log10padj, color = significant)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c("Upregulated" = "red", 
                                 "Downregulated" = "blue", 
                                 "Not significant" = "grey")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black") +
  labs(title = "Volcano Plot: Control vs Tumor (BRCA - Basal)",
       x = "log2 Fold Change",
       y = "-log10(adjusted p-value)") +
  theme_minimal() +
  theme(legend.position = "bottom")

pdf(file.path(output_dir, "FC2/volcanoPlot_FC2.pdf"), width = 8, height = 6)
print(volcano_plot)
dev.off()

png(file.path(output_dir, "FC2/volcanoPlot_FC2.png"), width = 800, height = 600)
print(volcano_plot)
dev.off()

cat("✓ Volcano plots generados\n")

##########################################################################################################################################
# PASO 9: EXPORTAR TABLAS DE DEGs
##########################################################################################################################################

cat("\nExportando tablas de DEGs...\n")

# Separar muestras Control y Tumor
control_samples <- col_data$sample[col_data$condition == "Control"]
tumor_samples <- col_data$sample[col_data$condition == "Tumor"]

# Tabla de genes normales (Control) con DEGs
normal_DEG_table <- normalized_counts[all_DEGs, control_samples, drop = FALSE]
normal_DEG_table <- cbind(attr_name = rownames(normal_DEG_table), normal_DEG_table)

write.table(normal_DEG_table, 
            file = file.path(output_dir, "FC3/normal_DEGs_FC3.csv"), 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Tabla de genes tumor con DEGs
tumor_DEG_table <- normalized_counts[all_DEGs, tumor_samples, drop = FALSE]
tumor_DEG_table <- cbind(attr_name = rownames(tumor_DEG_table), tumor_DEG_table)

write.table(tumor_DEG_table, 
            file = file.path(output_dir, "FC3/tumor_DEGs_FC3.csv"), 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Tabla completa de resultados DESeq2
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- res_df[order(res_df$padj), ]

write.csv(res_df, 
          file = file.path(output_dir, "FC3/DESeq2_results_complete_FC3.csv"), 
          row.names = FALSE)

# Listas de DEGs (solo nombres de genes)
write.table(data.frame(gene = activated_genes),
            file = file.path(output_dir, "FC3/activated_genes_FC3.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(data.frame(gene = repressed_genes),
            file = file.path(output_dir, "FC3/repressed_genes_FC3.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(data.frame(gene = all_DEGs),
            file = file.path(output_dir, "FC3/all_DEGs_FC3.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("✓ Tablas exportadas\n")

##########################################################################################################################################
# PASO 10: HEATMAP DE TOP DEGs
##########################################################################################################################################

cat("\nGenerando heatmap de top DEGs...\n")

# Seleccionar top 50 DEGs por p-valor ajustado
top_degs <- head(all_DEGs[order(adj_pval[all_DEGs])], 50)

# Matriz VST de top DEGs
vsd_mat <- assay(vsd)
top_deg_mat <- vsd_mat[top_degs, ]

# Anotaciones
ann_col <- data.frame(
  Condition = col_data$condition,
  row.names = col_data$sample
)

pdf(file.path(output_dir, "FC3/heatmap_top50_DEGs_FC3.pdf"), width = 10, height = 12)
pheatmap(top_deg_mat,
         annotation_col = ann_col,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_colnames = FALSE,
         show_rownames = TRUE,
         fontsize_row = 6,
         main = "Top 50 DEGs (by adjusted p-value)")
dev.off()

cat("✓ Heatmap generado\n") 

##########################################################################################################################################
# RESUMEN FINAL
##########################################################################################################################################

cat("\n", rep("=", 80), "\n", sep = "")
cat("RESUMEN FINAL\n")
cat(rep("=", 80), "\n", sep = "")

cat("\nMUESTRAS:\n")
cat(sprintf("  - Control (Solid): %d\n", sum(col_data$condition == "Control")))
cat(sprintf("  - Tumor (Primary, limpio): %d\n", sum(col_data$condition == "Tumor")))
cat(sprintf("  - Total: %d\n", nrow(col_data)))

cat("\nGENES:\n")
cat(sprintf("  - Genes comunes: %d\n", length(common_genes)))
cat(sprintf("  - Genes retenidos (filtrado): %d\n", sum(keep)))

cat("\nDEGs (|log2FC| > 2, padj < 0.05):\n")
cat(sprintf("  - Upregulated: %d\n", length(activated_genes)))
cat(sprintf("  - Downregulated: %d\n", length(repressed_genes)))
cat(sprintf("  - Total: %d\n", length(all_DEGs)))

cat("\nARCHIVOS GENERADOS EN:", output_dir, "\n")
cat("  - PCA_Control_vs_Tumor.pdf\n")
cat("  - volcanoPlot.pdf / .png / _base.png\n")
cat("  - heatmap_top50_DEGs.pdf\n")
cat("  - normal_DEGs.csv\n")
cat("  - tumor_DEGs.csv\n")
cat("  - DESeq2_results_complete.csv\n")
cat("  - activated_genes.txt\n")
cat("  - repressed_genes.txt\n")
cat("  - all_DEGs.txt\n")

cat("\n✓ Análisis completo exitosamente\n")
cat(rep("=", 80), "\n", sep = "")

