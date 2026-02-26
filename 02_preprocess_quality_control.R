suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)      # si no, usar ashr en lfcShrink
  library(edgeR)
  library(limma)
  library(matrixStats)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(RColorBrewer)
  library(pheatmap)
  library(ggrepel)
  library(data.table)
  library(cluster)
  library(vegan)
  library(tibble)
})


set.seed(123)
setwd("/Users/mriosc/Documents/NETWORKS/TCGA_RNAseq/TCGA-BRCA")

dir.create("./by_Subtype/LumB/initial_qc",     showWarnings = FALSE)



# ###########################
# 1) CARGA DE DATOS (dos CSV)
# ###########################

# Espera: 1ª columna = gene_id ; resto = muestras (conteos)
read_counts_csv <- function(f) {
  df <- readr::read_csv(f, show_col_types = FALSE)
  gene_ids <- df[[1]]
  mat <- as.matrix(df[,-1, drop=FALSE])
  rownames(mat) <- gene_ids
  colnames(mat) <- make.names(colnames(mat), unique=TRUE)
  mode(mat) <- "numeric"
  mat <- round(mat)
  storage.mode(mat) <- "integer"
  mat
}

# Rutas de CSV patologia y control
path_csv <- "./by_Subtype/LumB/counts_primary.csv"  
ctrl_csv <- "./by_Subtype/LumB/counts_control.csv"  

counts_case  <- read_counts_csv(path_csv)
counts_ctrl  <- read_counts_csv(ctrl_csv)

# #####################################################################
# 2) VERIFICACION QUE LOS DOS DATASETS SON COMPATIBLES A NIVEL DE GENES
# #####################################################################
common_genes <- intersect(rownames(counts_case), rownames(counts_ctrl))

# Si los dos datasets no tiene genes en comun, se aborta la ejecucion
if (length(common_genes) == 0) stop("No hay genes comunes entre ambos CSV.")

# Si el número de genes comunes es menor que el número total de genes en alguno de los conjuntos, imprime un mensaje de advertencia 
if (length(common_genes) < nrow(counts_case) || length(common_genes) < nrow(counts_ctrl)) { 
  message("Aviso: se usarán genes en intersección: ", length(common_genes))
}

# Reordena y filtra las filas de cada matriz para quedarse solo con los genes comunes.
counts_case <- counts_case[common_genes, , drop=FALSE]
counts_ctrl <- counts_ctrl[common_genes, , drop=FALSE]

# Creo la matriz final de conteos (total de genes en comun x muestras totales (patologia y control)
counts <- cbind(counts_case, counts_ctrl)

# Creo una matriz de metatados en el que guardo que muestras son de patologia y cuales son de control.
meta <- data.frame(
  sample    = colnames(counts),
  condition = c(rep("case",  ncol(counts_case)),
                rep("control", ncol(counts_ctrl))),
  stringsAsFactors = FALSE
)
rownames(meta) <- meta$sample

# Por ultimo, verificamos que todo el proceso es correcto.
# Comprueba que cada nombre de columna en la matriz counts coincide exactamente con el nombre de muestra en meta$sample
stopifnot(all(colnames(counts) == meta$sample))

# Verifica que la columna condition contenga solo los valores permitidos: "case" o "control"
stopifnot(all(meta$condition %in% c("case","control")))

# Separa los nombres de las muestras por grupo
case_ids <- meta$sample[meta$condition=="case"]
ctrl_ids <- meta$sample[meta$condition=="control"]

# ###############################################
# 3) FILTRADO BASICO PARA DETECTAR BAJA EXPRESION
# ###############################################
# Eliminamos aquellos genes que no estan expresados de forma fiable.

# Convertimos la columna 'condition' de la metadata en un factor con los dos grupos (case/control).
# Esto se usa para que el filtrado tenga en cuenta el diseño experimental (tamaños de grupo).
group <- factor(meta$condition)  # c("case","control")

# Filtrado recomendado de genes con baja expresión usando edgeR::filterByExpr:
# - Opera en el espacio CPM (ajustado por tamaño de biblioteca) y considera los grupos.
# - Retiene genes que muestran suficiente abundancia en un número razonable de muestras de AL MENOS uno de los grupos (apto para capturar genes on/off).
# - Ajusta automáticamente el umbral según tamaños de grupo y profundidades (robusto a downsampling).
keep  <- filterByExpr(counts, group = group)  # deja min.count por defecto (CPM-aware)

# Aplicamos el filtro para quedarnos solo con los genes “expresados” de forma fiable.
# drop=FALSE asegura que el resultado siga siendo matriz aunque quede 1 gen.
counts <- counts[keep, , drop=FALSE]

# ###############################
# 4) QUALITY CONTROL (QC) DE BASE
# ###############################
# 4.1) Preparacion del DESeq2
# Aseguramos que la variable del diseño sea un factor y fijamos la referencia.
# Poner "control" como referencia hace que los LFC sean "case vs control".
meta$condition <- factor(meta$condition)
meta$condition <- relevel(meta$condition, ref = "control")

# Construimos el objeto DESeq2:
# - countData: matriz de conteos CRUDOS (tras el filtrado previo)
# - colData : metadata con la columna 'condition'
# - design  : fórmula del modelo (comparación simple ~ condition)
dds_qc <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ condition)

# Normalización por profundidad de secuenciación (size factors) - sin entrenar el modelo.
dds_qc <- estimateSizeFactors(dds_qc)
vsd_qc <- vst(dds_qc, blind = TRUE)
vmat   <- assay(vsd_qc)

#PERMANOVA 
dist_mat <- dist(t(vmat))
set.seed(123)
adonis_result <- adonis2(dist_mat ~ condition, data=meta, permutations=999)

capture.output(adonis_result, file = "./by_Subtype/LumB/permanova_result_no_preprocess.txt")



# 4.2) Comprobar tamaños de librería (RAW COUNTS)
# Con esta prueba detectamos la variabilidad de los tamaños de biblioteca.
# Lo ideal es que no sea homogeneo entre patologia vs control para justificar que se necesita realizar una normalizacion de los datos y filtrar por expresion antes de realizar un modelado.
# Ademas, nos ayudaria si una muestra tiene un perfil bajo o muy alto de expresion que pueda comportarse como outlier.

libsizes <- colSums(counts)
pdf("./by_Subtype/LumB/initial_qc/01_library_sizes_barplot.pdf", width=8, height=4)
op <- par(no.readonly=TRUE); par(mar=c(5,6.5,4,2)+0.1)
barplot(libsizes, las=2, main="Library sizes (raw counts)",
        ylab="", cex.names=0.7, col="gray")
mtext("Total counts", side=2, line=4.2)
abline(h=median(libsizes), lty=2, col="red")
par(op); dev.off()

##########################################################################################################################################

# 4.3) PCA (VST)
# Detectamos batches, outliers, clusterización por grupo.

pdf("./by_Subtype/LumB/initial_qc/02_PCA_VST.pdf", width=6, height=5)
print(DESeq2::plotPCA(vsd_qc, intgroup="condition"))
dev.off()

##########################################################################################################################################

# 4.4) Heatmap de distancias (Euclídea sobre VST)
# Se realiza una similitud global entre muestras bajo métrica euclídea despues de estabilizar la varianza.
# Basicamente complementa el PCA y la idea es obtener bloques limpios por tipo de muestra. 

ann <- data.frame(condition = meta$condition); rownames(ann) <- colnames(vmat)
pdf("./by_Subtype/LumB/initial_qc/03_sample_distance_heatmap_VST_euclid.pdf", width=7, height=6)
pheatmap(as.matrix(dist(t(vmat))),
         annotation_col = ann,
         main = "Sample distances (Euclidean on VST)",
         clustering_method = "complete",
         show_colnames = FALSE,
         show_rownames = FALSE)
dev.off()

##########################################################################################################################################

# 4.5) Heatmap de 1 - Spearman (sobre VST)
# Spearman es mas robusto a outliers y relaciones no lineales monotónicas.  Esta medida ayuda a saber si las diferencias eran debidas a outliers.
# Ignoramos escala y nos centramos en el orden de la expresión.

cor_s <- cor(vmat, method = "spearman")
dist_s <- as.dist(1 - cor_s)

pdf("./by_Subtype/LumB/initial_qc/04_sample_heatmap_VST_1minusSpearman.pdf", 
    width = 7, height = 6)

s <- pheatmap(as.matrix(1 - cor_s),
              annotation_col = ann,
              main = "Sample distances = 1 - Spearman (VST)",
              clustering_method = "complete",
              clustering_distance_cols = dist_s,
              clustering_distance_rows = dist_s,
              show_colnames = FALSE,
              show_rownames = FALSE)

dev.off()

# Guardar orden del clustering
cluster_order_s <- data.frame(
  Position = 1:length(s$tree_col$order),
  Sample = colnames(vmat)[s$tree_col$order],
  Condition = meta$condition[s$tree_col$order]
)

write.csv(cluster_order_s,
          "./by_Subtype/LumB/initial_qc/04_sample_clustering_order_Spearman.csv",
          row.names = FALSE)



# Eliminar muestra conflictiva y nuevo Heatmap
# 4.5) Heatmap de 1 - Spearman (sobre VST) - sin muestras conflictivas
# Identificar y eliminar la muestra problemática

#samples_to_remove_s <- c("Control112","Control87","Control110","Control105","Control108","Control6","Control83","Control104","Control88",
                         "Control106","Control47","Control91","Control103","Control32","Control107","Control4","Control62","Control8",
                         "Control109","Control17","Control96","Control111","Control3","Control19","Control28","Control66","Primary405",
                         "Primary509","Primary213","Primary30","Primary177","Primary472","Primary425","Primary215","Primary419","Primary333",
                         "Primary423","Primary298","Primary502","Primary363","Primary451","Primary121","Primary206","Primary518","Primary261",
                         "Primary480","Primary548","Primary194","Primary393","Primary169","Primary431","Primary412","Primary335","Primary51",
                         "Primary406","Primary317","Primary535","Primary315","Primary5","Primary340","Primary493","Primary133","Primary62",
                         "Primary63","Primary131")

samples_to_keep_s <- !(colnames(vmat) %in% samples_to_remove_s)
vmat_filtered_s <- vmat[, samples_to_keep_s]
meta_filtered_s <- meta[colnames(vmat_filtered_s), , drop=FALSE]


#PERMANOVA 
dist_mat <- dist(t(vmat_filtered_s))
set.seed(123)
adonis_result <- adonis2(dist_mat ~ condition, data=meta_filtered_s, permutations=999)

capture.output(adonis_result, file = "./by_Subtype/LumB/initial_qc/permanova_result_preprocess_spearman.txt")


ann_filtered_s <- ann[samples_to_keep_s, , drop = FALSE]

cor_s_clean <- cor(vmat_filtered_s, method = "spearman")
dist_s_clean <- as.dist(1 - cor_s_clean)

pdf("./by_Subtype/LumB/initial_qc/04_sample_heatmap_VST_1minusSpearman_CLEAN.pdf", 
    width = 7, height = 6)

s_clean <- pheatmap(as.matrix(1 - cor_s_clean),
                    annotation_col = ann_filtered_s,
                    main = "Sample distances = 1 - Spearman (VST)\n(CLEAN: removed conflictive samples)",
                    clustering_method = "complete",
                    clustering_distance_cols = dist_s_clean,
                    clustering_distance_rows = dist_s_clean,
                    show_colnames = FALSE,
                    show_rownames = FALSE)

dev.off()

# Guardar orden del clustering
cluster_order_s_clean <- data.frame(
  Position = 1:length(s_clean$tree_col$order),
  Sample = colnames(vmat_filtered_s)[s_clean$tree_col$order],
  Condition = meta$condition[colnames(vmat_filtered_s)][s_clean$tree_col$order]
)

write.csv(cluster_order_s_clean, 
          "./by_Subtype/LumB/initial_qc/04_sample_clustering_order_Spearman_CLEAN.csv",
          row.names = FALSE)


##########################################################################################################################################

# 4.6) Heatmap de 1 - Pearson (sobre VST)
# Comprobamos relaciones lineales en niveles VST. Es sensible a outliers.
# Con este heatmap se puede inspeccionar genes con valores extremos o muestras con expresiones anómalas.

cor_p <- cor(vmat, method = "pearson")
dist_p <- as.dist(1 - cor_p)

pdf("./by_Subtype/LumB/initial_qc/05_sample_heatmap_VST_1minusPearson.pdf", 
    width = 7, height = 6)

p <- pheatmap(as.matrix(1 - cor_p),
              annotation_col = ann,
              main = "Sample distances = 1 - Pearson (VST)",
              clustering_method = "complete",
              clustering_distance_cols = dist_p,
              clustering_distance_rows = dist_p,
              show_colnames = FALSE,
              show_rownames = FALSE)

dev.off()

# Guardar orden del clustering
cluster_order_p <- data.frame(
  Position = 1:length(p$tree_col$order),
  Sample = colnames(vmat)[p$tree_col$order],
  Condition = meta$condition[p$tree_col$order]
)

write.csv(cluster_order_p,
          "./by_Subtype/LumB/initial_qc/05_sample_clustering_order_Pearson.csv",
          row.names = FALSE)



# Eliminar muestra conflictiva y nuevo Heatmap
# 4.6) Heatmap de 1 - Pearson (sobre VST) - sin muestras conflictivas
# Identificar y eliminar la muestra problemática

#samples_to_remove_p <- c("Control112","Control87","Control110","Control105","Control108","Control6","Control83","Control104","Control88",
                         "Control106","Control47","Control91","Control103","Control32","Control107","Control4","Control62","Control8",
                         "Control109","Control17","Control96","Control111","Control3","Control19","Control28","Control66","Primary405",
                         "Primary509","Primary213","Primary30","Primary177","Primary472","Primary425","Primary215","Primary419","Primary333",
                         "Primary423","Primary298","Primary502","Primary363","Primary451","Primary121","Primary206","Primary518","Primary261",
                         "Primary480","Primary548","Primary194","Primary393","Primary169","Primary431","Primary412","Primary335","Primary51",
                         "Primary406","Primary317","Primary535","Primary315","Primary5","Primary340","Primary493","Primary133","Primary62",
                         "Primary63","Primary131")



samples_to_keep_p <- !(colnames(vmat) %in% samples_to_remove_p)
vmat_filtered_p <- vmat[, samples_to_keep_p]
meta_filtered_p <- meta[colnames(vmat_filtered_p), , drop=FALSE]



#PERMANOVA 
dist_mat <- dist(t(vmat_filtered_p))
set.seed(123)
adonis_result <- adonis2(dist_mat ~ condition, data=meta_filtered_p, permutations=999)

capture.output(adonis_result, file = "./by_Subtype/LumB/permanova_result_preprocess_pearson.txt")


ann_filtered_p <- ann[samples_to_keep_p, , drop = FALSE]

cor_p_clean <- cor(vmat_filtered_p, method = "pearson")
dist_p_clean <- as.dist(1 - cor_p_clean)

pdf("./by_Subtype/LumB/initial_qc/05_sample_heatmap_VST_1minusPearson_CLEAN.pdf", 
    width = 7, height = 6)

p_clean <- pheatmap(as.matrix(1 - cor_p_clean),
                    annotation_col = ann_filtered_p,
                    main = "Sample distances = 1 - Pearson (VST)\n(CLEAN: removed conflictive samples)",
                    clustering_method = "complete",
                    clustering_distance_cols = dist_p_clean,
                    clustering_distance_rows = dist_p_clean,
                    show_colnames = FALSE,
                    show_rownames = FALSE)

dev.off()

# Guardar orden del clustering
cluster_order_p_clean <- data.frame(
  Position = 1:length(p_clean$tree_col$order),
  Sample = colnames(vmat_filtered_p)[p_clean$tree_col$order],
  Condition = meta$condition[colnames(vmat_filtered_p)][p_clean$tree_col$order]
)

write.csv(cluster_order_p_clean, 
          "./by_Subtype/LumB/initial_qc/05_sample_clustering_order_Pearson_CLEAN.csv",
          row.names = FALSE)



##########################################################################################################################################

# 4.7) MDS a partir de VST + 1 - Spearman

# CON TODAS LAS MUESTRAS
cor_s <- cor(vmat, method = "spearman")
d_s <- as.dist(1 - cor_s)
mds_s <- cmdscale(d_s, k=2)  # coordenadas MDS
mds_df_s <- data.frame(Dim1=mds_s[,1], Dim2=mds_s[,2],
                       condition=meta$condition, sample=colnames(vmat))

pdf("./by_Subtype/LumB/initial_qc/06_MDS_VST_1minusSpearman.pdf", width=6, height=5)
print(
  ggplot(mds_df_s, aes(Dim1, Dim2, color=condition, label=sample)) +
    geom_point(size=3) +
    labs(title="MDS on 1 - Spearman (VST)") +
    theme_minimal()
)
dev.off()

##########################################

# MDS VST + 1- Spearman SIN MUESTRAS CONFLICTIVAS HEATMAP

# Filtrar matriz y metadatos
samples_to_keep_s <- !(colnames(vmat) %in% samples_to_remove_s)
vmat_s_filtered_s <- vmat[, samples_to_keep_s]
meta_s_filtered_s <- meta[colnames(vmat_s_filtered_s), , drop = FALSE]

# Calcular correlaciones y distancias (Spearman)
cor_s <- cor(vmat_s_filtered_s, method = "spearman")
d_s <- as.dist(1 - cor_s)

# MDS
mds_s <- cmdscale(d_s, k = 2)
mds_df_s <- data.frame(
  Dim1 = mds_s[, 1],
  Dim2 = mds_s[, 2],
  condition = meta_s_filtered_s$condition,
  sample = colnames(vmat_s_filtered_s)
)

# Plot y guardado
pdf("./by_Subtype/LumB/initial_qc/06_MDS_VST_1minusSpearman_without_samplesremoved.pdf", width = 6, height = 5)
print(
  ggplot(mds_df_s, aes(Dim1, Dim2, color = condition, label = sample)) +
    geom_point(size = 3) +
    labs(title = "MDS on 1 - Spearman (VST)\n(samples_removed)") +
    theme_minimal()
)
dev.off()

##########################################################################################################################################

# 4.8) MDS a partir de VST + 1 - Pearson

# CON TODAS LAS MUESTRAS
cor_p <- cor(vmat, method = "pearson")
d_p <- as.dist(1 - cor_p)
mds_p <- cmdscale(d_p, k=2)  # coordenadas MDS
mds_df_p <- data.frame(Dim1=mds_p[,1], Dim2=mds_p[,2],
                       condition=meta$condition, sample=colnames(vmat))

pdf("./by_Subtype/LumB/initial_qc/07_MDS_VST_1minusPearson.pdf", width=6, height=5)
print(
  ggplot(mds_df_p, aes(Dim1, Dim2, color=condition, label=sample)) +
    geom_point(size=3) +
    labs(title="MDS on 1 - Pearson (VST)") +
    theme_minimal()
)
dev.off()

##########################################

# MDS VST + 1- Pearson SIN MUESTRAS CONFLICTIVAS HEATMAP

# Filtrar matriz y metadatos
samples_to_keep_p <- !(colnames(vmat) %in% samples_to_remove_p)
vmat_p_filtered_p <- vmat[, samples_to_keep_p]
meta_p_filtered_p <- meta[colnames(vmat_p_filtered_p), , drop = FALSE]

# Calcular correlaciones y distancias (Pearson)
cor_p <- cor(vmat_p_filtered_p, method = "pearson")
d_p <- as.dist(1 - cor_p)

# MDS
mds_p <- cmdscale(d_p, k = 2)
mds_df_p <- data.frame(
  Dim1 = mds_p[, 1],
  Dim2 = mds_p[, 2],
  condition = meta_p_filtered_p$condition,
  sample = colnames(vmat_p_filtered_p)
)

# Plot y guardado
pdf("./by_Subtype/LumB/initial_qc/07_MDS_VST_1minusPearson_samplesremoved.pdf", width = 6, height = 5)
print(
  ggplot(mds_df_p, aes(Dim1, Dim2, color = condition, label = sample)) +
    geom_point(size = 3) +
    labs(title = "MDS on 1 - Pearson (VST)\n(samples_removed)") +
    theme_minimal()
)
dev.off()


#########################################################################################################################################

# GENERAR ARCHIVO CON TODAS LAS MUESTRAS ELIMINADAS
# Combinar muestras de samples_to_remove_s y samples_to_remove_p

# Inicializar variables si no existen
if (!exists("samples_to_remove_s")) {
    samples_to_remove_s <- character(0)
  }

if (!exists("samples_to_remove_p")) {
    samples_to_remove_p <- character(0)
  }

# Combinar todas las muestras únicas
all_samples_to_remove <- unique(c(samples_to_remove_s, samples_to_remove_p))

# Guardar en archivo txt (vacío si no hay muestras)
if (length(all_samples_to_remove) > 0) {
    writeLines(all_samples_to_remove, 
                               "./by_Subtype/LumB/FINAL_samples_to_remove.txt")
    cat("\nArchivo generado: FINAL_samples_to_remove.txt\n")
    cat("Total de muestras eliminadas:", length(all_samples_to_remove), "\n")
    cat("Muestras:", paste(all_samples_to_remove, collapse = ", "), "\n")
  } else {
      # Crear archivo vacío
      file.create("./by_Subtype/LumB/FINAL_samples_to_remove.txt")
      cat("\nArchivo generado: FINAL_samples_to_remove.txt (vacío)\n")
      cat("No hay muestras para eliminar.\n")
    }

