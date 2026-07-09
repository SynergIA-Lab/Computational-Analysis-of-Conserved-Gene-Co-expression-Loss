#!/usr/bin/env python3
"""
deg_intersection_analysis.py

Análisis de intersección de DEGs (FC2) dentro de cada cohorte y entre las
cuatro cohortes de cáncer, para responder al Major Comment 2 del Reviewer 2
(Network Modeling Analysis in Health Informatics and Bioinformatics).

Uso:
  1. Revisa la sección CONFIGURACIÓN de abajo, en particular SUBCOHORT_FOLDER
     (solo el nombre de carpeta de BRCA/Basal está confirmado; el resto son
     mi mejor suposición y puede que necesiten ajuste).
  2. Ejecuta primero solo el chequeo de rutas:
         python3 deg_intersection_analysis.py --check
     Esto NO calcula nada, solo te dice qué carpetas/archivos encuentra y
     cuáles no, para que corrijas nombres antes de lanzar el análisis completo.
  3. Cuando todas las rutas den OK, ejecuta el análisis completo:
         python3 deg_intersection_analysis.py
  4. Los resultados se guardan en OUTPUT_DIR (por defecto ./deg_analysis_output/).
"""

import itertools
import os
import sys

# ---------------------------------------------------------------------------
# CONFIGURACIÓN
# ---------------------------------------------------------------------------

BASE_DIR = (
    "/Users/mriosc/Documents/Doctorado/papers/01_Ríos-Cadenas et al./"
    "Computational Analysis of Conserved Gene Co-expression Loss Reveals "
    "Prognostic Regulatory Network Disruption Across Solid Tumors/"
    "1_Data_and_Results"
)

OUTPUT_DIR = "./deg_analysis_output"

# Carpeta de proyecto TCGA por cohorte (confirmado por el ejemplo BRCA).
COHORT_TCGA_FOLDER = {
    "BRCA": "TCGA-BRCA",
    "LUAD": "TCGA-LUAD",
    "HNSC": "TCGA-HNSC",
    "STAD": "TCGA-STAD",
}

# Carpeta de estratificación por cohorte (confirmado por el usuario).
COHORT_BY_FOLDER = {
    "BRCA": "by_Subtype",
    "HNSC": "by_Anatomic_Region",
    "LUAD": "by_Tissue",
    "STAD": "by_Molecular_Subtype",
}

# Subcarpeta común dentro de cada subtipo (confirmado por el ejemplo BRCA/Basal).
DEG_SUBPATH = os.path.join("1DEGs_analysis", "FC2")

# label -> (cohorte, nombre exacto de la subcarpeta de subtipo)
# Confirmado a partir de 05_lost_interactions_consensus.r (líneas 16-27) y
# 07_survival_analysis.r (líneas 43-54). Ojo: "Upper_lobe__lung" y
# "Lower_lobe__lung" llevan doble guion bajo -- no es un error de tipeo.
SUBCOHORT_FOLDER = {
    "Basal":      ("BRCA", "Basal"),
    "Her2":       ("BRCA", "Her2"),
    "LumA":       ("BRCA", "LumA"),
    "LumB":       ("BRCA", "LumB"),
    "LUAD_Upper": ("LUAD", "Upper_lobe__lung"),
    "LUAD_Lower": ("LUAD", "Lower_lobe__lung"),
    "Larynx":     ("HNSC", "Larynx"),
    "OralCavity": ("HNSC", "Oral_Cavity"),
    "Oropharynx": ("HNSC", "Oropharynx"),
    "CIN":        ("STAD", "CIN"),
    "MSI":        ("STAD", "MSI"),
}

# Qué subcohortes forman cada cohorte (para el consenso intra-cohorte).
COHORT_MEMBERS = {
    "BRCA": ["Basal", "Her2", "LumA", "LumB"],
    "LUAD": ["LUAD_Upper", "LUAD_Lower"],
    "HNSC": ["Larynx", "OralCavity", "Oropharynx"],
    "STAD": ["CIN", "MSI"],
}

# El nombre "all_DEGs_FC2.txt" (deducido del script 03) no coincide con la
# realidad en disco. Se deja en None para que el script detecte automáticamente
# el .txt dentro de cada carpeta FC2 y, si no lo encuentra, te muestre qué hay
# realmente ahí dentro.
FIXED_FILENAME = None

# Recuentos reportados en la Tabla 4 del manuscrito (main_R1.tex), para
# verificación automática.
EXPECTED_COUNTS_TABLE4 = {
    "Basal": 2640, "Her2": 2617, "LumA": 2025, "LumB": 2790,
    "LUAD_Upper": 2848, "LUAD_Lower": 2345,
    "Larynx": 1867, "OralCavity": 2661, "Oropharynx": 2614,
    "CIN": 3551, "MSI": 2607,
}

# Genes únicos que aparecen en las 18 interacciones conservadas (Tabla 5,
# tab:lost_interactions_detailed). Actualizar a mano si la tabla cambia.
CONSERVED_INTERACTION_GENES = {
    "CST1", "MMP11", "COL10A1", "IGHV5-78", "TNFSF11", "IGKV3-11", "JSRP1",
    "FCRL4", "IGHV1OR15-2", "IGHV4-4", "B4GALNT4", "TMEM59L", "IGHV5-51",
    "IGKV1-9", "IGHV3-20", "IGHV4-55", "HOXC11", "ZIC2", "TNXB", "WNT2",
    "GRIN2D", "AC010789.1", "AL162413.1", "TRPM2-AS", "ADAMTS2",
}


# ---------------------------------------------------------------------------
# CONSTRUCCIÓN DE RUTAS
# ---------------------------------------------------------------------------

def build_fc2_dir(label):
    cohort, subfolder = SUBCOHORT_FOLDER[label]
    return os.path.join(
        BASE_DIR,
        COHORT_TCGA_FOLDER[cohort],
        COHORT_BY_FOLDER[cohort],
        subfolder,
        DEG_SUBPATH,
    )


def find_deg_file(fc2_dir):
    """Devuelve (ruta_completa, error). Si FIXED_FILENAME está fijado, lo usa
    directamente; si no, busca el único .txt dentro de la carpeta. Si algo
    falla, el mensaje de error incluye el contenido real de la carpeta."""
    if not os.path.isdir(fc2_dir):
        return None, f"carpeta no encontrada: {fc2_dir}"

    try:
        contents = sorted(os.listdir(fc2_dir))
    except OSError as e:
        return None, f"no se pudo leer la carpeta {fc2_dir}: {e}"

    if FIXED_FILENAME is not None:
        if FIXED_FILENAME in contents:
            return os.path.join(fc2_dir, FIXED_FILENAME), None
        return None, (f"'{FIXED_FILENAME}' no existe dentro de: {fc2_dir}\n"
                       f"          Contenido real de la carpeta: {contents}")

    txt_files = sorted(f for f in contents if f.lower().endswith(".txt"))
    if len(txt_files) == 0:
        return None, (f"no hay ningún .txt dentro de: {fc2_dir}\n"
                       f"          Contenido real de la carpeta: {contents}")
    if len(txt_files) == 1:
        return os.path.join(fc2_dir, txt_files[0]), None
    # más de un .txt: preferir uno que contenga 'deg'
    deg_like = [f for f in txt_files if "deg" in f.lower()]
    if len(deg_like) == 1:
        return os.path.join(fc2_dir, deg_like[0]), None
    return None, (f"hay varios .txt y no está claro cuál usar en: {fc2_dir}\n"
                   f"          Candidatos: {txt_files}\n"
                   f"          -> fija FIXED_FILENAME con el nombre correcto arriba")


def check_all_paths():
    """Comprueba las 11 rutas y devuelve dict label -> ruta_archivo (o None)."""
    print("=== Chequeo de rutas ===")
    resolved = {}
    any_missing = False
    for label in SUBCOHORT_FOLDER:
        fc2_dir = build_fc2_dir(label)
        filepath, error = find_deg_file(fc2_dir)
        if filepath is None:
            print(f"  [FALTA] {label:12s} -> {error}")
            any_missing = True
            resolved[label] = None
        else:
            print(f"  [OK]    {label:12s} -> {filepath}")
            resolved[label] = filepath
    if any_missing:
        print("\nHay rutas sin resolver. Ajusta SUBCOHORT_FOLDER (o FIXED_FILENAME) "
              "arriba y vuelve a ejecutar con --check.")
    else:
        print("\nTodas las rutas se han resuelto correctamente.")
    return resolved, any_missing


# ---------------------------------------------------------------------------
# ANÁLISIS (igual que en la v1)
# ---------------------------------------------------------------------------

def load_gene_set(filepath):
    with open(filepath, encoding="utf-8") as fh:
        lines = [line.strip() for line in fh if line.strip()]
    gene_set = set(lines)
    if len(lines) != len(gene_set):
        n_dup = len(lines) - len(gene_set)
        print(f"  [AVISO] {filepath}: {n_dup} línea(s) duplicada(s) ignorada(s).")
    return gene_set


def check_against_table4(data, expected):
    print("\n=== Verificación contra la Tabla 4 del manuscrito ===")
    all_ok = True
    for label, expected_n in expected.items():
        actual_n = len(data[label])
        status = "OK" if actual_n == expected_n else "DIFERENTE"
        if actual_n != expected_n:
            all_ok = False
        print(f"  {label:15s} paper={expected_n:5d}  archivo={actual_n:5d}  [{status}]")
    if all_ok:
        print("  Todos los recuentos coinciden con la Tabla 4.")
    else:
        print("  Hay diferencias respecto a la Tabla 4 (ver arriba). "
              "El análisis continúa con los archivos proporcionados.")


def build_cohort_consensus(data, cohort_members):
    cohorts = {}
    for cohort, members in cohort_members.items():
        sets_to_intersect = [data[m] for m in members]
        cohorts[cohort] = set.intersection(*sets_to_intersect)
    return cohorts


def write_gene_list(filepath, genes):
    with open(filepath, "w", encoding="utf-8") as fh:
        for g in sorted(genes):
            fh.write(g + "\n")


def run_analysis(resolved_paths):
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("\n=== Cargando los 11 archivos de DEGs (FC2) ===")
    data = {}
    for label, path in resolved_paths.items():
        genes = load_gene_set(path)
        data[label] = genes
        print(f"  {label:15s} {len(genes):5d} genes")

    check_against_table4(data, EXPECTED_COUNTS_TABLE4)

    print("\n=== DEG-set de consenso por cohorte (intersección de subtipos) ===")
    cohorts = build_cohort_consensus(data, COHORT_MEMBERS)
    with open(os.path.join(OUTPUT_DIR, "cohort_consensus_counts.csv"), "w") as fh:
        fh.write("Cohort,N_genes\n")
        for cohort, genes in cohorts.items():
            print(f"  {cohort:6s} {len(genes):5d} genes")
            fh.write(f"{cohort},{len(genes)}\n")
            write_gene_list(os.path.join(OUTPUT_DIR, f"cohort_consensus_{cohort}.txt"), genes)

    print("\n=== Intersecciones por pares entre cohortes ===")
    with open(os.path.join(OUTPUT_DIR, "pairwise_intersection_counts.csv"), "w") as fh:
        fh.write("Cohort_A,Cohort_B,N_shared_genes\n")
        for a, b in itertools.combinations(cohorts.keys(), 2):
            inter = cohorts[a] & cohorts[b]
            print(f"  {a}-{b:6s} {len(inter):5d} genes")
            fh.write(f"{a},{b},{len(inter)}\n")
            write_gene_list(os.path.join(OUTPUT_DIR, f"pairwise_{a}_{b}.txt"), inter)

    print("\n=== Intersección a 4 vías (BRCA ∩ LUAD ∩ HNSC ∩ STAD) ===")
    four_way = set.intersection(*cohorts.values())
    print(f"  {len(four_way)} genes comunes a las 4 cohortes")
    write_gene_list(os.path.join(OUTPUT_DIR, "four_way_intersection.txt"), four_way)

    print("\n=== Cruce con los genes de las 18 interacciones conservadas (Tabla 5) ===")
    overlap = sorted(four_way & CONSERVED_INTERACTION_GENES)
    print(f"  {len(overlap)} de los genes de la Tabla 5 están en la intersección a 4 vías:")
    for g in overlap:
        print(f"    - {g}")
    with open(os.path.join(OUTPUT_DIR, "four_way_x_conserved_interactions.txt"), "w") as fh:
        fh.write("\n".join(overlap) + "\n")

    print(f"\nTodos los resultados se han guardado en: {os.path.abspath(OUTPUT_DIR)}")


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

def main():
    resolved, any_missing = check_all_paths()

    if "--check" in sys.argv:
        return  # solo chequeo de rutas, no calcular nada

    if any_missing:
        sys.exit("\nERROR: hay rutas sin resolver (ver arriba). "
                  "Corrige la configuración o ejecuta primero con --check.")

    run_analysis(resolved)


if __name__ == "__main__":
    main()