#!/usr/bin/env Rscript
# ============================================================================
# 02_process_and_visualize.R
# ----------------------------------------------------------------------------
# Módulo de BioConvertR para búsqueda por HOMOLOGÍA DE SECUENCIA (BLAST)
# contra PlasticDB, como complemento/reemplazo del filtro por palabras clave.
#
# Este script:
#   1. Lee TODOS los resultados de BLAST (uno o más archivos, una fila por
#      hit) generados por 01_run_blast.sh en la carpeta results/
#   2. Filtra por umbrales de homología real (identidad, e-value, longitud)
#   3. Separa la descripción de PlasticDB en columnas (enzima, organismo,
#      tipo de plástico) para su análisis
#   4. (Opcional) Cruza cada muestra con una tabla de taxonomía/metadatos
#      propia, si existe en data/taxonomy.tsv (columnas: sample_id, taxon,
#      y cualquier otra métrica de calidad del ensamblaje)
#   5. Genera visualizaciones (barras, bubble plot, heatmap)
#   6. Exporta un FASTA con la secuencia EXACTA de cada gen validado
#
# R no alinea secuencias a esta escala -- eso lo hace BLAST en los scripts
# 00 y 01 -- pero es apropiado para cruzar, filtrar, visualizar y exportar
# los resultados una vez que existen.
#
# CONFIGURACIÓN: ajustar la sección "0. CONFIGURACIÓN" según el análisis.
# ============================================================================

suppressMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
  library(Biostrings)
})

# ----------------------------------------------------------------------------
# 0. CONFIGURACIÓN
# ----------------------------------------------------------------------------

DIR_RESULTS       <- "../results"       # carpeta con los *_vs_plasticdb.tsv
DIR_DATA          <- "../data"          # carpeta con PlasticDB.fasta, taxonomy.tsv (opcional)
UMBRAL_IDENTIDAD  <- 40                 # % identidad mínimo (homología real, no ruido)
UMBRAL_EVALUE     <- 1e-5               # e-value máximo (significancia estadística)
LONGITUD_MINIMA   <- 50                 # aa mínimos de alineamiento

# Si un archivo de resultados combina proteínas de varios bins/genomas en
# una sola muestra (p. ej. "genes.faa" con encabezados "Bin.005.NODE_..."),
# un solo sample_id no alcanza para diferenciarlos. Este regex, aplicado a
# qseqid, extrae ese identificador por gen individual. Poner en NULL si no
# aplica (todas las muestras son de un solo genoma/bin cada una).
BIN_ID_REGEX      <- "^Bin\\.[0-9]+"

COLUMNAS_BLAST <- c(
  "qseqid", "sseqid", "stitle", "pident", "length", "mismatch", "gapopen",
  "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qseq", "sseq"
)

# ----------------------------------------------------------------------------
# 1. LEER TODOS LOS RESULTADOS DE BLAST (una o varias muestras)
# ----------------------------------------------------------------------------

archivos_blast <- list.files(DIR_RESULTS, pattern = "_vs_plasticdb\\.tsv$", full.names = TRUE)

if (length(archivos_blast) == 0) {
  stop("No se encontraron archivos *_vs_plasticdb.tsv en ", DIR_RESULTS,
       ". Corre primero 01_run_blast.sh")
}

leer_uno <- function(ruta) {
  muestra <- str_remove(basename(ruta), "_vs_plasticdb\\.tsv$")
  read_tsv(ruta, col_names = COLUMNAS_BLAST, show_col_types = FALSE) |>
    mutate(sample_id = muestra)
}

blast_todos <- map_dfr(archivos_blast, leer_uno)
cat(sprintf(">> %d hits crudos leídos de %d archivo(s)\n",
            nrow(blast_todos), length(archivos_blast)))

# ----------------------------------------------------------------------------
# 2. FILTRAR POR UMBRALES DE HOMOLOGÍA REAL
# ----------------------------------------------------------------------------

blast_filtrado <- blast_todos |>
  filter(pident >= UMBRAL_IDENTIDAD,
         evalue <= UMBRAL_EVALUE,
         length  >= LONGITUD_MINIMA)

cat(sprintf(">> %d hits pasaron los umbrales (identidad >= %g%%, e-value <= %g, longitud >= %d aa)\n",
            nrow(blast_filtrado), UMBRAL_IDENTIDAD, UMBRAL_EVALUE, LONGITUD_MINIMA))

# ----------------------------------------------------------------------------
# 3. EXTRAER bin_id POR GEN INDIVIDUAL (si aplica)
#    Necesario cuando un archivo de resultados combina proteínas de varios
#    bins/genomas en una sola muestra (ver BIN_ID_REGEX en CONFIGURACIÓN).
# ----------------------------------------------------------------------------

if (!is.null(BIN_ID_REGEX)) {
  blast_filtrado <- blast_filtrado |>
    mutate(bin_id = str_extract(qseqid, BIN_ID_REGEX))

  if (all(is.na(blast_filtrado$bin_id))) {
    message(">> Aviso: BIN_ID_REGEX no encontró coincidencias en qseqid; ",
            "se omite la columna bin_id (revisa el patrón si esperabas bins).")
    blast_filtrado <- blast_filtrado |> select(-bin_id)
  } else {
    cat(sprintf(">> bin_id extraído para %d de %d hits.\n",
                sum(!is.na(blast_filtrado$bin_id)), nrow(blast_filtrado)))
  }
}

# ----------------------------------------------------------------------------
# 4. SEPARAR LA DESCRIPCIÓN DE PLASTICDB EN COLUMNAS
#    Formato observado en PlasticDB: codigo|enzima|organismo|tipo_plastico|accesion
# ----------------------------------------------------------------------------

blast_filtrado <- blast_filtrado |>
  separate(stitle,
           into = c("plasticdb_code", "enzyme", "reference_organism", "plastic_type", "accession"),
           sep = "\\|", extra = "merge", fill = "right", remove = FALSE)

# ----------------------------------------------------------------------------
# 5. (OPCIONAL) CRUZAR CON TAXONOMÍA/METADATOS PROPIOS
#    Si existe data/taxonomy.tsv, se cruza por "bin_id" cuando esa columna
#    exista en ambos lados (caso de bins individuales combinados en un
#    archivo), o por "sample_id" si no hay bin_id (caso de una muestra =
#    un genoma/metagenoma completo). Si no existe taxonomy.tsv, el script
#    sigue funcionando sin esa información adicional.
# ----------------------------------------------------------------------------

ruta_taxonomia <- file.path(DIR_DATA, "taxonomy.tsv")
if (file.exists(ruta_taxonomia)) {
  taxonomia <- read_tsv(ruta_taxonomia, show_col_types = FALSE)

  if ("bin_id" %in% names(taxonomia) && "bin_id" %in% names(blast_filtrado)) {
    blast_filtrado <- blast_filtrado |> left_join(taxonomia, by = "bin_id")
    cat(">> Metadatos de taxonomy.tsv incorporados (cruce por bin_id).\n")
  } else if ("sample_id" %in% names(taxonomia)) {
    blast_filtrado <- blast_filtrado |> left_join(taxonomia, by = "sample_id")
    cat(">> Metadatos de taxonomy.tsv incorporados (cruce por sample_id).\n")
  } else {
    message(">> Aviso: taxonomy.tsv no tiene columna 'bin_id' ni 'sample_id'; se omite el cruce.")
  }
}

# ----------------------------------------------------------------------------
# 6. GUARDAR TABLA RESUMEN
# ----------------------------------------------------------------------------

dir.create(file.path(DIR_RESULTS, "summary"), showWarnings = FALSE)
write_tsv(blast_filtrado, file.path(DIR_RESULTS, "summary", "validated_homology_hits.tsv"))
cat(">> Summary table saved to: results/summary/validated_homology_hits.tsv\n")

# ----------------------------------------------------------------------------
# 7. VISUALIZACIÓN 1 -- Barras: número de genes homólogos por tipo de plástico
# ----------------------------------------------------------------------------

plot_plastic_types <- blast_filtrado |>
  count(plastic_type, sort = TRUE) |>
  ggplot(aes(x = reorder(plastic_type, n), y = n)) +
  geom_col(fill = "#2166ac") +
  coord_flip() +
  labs(
    title = "Genes with homology to PlasticDB",
    subtitle = paste0("Threshold: identity >= ", UMBRAL_IDENTIDAD,
                       "%, e-value <= ", UMBRAL_EVALUE),
    x = "Plastic type (per PlasticDB)",
    y = "Number of homologous genes"
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(DIR_RESULTS, "summary", "plot_plastic_types.png"),
       plot = plot_plastic_types, width = 9, height = 7, dpi = 300)

# ----------------------------------------------------------------------------
# 8. VISUALIZACIÓN 2 -- Bubble plot: identidad y significancia por muestra
# ----------------------------------------------------------------------------

plot_bubble <- blast_filtrado |>
  ggplot(aes(x = plastic_type, y = reorder(sseqid, pident),
             size = pident, color = -log10(evalue))) +
  geom_point(alpha = 0.85) +
  scale_color_viridis_c(name = "-log10(e-value)") +
  scale_size_continuous(name = "% identity", range = c(2, 10)) +
  labs(
    title = "Sequence homology: query genes vs PlasticDB",
    subtitle = paste0(nrow(blast_filtrado), " homologous genes (identity >= ",
                       UMBRAL_IDENTIDAD, "%, e-value <= ", UMBRAL_EVALUE, ")"),
    x = "Plastic type",
    y = "Homologous protein in PlasticDB (ID)"
  ) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7))

ggsave(file.path(DIR_RESULTS, "summary", "plot_bubble.png"),
       plot = plot_bubble, width = 10, height = 8, dpi = 300)

# ----------------------------------------------------------------------------
# 9. VISUALIZACIÓN 3 -- Heatmap: identidad media por muestra x tipo de plástico
# ----------------------------------------------------------------------------

# Elegir el agrupador de filas: si hay taxonomía cruzada usamos "taxon"
# (más informativo); si no, "bin_id" cuando exista (gen individual dentro
# de una muestra combinada); si no, "sample_id" como antes. Agrupar por
# sample_id cuando hay bin_id colapsaría todos los bins en una sola fila.
columna_agrupadora <- if ("taxon" %in% names(blast_filtrado)) {
  "taxon"
} else if ("bin_id" %in% names(blast_filtrado)) {
  "bin_id"
} else {
  "sample_id"
}

matriz_heatmap <- blast_filtrado |>
  filter(!is.na(plastic_type)) |>
  group_by(across(all_of(columna_agrupadora)), plastic_type) |>
  summarise(mean_identity = mean(pident), .groups = "drop") |>
  pivot_wider(names_from = plastic_type, values_from = mean_identity, values_fill = 0) |>
  tibble::column_to_rownames(columna_agrupadora) |>
  as.matrix()

if (nrow(matriz_heatmap) >= 2 && ncol(matriz_heatmap) >= 2) {
  png(file.path(DIR_RESULTS, "summary", "heatmap_identity_by_plastic_type.png"),
      width = 1800, height = 1200, res = 200)
  pheatmap(matriz_heatmap,
           main = "Mean identity (%) by sample and homologous plastic type",
           color = colorRampPalette(c("white", "#2166ac"))(100),
           display_numbers = TRUE,
           number_format = "%.1f",
           na_col = "grey95")
  dev.off()
} else {
  message(">> Aviso: muy pocas muestras/categorías para un heatmap significativo (se omite).")
}

# ----------------------------------------------------------------------------
# 10. EXTRACCIÓN DE SECUENCIAS EXACTAS (para heterología u otro downstream)
#    Se exporta la secuencia de la proteína consultada (qseq) SIN los gaps
#    de alineamiento, lista para síntesis, clonación, o cualquier análisis
#    posterior.
# ----------------------------------------------------------------------------

columnas_encabezado <- c("qseqid", "sseqid", "enzyme", "reference_organism",
                          "plastic_type", "pident", "evalue")
if ("taxon" %in% names(blast_filtrado)) {
  columnas_encabezado <- c(columnas_encabezado, "taxon")
}

secuencias_validadas <- blast_filtrado |>
  mutate(
    clean_sequence = str_remove_all(qseq, "-"),
    fasta_header = sprintf(
      "%s|sample=%s|homolog=%s|enzyme=%s|reference_organism=%s|plastic_type=%s|identity=%.1f|evalue=%.2e%s%s",
      qseqid, sample_id, sseqid, enzyme, reference_organism, plastic_type, pident, evalue,
      if ("bin_id" %in% names(blast_filtrado)) paste0("|bin=", bin_id) else "",
      if ("taxon" %in% names(blast_filtrado)) paste0("|taxon=", taxon) else ""
    )
  )

set_fasta <- AAStringSet(secuencias_validadas$clean_sequence)
names(set_fasta) <- secuencias_validadas$fasta_header

writeXStringSet(set_fasta, filepath = file.path(DIR_RESULTS, "summary", "validated_homologous_sequences.fasta"))

cat(sprintf(">> FASTA exportado: results/summary/validated_homologous_sequences.fasta (%d secuencias)\n",
            length(set_fasta)))
cat(">> Proceso completo.\n")
