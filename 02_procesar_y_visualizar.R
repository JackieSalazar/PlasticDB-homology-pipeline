# ==============================================================================
# Pipeline de homología de secuencia vs PlasticDB (v2)
# Hermetia illucens gut metagenome - deteccion de genes candidatos de
# biodegradacion de plastico
#
# CAMBIOS RESPECTO A LA VERSION ANTERIOR (v1), en respuesta a revision por
# pares (CITIS 2026):
#   1. Se agrego la columna qcovs (cobertura de la consulta) al -outfmt de
#      blastp, no solicitada en la corrida original.
#   2. Se aplica DEDUPLICACION por gen unico (qseqid): cuando una misma
#      proteina de consulta produce mas de un hit contra PlasticDB, se
#      conserva unicamente el mejor (menor e-value, y en empate, mayor
#      bitscore). Anteriormente se contaban filas de alineamiento, no genes
#      unicos, lo que sobreestimaba el numero real de candidatos.
#   3. Se verifica el SOLAPAMIENTO entre los candidatos identificados a
#      nivel de metagenoma completo y los identificados a nivel de bins
#      individuales (mismo metagenoma), para evitar sumar el mismo gen dos
#      veces al reportar el total combinado.
#   4. La extraccion de secuencias en FASTA ahora recupera la PROTEINA
#      COMPLETA desde el archivo de entrada original (usando qseqid), en
#      lugar de reconstruirla a partir del fragmento alineado (qseq), que
#      solo cubre la porcion homologa y no la proteina completa.
#   5. Se agrego un analisis de sensibilidad con umbrales mas estrictos
#      (60% identidad / 70% cobertura / e-value <=1e-30) y un umbral
#      intermedio de respaldo (50% identidad / e-value <=1e-20), para
#      contextualizar la confianza de los candidatos reportados bajo el
#      umbral permisivo principal (40% identidad / 50 aa / e-value <=1e-5).
#
# Requiere como entrada los archivos generados por blastp con el siguiente
# -outfmt (16 columnas, en este orden):
#   qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart
#   send evalue bitscore qcovs qseq sseq
#
# Comando blastp de referencia (ver comandos_ejecutados_v2.txt para el
# historial completo):
#   blastp -query <proteinas.fasta> -db <plasticdb_blastdb> -evalue 1e-5 \
#     -max_target_seqs 329 \
#     -outfmt "6 qseqid sseqid stitle pident length mismatch gapopen qstart
#              qend sstart send evalue bitscore qcovs qseq sseq" \
#     -out <resultado.tsv>
#
# Entorno: R 4.6.1, BLAST+ v2.17.0
# ==============================================================================

library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)
library(ggplot2)
library(Biostrings)

# ------------------------------------------------------------------------
# 0. Configuracion: ajustar esta ruta a la carpeta de trabajo local
# ------------------------------------------------------------------------
setwd("C:/Users/Jacque/Documents/BioConvertR_PlasticDB2")

col_names_v2 <- c("qseqid", "sseqid", "stitle", "pident", "length", "mismatch",
                   "gapopen", "qstart", "qend", "sstart", "send", "evalue",
                   "bitscore", "qcovs", "qseq", "sseq")

# ------------------------------------------------------------------------
# 1. Carga de datos crudos (metagenoma completo y bins individuales)
# ------------------------------------------------------------------------
meta_raw <- read_tsv("results/metagenoma_vs_plasticdb_v2.tsv",
                      col_names = col_names_v2, col_types = cols(.default = "c"))

bins_raw <- read_tsv("results/bins_vs_plasticdb_v2.tsv",
                      col_names = col_names_v2, col_types = cols(.default = "c"))

meta <- meta_raw %>%
  mutate(pident = as.numeric(pident), length = as.numeric(length),
         evalue = as.numeric(evalue), bitscore = as.numeric(bitscore),
         qcovs = as.numeric(qcovs))

bins <- bins_raw %>%
  mutate(pident = as.numeric(pident), length = as.numeric(length),
         evalue = as.numeric(evalue), bitscore = as.numeric(bitscore),
         qcovs = as.numeric(qcovs))

cat("Metagenoma - alineamientos crudos:", nrow(meta), "\n")
cat("Bins - alineamientos crudos:", nrow(bins), "\n\n")

# ------------------------------------------------------------------------
# 2. Analisis de sensibilidad: umbral ESTRICTO y umbral INTERMEDIO
#    (se reportan junto con el umbral principal como contexto de confianza)
# ------------------------------------------------------------------------
meta_estricto <- meta %>% filter(pident >= 60, qcovs >= 70, evalue <= 1e-30)
bins_estricto <- bins %>% filter(pident >= 60, qcovs >= 70, evalue <= 1e-30)

cat("=== UMBRAL ESTRICTO (identidad>=60%, qcovs>=70%, e-value<=1e-30) ===\n")
cat("Metagenoma:", nrow(meta_estricto), "| Bins:", nrow(bins_estricto), "\n\n")

meta_intermedio <- meta %>% filter(pident >= 50, evalue <= 1e-20)
bins_intermedio <- bins %>% filter(pident >= 50, evalue <= 1e-20)

cat("=== UMBRAL INTERMEDIO DE RESPALDO (identidad>=50%, e-value<=1e-20) ===\n")
cat("Metagenoma:", nrow(meta_intermedio), "| Bins:", nrow(bins_intermedio), "\n\n")

# ------------------------------------------------------------------------
# 3. Umbral PERMISIVO principal (identidad>=40%, longitud>=50aa, e-value<=1e-5)
#    + DEDUPLICACION por mejor hit (unico gen por qseqid)
# ------------------------------------------------------------------------
meta_permisivo <- meta %>% filter(pident >= 40, length >= 50, evalue <= 1e-5)
bins_permisivo <- bins %>% filter(pident >= 40, length >= 50, evalue <= 1e-5)

meta_dedup <- meta_permisivo %>%
  arrange(qseqid, evalue, desc(bitscore)) %>%
  distinct(qseqid, .keep_all = TRUE)

bins_dedup <- bins_permisivo %>%
  arrange(qseqid, evalue, desc(bitscore)) %>%
  distinct(qseqid, .keep_all = TRUE)

cat("=== DESPUES DE DEDUPLICAR (1 mejor hit por gen) ===\n")
cat("Metagenoma - genes unicos:", nrow(meta_dedup), "\n")
cat("Bins - genes unicos:", nrow(bins_dedup), "\n\n")

# ------------------------------------------------------------------------
# 4. Verificacion de solapamiento metagenoma vs bins (mismo gen, dos
#    conjuntos distintos). El core_id se construye quitando el prefijo
#    especifico de cada archivo (unico dato que difiere entre ambos).
# ------------------------------------------------------------------------
meta_dedup <- meta_dedup %>%
  mutate(core_id = str_remove(qseqid, "^BSFmerged_metaSPAdes\\.Assembly_"))

bins_dedup <- bins_dedup %>%
  mutate(core_id = str_remove(qseqid, "^Bin\\.[0-9]+\\.fasta_assembly_"))

genes_solapados <- intersect(meta_dedup$core_id, bins_dedup$core_id)

cat("=== SOLAPAMIENTO metagenoma vs bins ===\n")
cat("Genes de bins ya presentes en el metagenoma completo:",
    length(genes_solapados), "de", nrow(bins_dedup), "\n")
cat("Genes de bins EXCLUSIVOS (no detectados a nivel de metagenoma):",
    nrow(bins_dedup) - length(genes_solapados), "\n\n")

total_unico <- length(union(meta_dedup$core_id, bins_dedup$core_id))
cat("=== TOTAL DE GENES CANDIDATOS UNICOS (metagenoma + bins, sin duplicar) ===\n")
cat(total_unico, "\n\n")

# ------------------------------------------------------------------------
# 5. Anotacion funcional (parseo de stitle) y cruce con taxonomia de bins
# ------------------------------------------------------------------------
meta_dedup <- meta_dedup %>%
  separate(stitle, into = c("codigo_plasticdb", "enzima", "organismo",
                             "tipo_plastico", "accesion"),
           sep = "\\|", extra = "merge", fill = "right", remove = FALSE)

bins_dedup <- bins_dedup %>%
  separate(stitle, into = c("codigo_plasticdb", "enzima", "organismo",
                             "tipo_plastico", "accesion"),
           sep = "\\|", extra = "merge", fill = "right", remove = FALSE) %>%
  mutate(bin_id = str_extract(qseqid, "^Bin\\.[0-9]+"))

# Tabla de taxonomia y calidad de los 13 bins de alta calidad (CheckM + GTDB-Tk)
taxonomia_bins <- tribble(
  ~bin_id_original, ~taxon, ~completeness, ~contamination, ~size_bp, ~gc,
  "bin.001", "g__Dysgonomonas",       98.91, 0.00, 3507342, 0.33555,
  "bin.002", "f__Lachnospiraceae",    99.35, 2.22, 4933684, 0.50164,
  "bin.003", "g__Vagococcus_E",       92.40, 1.52, 2754790, 0.36482,
  "bin.005", "g__Scrofimicrobium",    99.33, 2.70, 3404153, 0.47299,
  "bin.006", "g__CHH4-2",             98.42, 0.63, 3866424, 0.40742,
  "bin.007", "f__Enterobacteriaceae", 93.93, 1.13, 1809183, 0.36655,
  "bin.008", "g__Dysgonomonas",      100.00, 0.27, 3059482, 0.34031,
  "bin.010", "g__Dendrosporobacter",  99.28, 1.90, 3154279, 0.46163,
  "bin.011", "g__Orbus",              99.44, 0.00, 3115942, 0.35738,
  "bin.012", "g__UBA5962",            99.33, 0.00, 3484036, 0.36920,
  "bin.013", "g__Entomomonas",        94.64, 0.94, 3267125, 0.37342,
  "bin.014", "f__Lachnospiraceae",    98.42, 2.53, 3893934, 0.39396,
  "bin.015", "g__Campylobacter_B",    90.15, 0.49, 1724306, 0.28780
) %>%
  mutate(bin_id = str_to_title(bin_id_original))

bin_lookup <- bins_dedup %>%
  select(core_id, bin_id) %>%
  left_join(taxonomia_bins, by = "bin_id") %>%
  select(core_id, bin_id, taxon)

tabla_final <- meta_dedup %>%
  left_join(bin_lookup, by = "core_id") %>%
  dplyr::select(qseqid, core_id, bin_id, taxon, sseqid, accesion, enzima,
                 organismo, tipo_plastico, pident, length, bitscore, evalue, qcovs)

cat("Genes con bin/taxon asignado:", sum(!is.na(tabla_final$bin_id)),
    "de", nrow(tabla_final), "\n")
cat("Genes SIN bin asignado (solo metagenoma completo):",
    sum(is.na(tabla_final$bin_id)), "\n\n")

# ------------------------------------------------------------------------
# 6. Extraccion de la secuencia PROTEICA COMPLETA (no el fragmento
#    alineado) desde el FASTA original de entrada, usando qseqid
# ------------------------------------------------------------------------
proteinas_originales <- readAAStringSet("genomas/metagenoma_kbase/metagenoma_proteins.fasta")
names(proteinas_originales) <- word(names(proteinas_originales), 1)

tabla_final <- tabla_final %>%
  mutate(secuencia_completa = as.character(proteinas_originales[qseqid]))

cat("Secuencias no encontradas en el FASTA original (deberia ser 0):",
    sum(is.na(tabla_final$secuencia_completa)), "\n\n")

encabezados <- with(tabla_final, sprintf(
  "%s|bin=%s|taxon=%s|homologo=%s|enzima=%s|plastico=%s|identidad=%.1f|evalue=%.2e|cobertura=%.1f",
  qseqid, ifelse(is.na(bin_id), "NA", bin_id), ifelse(is.na(taxon), "NA", taxon),
  sseqid, enzima, tipo_plastico, pident, evalue, qcovs
))

set_fasta_final <- AAStringSet(tabla_final$secuencia_completa)
names(set_fasta_final) <- encabezados
writeXStringSet(set_fasta_final, filepath = "results/secuencias_homologas_32genes_COMPLETAS.fasta")

# ------------------------------------------------------------------------
# 7. Exportar tabla final (candidatos, para Tabla 2 del articulo)
# ------------------------------------------------------------------------
tabla_final_export <- tabla_final %>%
  dplyr::select(qseqid, bin_id, taxon, sseqid, accesion, enzima, tipo_plastico,
                 pident, length, bitscore, evalue, qcovs) %>%
  dplyr::rename(gen_id = qseqid, bin = bin_id, taxon_GTDBTk = taxon,
                plasticdb_id = sseqid, plasticdb_accesion = accesion,
                enzima_homologa = enzima, identidad_pct = pident,
                longitud_alineamiento = length, cobertura_pct = qcovs) %>%
  dplyr::arrange(desc(identidad_pct))

write_csv(tabla_final_export, "results/tabla_final_genes_candidatos.csv")

# ------------------------------------------------------------------------
# 8. Bubble plot (metagenoma completo, genes deduplicados)
# ------------------------------------------------------------------------
datos <- read_csv("results/tabla_final_genes_candidatos.csv", show_col_types = FALSE) %>%
  mutate(log_evalue = -log10(evalue),
         etiqueta_y = paste0(plasticdb_id, " | ", enzima_homologa))

plot_bubble <- ggplot(datos, aes(x = tipo_plastico, y = etiqueta_y)) +
  geom_point(aes(size = identidad_pct, color = log_evalue), alpha = 0.85) +
  scale_color_viridis_c(name = "-log10(e-value)", option = "viridis") +
  scale_size_continuous(name = "% identity", range = c(3, 10)) +
  labs(
    title = "Sequence homology: metagenome genes vs PlasticDB",
    subtitle = sprintf("n = %d unique genes | identity >= 40%%, length >= 50 aa,\ne-value <= 1e-05 (deduplicated)", nrow(datos)),
    x = "Plastic type",
    y = "Homologous protein in PlasticDB"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold"))

ggsave("results/bubble_plot_genes_metagenoma.png", plot = plot_bubble,
       width = 12, height = 8, dpi = 300)

cat("Pipeline completo. Archivos generados en results/:\n",
    "- secuencias_homologas_32genes_COMPLETAS.fasta\n",
    "- tabla_final_genes_candidatos.csv\n",
    "- bubble_plot_genes_metagenoma.png\n")
