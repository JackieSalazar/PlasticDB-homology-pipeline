#!/usr/bin/env bash
#
# 01_run_blast.sh
# --------------------------------------------------------------------------
# Corre blastp de un archivo de proteínas (.faa/.fasta) contra la base local
# de PlasticDB creada por 00_build_plasticdb_database.sh.
#
# Reporta explícitamente las columnas qseq/sseq: la secuencia EXACTA que
# hizo match en el alineamiento, para permitir su extracción posterior
# (p. ej. para diseño de primers, síntesis génica, o estudios de expresión
# heteróloga).
#
# IMPORTANTE (validado empíricamente): max_target_seqs limita cuántos de
# los MEJORES resultados se MUESTRAN por cada proteína consultada -- no
# limita contra cuántas secuencias de PlasticDB se compara (eso siempre es
# contra la base completa). Este script detecta automáticamente el número
# total de secuencias en PlasticDB y lo usa como techo, para capturar TODOS
# los hits que superen el umbral de e-value, no solo un subconjunto arbitrario.
#
# Uso:
#   ./01_run_blast.sh ruta/a/proteinas.faa nombre_muestra
#
# Ejemplo:
#   ./01_run_blast.sh /ruta/genoma/proteins.faa mi_genoma
# --------------------------------------------------------------------------

set -euo pipefail

if [ $# -lt 1 ]; then
  echo "Uso: $0 ruta/a/proteinas.faa [nombre_muestra]"
  exit 1
fi

QUERY_FAA="$1"
SAMPLE_NAME="${2:-$(basename "$QUERY_FAA" | sed 's/\.[^.]*$//')}"

DIR_SCRIPT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DIR_BASE="$(cd "$DIR_SCRIPT/.." && pwd)"
DIR_DATA="$DIR_BASE/data"
DIR_RESULTS="$DIR_BASE/results"
mkdir -p "$DIR_RESULTS"

BLASTDB="$DIR_DATA/plasticdb_blastdb"
if [ ! -f "${BLASTDB}.pin" ] && [ ! -f "${BLASTDB}.pdb" ]; then
  echo "ERROR: no se encontró la base de BLAST en $BLASTDB"
  echo "Corre primero: 00_build_plasticdb_database.sh"
  exit 1
fi

# Detectar automáticamente cuántas secuencias tiene PlasticDB, para usarlo
# como techo de max_target_seqs (evita perder hits significativos por un
# límite arbitrario).
if [ -f "$DIR_DATA/.plasticdb_n_sequences.txt" ]; then
  MAX_HITS=$(cat "$DIR_DATA/.plasticdb_n_sequences.txt")
else
  MAX_HITS=$(grep -c "^>" "$DIR_DATA/PlasticDB.fasta" 2>/dev/null || echo 500)
fi

OUTPUT_TSV="$DIR_RESULTS/${SAMPLE_NAME}_vs_plasticdb.tsv"

echo ">> Corriendo blastp: $SAMPLE_NAME vs PlasticDB ($MAX_HITS secuencias)..."

# Columnas de salida (en este orden):
#   qseqid   = ID de la proteína consultada
#   sseqid   = ID de la proteína homóloga en PlasticDB
#   stitle   = descripción/organismo/tipo de plástico de esa proteína
#   pident   = % identidad
#   length   = longitud del alineamiento
#   mismatch, gapopen = calidad del alineamiento
#   qstart/qend, sstart/send = coordenadas del alineamiento
#   evalue, bitscore = significancia estadística
#   qseq, sseq = secuencia exacta alineada (consulta y PlasticDB)
blastp \
  -query "$QUERY_FAA" \
  -db "$BLASTDB" \
  -evalue 1e-5 \
  -max_target_seqs "$MAX_HITS" \
  -outfmt "6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq" \
  -num_threads 4 \
  -out "$OUTPUT_TSV"

N_HITS=$(wc -l < "$OUTPUT_TSV" || echo 0)
echo ">> Listo: $N_HITS hits guardados en $OUTPUT_TSV"
