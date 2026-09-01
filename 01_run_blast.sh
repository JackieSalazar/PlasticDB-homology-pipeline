#!/bin/bash
# ==============================================================================
# 01_run_blast.sh
# Busqueda de homologia (blastp) del proteoma de entrada contra la base local
# de PlasticDB.
#
# ACTUALIZADO: se agrega qcovs (cobertura de la consulta) al -outfmt, no
# incluida en la version original de este script. Esta columna es necesaria
# para el analisis de sensibilidad (umbral estricto) implementado en
# 02_procesar_y_visualizar.R.
#
# Columnas del -outfmt (16, en este orden):
#   qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart
#   send evalue bitscore qcovs qseq sseq
# ==============================================================================

set -euo pipefail

QUERY="$1"        # archivo FASTA de proteinas de entrada (.faa)
DB="$2"           # ruta a la base local de PlasticDB (sin extension)
OUT="$3"          # archivo de salida .tsv

EVALUE="1e-5"
MAX_TARGET_SEQS=329   # total de secuencias en PlasticDB al momento de la
                       # descarga; verificar con: grep -c "^>" PlasticDB.fasta

blastp \
  -query "$QUERY" \
  -db "$DB" \
  -evalue "$EVALUE" \
  -max_target_seqs "$MAX_TARGET_SEQS" \
  -outfmt "6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qseq sseq" \
  -out "$OUT"

echo "blastp completado. Resultado: $OUT"
