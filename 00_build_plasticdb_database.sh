#!/usr/bin/env bash
#
# 00_build_plasticdb_database.sh
# --------------------------------------------------------------------------
# Descarga las secuencias de proteínas curadas de PlasticDB (plasticdb.org)
# y construye una base de datos local de BLAST (proteína), para permitir
# búsquedas de HOMOLOGÍA DE SECUENCIA contra genomas
# o metagenomas ya anotados (p. ej. proteínas .faa de Prokka, DRAM, PGAP,
# o exportaciones de KBase).
#
# Requiere: curl (o wget), BLAST+ (makeblastdb)
#
# Uso:
#   ./00_build_plasticdb_database.sh
# --------------------------------------------------------------------------

set -euo pipefail

DIR_SCRIPT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$DIR_SCRIPT/../data"
DIR_DATA="$(cd "$DIR_SCRIPT/../data" && pwd)"
cd "$DIR_DATA"

echo ">> Descargando FASTA de proteínas de PlasticDB..."
URL_FASTA="https://plasticdb.org/static/PlasticDB.fasta"
URL_METADATA="https://plasticdb.org/static/degraders_list.tsv"

if command -v curl >/dev/null 2>&1; then
  curl -L -o PlasticDB.fasta "$URL_FASTA"
  curl -L -o degraders_list.tsv "$URL_METADATA" || \
    echo "!! No se pudo descargar la metadata automáticamente; descárgala manualmente desde https://plasticdb.org/downloaddata"
else
  wget -O PlasticDB.fasta "$URL_FASTA"
  wget -O degraders_list.tsv "$URL_METADATA" || \
    echo "!! No se pudo descargar la metadata automáticamente; descárgala manualmente desde https://plasticdb.org/downloaddata"
fi

if [ ! -s PlasticDB.fasta ]; then
  echo "ERROR: PlasticDB.fasta está vacío o no se descargó."
  echo "Descárgalo manualmente desde https://plasticdb.org/downloaddata"
  echo "y colócalo en: $DIR_DATA/PlasticDB.fasta"
  exit 1
fi

N_SECUENCIAS=$(grep -c "^>" PlasticDB.fasta || true)
echo ">> PlasticDB.fasta descargado: $N_SECUENCIAS secuencias de proteína."
echo "$N_SECUENCIAS" > .plasticdb_n_sequences.txt

echo ">> Construyendo base de datos local de BLAST (proteína)..."
# NOTA (validado empíricamente): NO usar -parse_seqids aquí.
# Algunos encabezados de PlasticDB superan los 50 caracteres permitidos
# por esa opción, lo cual hace fallar la construcción de la base de datos
# con el error "the local id is too long". -parse_seqids no es necesario
# para blastp/tblastn estándar, solo para búsquedas por ID exacto.
makeblastdb -in PlasticDB.fasta -dbtype prot -out plasticdb_blastdb

echo ">> Listo. Base de datos creada en: $DIR_DATA/plasticdb_blastdb.*"
echo ">> Total de secuencias en la base: $N_SECUENCIAS"
echo ">> Siguiente paso: correr 01_run_blast.sh apuntando a tus proteínas (.faa)"
