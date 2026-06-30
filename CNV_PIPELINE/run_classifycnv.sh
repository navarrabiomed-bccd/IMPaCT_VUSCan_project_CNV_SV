#!/bin/bash
set -euo pipefail

# Verificar que se pase el argumento
if [[ $# -ne 1 ]]; then
    echo "Uso: $0 <directorio_trabajo>"
    exit 1
fi


# Rutas a las carpetas de entrada y salida
WORKDIR="$1"
FOLDER_INTERMEDIATES="${WORKDIR}CNV/INTERMEDIATES/"
RESULTS_DIR="${FOLDER_INTERMEDIATES}ClassifyCNV"

# Crear carpeta de resultados si no existe
mkdir -p "$RESULTS_DIR"

# Procesar todos los archivos BED
for BED_FILE in "${FOLDER_INTERMEDIATES}"/*.bed; do
    # Extraer el ID del archivo
    ID=$(basename "$BED_FILE" | sed 's/-.*//')
    echo "Procesando $ID ..."

    # Ejecutar ClassifyCNV desde un directorio temporal (evitar colisiones al paralelizar)
    TMPDIR_SAMPLE=$(mktemp -d)
    pushd "$TMPDIR_SAMPLE" > /dev/null
    python3 /home/ehuergo/ClassifyCNV/ClassifyCNV.py --infile "$BED_FILE" --GenomeBuild hg38
    popd > /dev/null

    # Renombrar archivo con el ID
    SCOREFILE=$(find "${TMPDIR_SAMPLE}/ClassifyCNV_results/Result"* -name "Scoresheet.txt" | head -n1)
    mv "$SCOREFILE" "$(dirname "$SCOREFILE")/${ID}_Scoresheet.txt"

    # Copiar archivo a la carpeta de la familia
    cp "$(dirname "$SCOREFILE")/${ID}_Scoresheet.txt" "${FOLDER_INTERMEDIATES}/ClassifyCNV/"

    # Mover resultados completos a la carpeta de la familia
    mv "${TMPDIR_SAMPLE}/ClassifyCNV_results/Result"* "${RESULTS_DIR}/${ID}_ClassifyCNV_result"

    # Borrar directorio temporal
    rm -rf "$TMPDIR_SAMPLE"

    echo "$ID completado"
done

echo "Todos los archivos procesados"
