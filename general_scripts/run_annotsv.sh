#!/bin/bash
set -euo pipefail

# Script de anotacion de variantes somaticas con AnnotSV.
# Tipo de variante según nodo:
# - NASERTIC: CNVs (TO y PON) y SVs (TO)
# - FPGMX: CNVs y SVs
# - CNAG: No necesario (archivos originales anotados)

# Verificar argumentos de entrada
if [[ $# -ne 2 ]]; then
    echo "Uso: $0 <NODE> <working_directory>"
    exit 1
fi

# Parametros
NODE="$1"
WORKDIR="$2"

# Ruta al ejecutable de AnnotSV
ANNOTSV_BIN="/home/ehuergo/AnnotSV/bin/AnnotSV"

if [[ ! -x "$ANNOTSV_BIN" ]]; then
    echo "[ERROR] No existe o no es ejecutable: $ANNOTSV_BIN"
    exit 1
fi

# Función para filtrar PASS con bcftools -> ejecutar AnnotSV -> comprimir
run_annotsv_tool() {
    local VCF_FILE="$1"  # Archivo VCF de entrada
    local OUTPUT_TSV_GZ="$2"  # Ruta esperada del archivo comprimido (.tsv.gz)
    local TMP_DIR="$3"  # Carpeta temporal para los intermedios
    local ID PREPARED_VCF ANNOT_TSV  # Variables locales auxiliares

    # Nombre del archivo sin la extensión (.vcf o .vcf.gz)
    ID=$(basename "$VCF_FILE" | sed -E 's/\.vcf(\.gz)?$//')
    # Nombre del VCF temporal filtrado por PASS y comprimido
    PREPARED_VCF="${TMP_DIR}/${ID}.prepared.vcf.gz"
    # Nombre del TSV temporal generado por AnnotSV (sin comprimir)
    ANNOT_TSV="${TMP_DIR}/${ID}.tmp.tsv"

    # Preparar VCF filtrando por PASS
    bcftools view -f PASS "$VCF_FILE" -Oz -o "$PREPARED_VCF"

    # Indexar VCF preparado para AnnotSV
    bcftools index "$PREPARED_VCF"

    # Ejecutar AnnotSV
    "$ANNOTSV_BIN" -SVinputFile "$PREPARED_VCF" -SVinputInfo 1 \
        -outputFile "$ANNOT_TSV" -annotationMode full

    # Comprimir y mover el archivo final
    gzip -c "$ANNOT_TSV" > "${OUTPUT_TSV_GZ}.tmp"
    mv -f "${OUTPUT_TSV_GZ}.tmp" "$OUTPUT_TSV_GZ"

    # Eliminar archivos temporales
    rm -f "$PREPARED_VCF" "${PREPARED_VCF}.csi" "$ANNOT_TSV"
}


# ---- NASERTIC ----
# 1) .tsv.gz ya existe -> omitir
# 2) .tsv sin comprimir -> comprimir
# 3) Solo VCF -> ejecutar AnnotSV y comprimir
annotate_or_compress_nasertic_file() {
    local FILE="$1"       # Archivo de entrada
    local FOLDER_PATH="$2"  # Carpeta de análisis
    local MODE="$3"       # CNV o SV
    local TMP_DIR="$4"    # Carpeta temporal para intermedios
    local ID FINAL_TSV_GZ

    # 1) Existe archivo TSV comprimido: omitir
    if [[ "$(basename "$FILE")" == *.tsv.gz ]]; then
        echo "[INFO] Archivo ya anotado y comprimido: $(basename "$FILE")"
        return 0

    # 2) Existe archivo TSV sin comprimir: comprimir
    elif [[ "$(basename "$FILE")" == *.tsv ]]; then
        # Si existe también el archivo comprimido, omitir para evitar fallo al comprimir
        if [[ -e "${FILE}.gz" ]]; then
            echo "[WARNING] Ya existe $(basename "${FILE}.gz"), se omite compresión de $(basename "$FILE")"
            return 0
        fi

        echo "[INFO] Comprimiendo $(basename "$FILE")"
        gzip "$FILE"
        return 0

    # 3) No existe, solo está el VCF: ejecutar AnnotSV y comprimir
    elif [[ "$(basename "$FILE")" == *.vcf || "$(basename "$FILE")" == *.vcf.gz ]]; then
        # Obtener ID de la muestra (ej. 3984-3984-4impact-04 a partir de 3984-3984-4impact-04.cnv.pon.vcf)
        ID=$(basename "$FILE" | sed -E 's/\.cnv\.pon\.vcf(\.gz)?$//; s/\.cnv\.to\.vcf(\.gz)?$//; s/\.sv\.vcf(\.gz)?$//; s/\.vcf(\.gz)?$//')

        # Ruta final esperada para el archivo comprimido según el tipo de variante
        if [[ "$MODE" == "CNV" ]]; then
            # TO: *.cnv.to.vcf / PON: *.cnv.pon.vcf
            if [[ "$FILE" == *.cnv.pon.vcf.gz ]]; then
                FINAL_TSV_GZ="${FOLDER_PATH}/${ID}.CNVs.PON.annotated.PASSfiltered.tsv.gz"
            else
                FINAL_TSV_GZ="${FOLDER_PATH}/${ID}.CNVs.TO.annotated.PASSfiltered.tsv.gz"
            fi
        else
            FINAL_TSV_GZ="${FOLDER_PATH}/${ID}.SVs.TO.annotated.PASSfiltered.tsv.gz"
        fi

        # Si ya existe el archivo final esperado, omitir este archivo original
        if [[ -s "$FINAL_TSV_GZ" ]]; then
            echo "[WARNING] Ya existe el archivo anotado y comprimido para $(basename "$FILE"), se omite"
            return 0
        fi

        # Ejecutar AnnotSV (filtrando por PASS, anotando y comprimiendo)
        echo "[INFO] Ejecutando AnnotSV para $(basename "$FILE" | sed -E 's/\.vcf(\.gz)?$//') ..."
        run_annotsv_tool "$FILE" "$FINAL_TSV_GZ" "$TMP_DIR"
        echo "[INFO] $ID completado"
        return 0
    fi

    echo "[WARNING] Archivo ignorado (extensión no reconocida): $FILE"
    return 0
}

run_nasertic() {
    local FOLDER MODE
    local FOUND_ANY=0

    # Procesar CNV TO y SV TO
    for FOLDER in "${WORKDIR%/}/CNV/INPUTS/SOMATIC" "${WORKDIR%/}/SV/INPUTS/SOMATIC"; do
        if [[ "$FOLDER" == *"/CNV/INPUTS/SOMATIC" ]]; then
            MODE="CNV"
        else
            MODE="SV"
        fi

        if [[ ! -d "$FOLDER" ]]; then
            echo "[WARNING] No existe la carpeta $FOLDER"
            continue
        fi

        FOUND_ANY=1

        shopt -s nullglob
        local FILES=("$FOLDER"/*)

        if [[ ${#FILES[@]} -eq 0 ]]; then
            echo "[INFO] No hay archivos en $FOLDER"
            continue
        fi

        local TMP_DIR
        TMP_DIR=$(mktemp -d "${FOLDER}/.annotsv_tmp.XXXXXX")

        for FILE in "${FILES[@]}"; do
            annotate_or_compress_nasertic_file "$FILE" "$FOLDER" "$MODE" "$TMP_DIR"
        done

        rm -rf "$TMP_DIR"
    done

    if [[ "$FOUND_ANY" -eq 0 ]]; then
        echo "[WARNING] No existen carpetas somaticas para NASERTIC"
        exit 1
    fi
}


# ---- FPGMX ----
# 1) .tsv.gz ya existe -> omitir
# 2) .tsv sin comprimir -> comprimir
# 3) Solo VCF -> ejecutar AnnotSV y comprimir
annotate_or_compress_fpgmx_file() {
    local FILE="$1"  # Archivo VCF de entrada
    local FOLDER_PATH="$2"  # Carpeta de análisis
    local TMP_DIR="$3"  # Carpeta temporal para archivos intermedios
    local ID OUTPUT_BASENAME FINAL_TSV_GZ  # Variables auxiliares para nombres de salida

    # 1) Existe archivo TSV comprimido: omitir
    if [[ "$(basename "$FILE")" == *.tsv.gz ]]; then
        echo "[INFO] Archivo ya anotado y comprimido: $(basename "$FILE")"
        return 0

    # 2) Existe archivo TSV sin comprimir: comprimir
    elif [[ "$(basename "$FILE")" == *.tsv ]]; then
        # Si existe también el archivo comprimido, omitir para evitar fallo al comprimir
        if [[ -e "${FILE}.gz" ]]; then
            echo "[WARNING] Ya existe $(basename "${FILE}.gz"), se omite compresión de $(basename "$FILE")"
            return 0
        fi

        echo "[INFO] Comprimiendo $(basename "$FILE")"
        gzip "$FILE"
        return 0

    # 3) No existe, solo está el VCF: ejecutar AnnotSV y comprimir
    elif [[ "$(basename "$FILE")" == *.vcf || "$(basename "$FILE")" == *.vcf.gz ]]; then
        # Obtener nombre del archivo sin .vcf o .vcf.gz
        ID=$(basename "$FILE" | sed -E 's/\.vcf(\.gz)?$//')
        # Nombre del archivo de salida TSV
        OUTPUT_BASENAME="${ID}.annotated.tsv"
        # Ruta final esperada  para el archivo comprimido
        FINAL_TSV_GZ="${FOLDER_PATH}/${OUTPUT_BASENAME}.gz"

        # Si ya existe el archivo final esperado, omitir este archivo original
        if [[ -s "$FINAL_TSV_GZ" ]]; then
            echo "[INFO] Ya existe el archivo anotado y comprimido para $(basename "$FILE"), se omite"
            return 0
        fi

        # Ejecutar AnnotSV (filtrando por PASS, anotando y comprimiendo)
        echo "[INFO] Ejecutando AnnotSV para $ID ..."
        run_annotsv_tool "$FILE" "$FINAL_TSV_GZ" "$TMP_DIR"
        echo "[INFO] $ID completado"

    else
        echo "[WARNING] Archivo ignorado (extensión no reconocida): $FILE"
        return 0

    fi
}

# Comprobar/procesar archivos PON y TN
check_fpgmx_somatic_files() {
    local FOLDER VCF_FILE TMP_DIR  # Variables locales para iterar carpetas, entradas y temporal

    # Recorrer carpetas de CNV y SV para seleccionar archivos PON y TN
    for FOLDER in "${WORKDIR%/}/INPUTS/SOMATIC_CNV" "${WORKDIR%/}/INPUTS/SOMATIC_SV"; do
        [[ ! -d "$FOLDER" ]] && continue

        shopt -s nullglob
        local FILES=()
        for FILE in "$FOLDER"/*; do
            local BN
            BN=$(basename "$FILE")
            if [[ "$BN" == *_D_* || "$BN" == *DTN_* ]]; then
                FILES+=("$FILE")
            fi
        done

        [[ ${#FILES[@]} -eq 0 ]] && continue

        # Crear carpeta temporal
        TMP_DIR=$(mktemp -d "${FOLDER}/.annotsv_ref_tmp.XXXXXX")
        # Procesa archivos TN y PON detectados
        for FILE in "${FILES[@]}"; do
            annotate_or_compress_fpgmx_file "$FILE" "$FOLDER" "$TMP_DIR"
        done
        # Borrar carpeta temporal
        rm -rf "$TMP_DIR"
    done
}

run_fpgmx() {
    # Ruta a la carpeta de entrada de CNVs y SVs
    local FOLDERS=("${WORKDIR%/}/INPUTS/SOMATIC_SV" "${WORKDIR%/}/INPUTS/SOMATIC_CNV")

    # Comprobar/procesar archivos PON y TN
    check_fpgmx_somatic_files

    # Variables de control para procesar archivos de TO
    local failures=0 processed_any_folder=0 pending_any_file=0

    # Recorre las carpetas de CNVs y SVs para procesar los archivos de TO
    for FOLDER_PATH in "${FOLDERS[@]}"; do
        if [[ ! -d "$FOLDER_PATH" ]]; then
            echo "[WARNING] No existe la carpeta $FOLDER_PATH"
            exit 1
        fi

        processed_any_folder=1

        # Seleccionar archivos VCF
        shopt -s nullglob
        local VCF_FILES=("$FOLDER_PATH"/*.vcf "$FOLDER_PATH"/*.vcf.gz)

        if [[ ${#VCF_FILES[@]} -eq 0 ]]; then
            echo "[ERROR] Archivos .vcf o .vcf.gz no encontrados en $FOLDER_PATH"
            continue
        fi

        local VCF_FILE ID OUTPUT_BASENAME EXPECTED_OUTPUT  # Variables locales para evaluar entradas/salidas
        local SKIPPED=0  # Contador de archivos anotados
        local PENDING=()  # Lista de archivos pendientes para AnnotSV

        # Seleccionar archivos para anotar
        for VCF_FILE in "${VCF_FILES[@]}"; do
            ID=$(basename "$VCF_FILE" | sed -E 's/\.sv\.vcf(\.gz)?$//; s/\.vcf(\.gz)?$//')

            # Procesar muestras de Tumor Only
            [[ "$ID" != *DTO_* ]] && continue

            # Definir nombre de salida segun si es CNV o SV
            if [[ "$VCF_FILE" == *.sv.vcf || "$VCF_FILE" == *.sv.vcf.gz ]]; then
                OUTPUT_BASENAME="${ID}.sv.filtered_annotated.tsv"
            else
                OUTPUT_BASENAME="${ID}.filtered_annotated.tsv"
            fi

            # Ruta al archivo comprimido
            EXPECTED_OUTPUT="${FOLDER_PATH}/${OUTPUT_BASENAME}.gz"

            # Omitir el archivo si ya existe, añadir en caso contrario
            if [[ -s "$EXPECTED_OUTPUT" ]]; then
                SKIPPED=$((SKIPPED + 1))
            else
                PENDING+=("$VCF_FILE")
            fi
        done

        [[ "$SKIPPED" -gt 0 ]] && echo "[INFO] $SKIPPED archivos ya anotados, omitidos"

        if [[ ${#PENDING[@]} -eq 0 ]]; then
            echo "[INFO] No hay archivos TO pendientes de anotar"
            continue
        fi

        pending_any_file=1

        # Carpeta temporal para archivos intermedios de AnnotSV
        local TMP_DIR
        TMP_DIR=$(mktemp -d "${FOLDER_PATH}/.annotsv_tmp.XXXXXX")

        echo "[INFO] Procesando ${#PENDING[@]} archivo(s) TO con AnnotSV..."

        # Ejecutar AnnotSV para cada archivo TO pendiente
        for VCF_FILE in "${PENDING[@]}"; do
            ID=$(basename "$VCF_FILE" | sed -E 's/\.sv\.vcf(\.gz)?$//; s/\.vcf(\.gz)?$//')
            # Definir nombre de salida según si es CNV o SV
            if [[ "$VCF_FILE" == *.sv.vcf || "$VCF_FILE" == *.sv.vcf.gz ]]; then
                OUTPUT_BASENAME="${ID}.sv.filtered_annotated.tsv"
            else
                OUTPUT_BASENAME="${ID}.filtered_annotated.tsv"
            fi
            # Nombre del archivo generado comprimido
            EXPECTED_OUTPUT="${FOLDER_PATH}/${OUTPUT_BASENAME}.gz"

            echo "[INFO] Ejecutando AnnotSV para $ID ..."
            run_annotsv_tool "$VCF_FILE" "$EXPECTED_OUTPUT" "$TMP_DIR"
            echo "[INFO] $ID completado"
        done

        # Borrar carpeta temporal
        rm -rf "$TMP_DIR"
    done
}

# ---- Ejecutar según nodo ----
case "$NODE" in
    NASERTIC) run_nasertic ;;
    FPGMX)    run_fpgmx ;;
    CNAG)     echo "[INFO] Archivos somáticos de CNAG: No se requiere ejecutar AnnotSV"; exit 0 ;;
esac
