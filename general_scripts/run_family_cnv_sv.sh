#!/bin/bash

# Orchestrator script to run CNV and SV pipelines for VUSCan families.
# - Receives a node and one or more families (or --all).
# - Runs the CNV/SV pipelines and SNV input preparation per family.
# - Processes families in parallel.
# - Saves logs per family and moves the completed folders to "PENDIENTES_SUBIR".


# Argument validation: node and at least one family ID or --all
if [ "$#" -lt 2 ]; then
	echo "Uso: $0 <NODO> <FAMILIA> [FAMILIA2 ...]"
	echo "      $0 <NODO> --all"
	echo "Ejemplo (una familia):   $0 FPGMX 1042"
	echo "Ejemplo (varias):        $0 FPGMX 1042 1043 1044"
	echo "Ejemplo (todas):         $0 FPGMX --all"
	exit 1
fi

NODO="$1"
BASEDIR="/mnt/c/Users/bioinfor/Desktop/VUSCAN_PIPELINES/DATA/SAMPLES/VUSCan_families/${NODO}"
ANALIZADA_DIR="${BASEDIR}/PENDIENTES_SUBIR"
shift

# If --all is requested, find all families in the directory; otherwise, use provided family IDs
if [ "$1" = "--all" ]; then
	if [ ! -d "$BASEDIR" ]; then
		echo "Error: directorio del nodo no existe: $BASEDIR"
		exit 1
	fi
	mapfile -t FAMILIAS < <(find "$BASEDIR" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort)
	if [ ${#FAMILIAS[@]} -eq 0 ]; then
		echo "Error: no se encontraron familias en $BASEDIR"
		exit 1
	fi
	echo "Familias encontradas en $BASEDIR: ${FAMILIAS[*]}"
else
	FAMILIAS=("$@")
fi

# Keep only family folder names that are numeric
FAMILIAS_FILTERED=()
for F in "${FAMILIAS[@]}"; do
	if [[ "$F" =~ ^[0-9]+$ ]]; then
		FAMILIAS_FILTERED+=("$F")
	else
		echo "Carpeta ignorada (nombre no numérico): $F"
	fi
done
FAMILIAS=("${FAMILIAS_FILTERED[@]}")

if [ ${#FAMILIAS[@]} -eq 0 ]; then
	echo "Error: no hay familias en el directorio para procesar."
	exit 1
fi

# Ensure destination folder exists for processed families
mkdir -p "$ANALIZADA_DIR"

# Global error flag: if any family fails, final exit code will be 1
FAILED=0

# Maximum number of families to process in parallel
MAX_PARALLEL=2

# Function to run the full workflow for one family
run_family() {
	local FAMILIA="$1"
	local WORKDIR="${BASEDIR}/${FAMILIA}/"

	if [ ! -d "$WORKDIR" ]; then
		echo "[${FAMILIA}] Error: WORKDIR no existe: $WORKDIR"
		exit 1
	fi

	local LOGFILE="${WORKDIR}/${FAMILIA}_family.log"

	{
		echo "Iniciando análisis de la familia: ${FAMILIA}"
		set -e

		# STEP 0: Annotate somatic files when necessary
		echo "[${FAMILIA}] Ejecutando anotación de archivos de somático ..."
		conda run --no-capture-output -n annotsv_env bash run_annotsv.sh "$NODO" "$WORKDIR" 2>&1 | grep -E "ERROR|WARNING|No hay|Comprimiendo|Ejecutando|completado|CNAG" || true

		# ------ CNVs ------
		echo "[${FAMILIA}] CNVs: Iniciando ejecución ..."

		# STEP 1: CNV preprocessing
		echo "[${FAMILIA}] CNVs: Ejecutando preprocesado ..."
		conda run --no-capture-output -n pipeline python3 ../CNV_PIPELINE/cnv_preprocess.py "$WORKDIR"

		# STEP 2.0: ClassifyCNV execution
		echo "[${FAMILIA}] CNVs: Ejecutando ClassifyCNV ..."
		conda run --no-capture-output -n classifycnv_env bash ../CNV_PIPELINE/run_classifycnv.sh "$WORKDIR"

		# STEPS 2-3: CNV annotation and report generation
		echo "[${FAMILIA}] CNVs: Ejecutando priorización y generación del informe ..."
		conda run --no-capture-output -n pipeline python3 ../CNV_PIPELINE/cnv_annotation_reporting.py "$WORKDIR"

		echo "[${FAMILIA}] Ejecución finalizada de CNVs"

		# ------ SV ------
		echo "[${FAMILIA}] SVs: Iniciando ejecución ..."
		conda run --no-capture-output -n pipeline python3 ../SV_PIPELINE/sv_pipeline.py "$WORKDIR"

		echo "[${FAMILIA}] Ejecución finalizada de SVs"

		# ------ SNV ------
		# Prepare files required by the SNV pipeline
		echo "[${FAMILIA}] CNVs y SVs: Preparando archivos para el pipeline de SNVs ..."
		conda run --no-capture-output -n pipeline python3 preprocess_cnv_sv_vuscan.py "$WORKDIR"
		echo "[${FAMILIA}] CNVs y SVs: Archivos preparados para el pipeline de SNVs"

		# Duplicate output: terminal + log file
	} > >(tee "$LOGFILE") 2>&1
}

# Function to wait for a family job, report status, and move folder on success
finalize_family() {
	local PID="$1"
	local FAMILIA="$2"

	if wait "$PID"; then
		echo "Familia ${FAMILIA}: completada."

		SRC_DIR="${BASEDIR}/${FAMILIA}"
		DST_DIR="${ANALIZADA_DIR}/${FAMILIA}"

		if [ -d "$SRC_DIR" ]; then
			if [ -e "$DST_DIR" ]; then
				# Avoid overwriting an existing destination by appending a timestamp
				TS=$(date +%Y%m%d_%H%M%S)
				DST_DIR="${ANALIZADA_DIR}/${FAMILIA}_${TS}"
				echo "Familia ${FAMILIA}: destino ya existe, moviendo a ${DST_DIR}"
			fi

			if mv "$SRC_DIR" "$DST_DIR"; then
				echo "Familia ${FAMILIA}: carpeta movida a ${DST_DIR}"
			else
				echo "Familia ${FAMILIA}: Error al mover carpeta a 'PENDIENTES_SUBIR'"
				FAILED=1
			fi
		else
			echo "Familia ${FAMILIA}: No se encontró la carpeta para mover"
			FAILED=1
		fi
	else
		echo "Familia ${FAMILIA}: FALLÓ."
		FAILED=1
	fi
}


# Run families in parallel
PIDS=()
PENDING_FAMILIAS=()
for FAMILIA in "${FAMILIAS[@]}"; do
	# If the concurrency limit is reached, wait for the oldest running family
	while (( ${#PIDS[@]} >= MAX_PARALLEL )); do
		finalize_family "${PIDS[0]}" "${PENDING_FAMILIAS[0]}"
		PIDS=("${PIDS[@]:1}")
		PENDING_FAMILIAS=("${PENDING_FAMILIAS[@]:1}")
	done
	run_family "$FAMILIA" &
	PIDS+=($!)
	PENDING_FAMILIAS+=("$FAMILIA")
done

# Wait for remaining families and propagate global failure state
while (( ${#PIDS[@]} > 0 )); do
	finalize_family "${PIDS[0]}" "${PENDING_FAMILIAS[0]}"
	PIDS=("${PIDS[@]:1}")
	PENDING_FAMILIAS=("${PENDING_FAMILIAS[@]:1}")
done

exit $FAILED