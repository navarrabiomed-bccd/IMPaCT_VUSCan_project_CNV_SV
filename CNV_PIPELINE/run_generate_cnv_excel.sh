#!/bin/bash
# Regenerate the CNV final report (step 3) for families with
# intermediate files available.
# Runs in parallel with a configurable limit.


# Validate arguments
if [ "$#" -lt 2 ]; then
	echo "Uso: $0 <NODO> <FAMILIA> [FAMILIA2 ...]"
	echo "      $0 <NODO> --all"
	echo "Ejemplo (una familia):   $0 CNAG 5501"
	echo "Ejemplo (varias):        $0 CNAG 5501 5502 5503"
	echo "Ejemplo (todas):         $0 CNAG --all"
	exit 1
fi

NODO="$1"
BASEDIR="/mnt/c/Users/bioinfor/Desktop/VUSCAN_PIPELINES/DATA/SAMPLES/${NODO}/"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
shift

# Maximum number of families to process in parallel
MAX_PARALLEL=4

# Collect family list
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

# Keep only numeric folder names
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
	echo "Error: no hay familias válidas para procesar."
	exit 1
fi

echo "Familias a procesar (${#FAMILIAS[@]}): ${FAMILIAS[*]}"

FAILED=0

run_family() {
	local FAMILIA="$1"
	local WORKDIR="${BASEDIR}/${FAMILIA}/"
	local LOGFILE="${BASEDIR}/${FAMILIA}/${FAMILIA}_excel_adding_freq.log"

	if [ ! -d "$WORKDIR" ]; then
		echo "[${FAMILIA}] Error: directorio no existe: $WORKDIR"
		return 1
	fi

	{
		echo "[${FAMILIA}] Iniciando generación del Excel final..."
		conda run --no-capture-output -n pipeline \
			python3 "${SCRIPT_DIR}/generate_cnv_excel_only.py" "$WORKDIR"
		echo "[${FAMILIA}] Excel generado correctamente."
	} > >(tee "$LOGFILE") 2>&1
}

wait_oldest() {
	local PID="$1"
	local FAMILIA="$2"

	if wait "$PID"; then
		echo "Familia ${FAMILIA}: completada."
	else
		echo "Familia ${FAMILIA}: FALLÓ (revisa ${BASEDIR}/${FAMILIA}/${FAMILIA}_excel_adding_freq.log)."
		FAILED=1
	fi
}

# Run families in parallel
PIDS=()
PENDING_FAMILIAS=()
for FAMILIA in "${FAMILIAS[@]}"; do
	while (( ${#PIDS[@]} >= MAX_PARALLEL )); do
		wait_oldest "${PIDS[0]}" "${PENDING_FAMILIAS[0]}"
		PIDS=("${PIDS[@]:1}")
		PENDING_FAMILIAS=("${PENDING_FAMILIAS[@]:1}")
	done
	echo "Iniciando familia: $FAMILIA"
	run_family "$FAMILIA" &
	PIDS+=($!)
	PENDING_FAMILIAS+=("$FAMILIA")
done

# Wait for any remaining running families
while (( ${#PIDS[@]} > 0 )); do
	wait_oldest "${PIDS[0]}" "${PENDING_FAMILIAS[0]}"
	PIDS=("${PIDS[@]:1}")
	PENDING_FAMILIAS=("${PENDING_FAMILIAS[@]:1}")
done

if [ "$FAILED" -eq 0 ]; then
	echo "Todas las familias generadas correctamente."
else
	echo "Algunas familias fallaron. Revisa los logs en cada carpeta."
fi

exit $FAILED
