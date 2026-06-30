# -*- coding: utf-8 -*-
"""
CNV pipeline - STEP 1: Process and parse CNV data from different sequencing
nodes (CNAG, FPGMX, NASERTIC).
Inputs: Files with germinal CNVs annotated with AnnotSV for a family.

Last modification: April 2026
"""

import os
import csv
import gzip
import pickle
import sys
import logging
from natsort import natsorted
import pandas as pd


# -----------------------------------------------
# PARAMETERS CONFIGURATION
# -----------------------------------------------

OVERLAP_THRESHOLD = 0.7   # 70% reciprocal overlap
QUAL_MIN = 30             # Minimum quality QUAL
INTERNAL_FREQ_MAX = 10    # Maximum internal frequency count
SAMPLE_PHENO_FILE = '/mnt/c/Users/bioinfor/Desktop/VUSCAN_PIPELINES/DATA/muestras_fenotipo.txt'

# Patrones de nombre para archivos CNAG
CNAG_FILE_PATTERNS = {
    '20k': ('fixed.tab.gz', '20000.CNVs.p.value.annotated.IntFreq.tab.gz'),
    '5k': ('5000.CNVs.p.value.annotated.tsv','5000.CNVs.p.value.annotated.IntFreq.tab.gz')
}


# -----------------------------------------------
# LOGGING SETUP
# -----------------------------------------------

logging.basicConfig(
    level=logging.DEBUG,  # Minimum level to display
    format='%(asctime)s - [ %(levelname)s ] - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)


# -----------------------------------------------
# Funciones de lectura y guardado
# -----------------------------------------------

def set_csv_field_size():
    """Configura el límite de tamaño de campo CSV de forma segura."""
    max_size = sys.maxsize
    while True:
        try:
            csv.field_size_limit(max_size)
            break
        except OverflowError:
            max_size = int(max_size / 10)


def load_phenotype_data(file_path, family_id):
    """Carga un archivo con información de fenotipos, filtrando por familia y
    devuelve un diccionario indexado por ID de muestra."""
    try:
        df = pd.read_csv(file_path, sep='\t', encoding='utf-8')
    except UnicodeDecodeError:
        df = pd.read_csv(file_path, sep='\t', encoding='ISO-8859-1')
    
    pattern = f'-{family_id}-'
    df = df[df['Muestra ID'].astype(str).str.contains(pattern)]

    if df.empty:
        raise ValueError('No se encontraron muestras para esta familia')
    
    # Convertir a diccionario (clave: ID IMPaCT)
    pheno_dict = df.set_index('Muestra ID').T.to_dict()
    return pheno_dict


def build_sample_id_maps(pheno_dict):
    """
    Construye los mapas de conversión entre IDs de muestras de CNAG e IMPaCT.

    Args:
        pheno_dict (dict): Diccionario de fenotipos de las muestras.

    Returns:
        tuple: imp_to_cnag (IMPaCT a CNAG), cnag_to_imp (CNAG a IMPaCT)
    """
    imp_to_cnag, cnag_to_imp = {}, {}

    for imp_id, data in pheno_dict.items():
        cnag_id = data.get('ID_CNAG')
        # Construir diccionarios con el mapa de IDs
        if cnag_id:
            imp_to_cnag[imp_id] = cnag_id
            cnag_to_imp[cnag_id] = imp_id

    return imp_to_cnag, cnag_to_imp


def save_pickle(data, file_path):
    """Guarda un archivo en formato pickle en la ruta indicada."""
    with open(file_path, 'wb') as f:
        pickle.dump(data, f)


def assert_full_split_id_match(full_dict, split_dict, sample_id):
    """Comprueba si las CNVs 'full' y 'split' son las mismas. Detiene la
    ejecución si no coinciden, excepto si son CNVs sin genes asociados."""
    full_ids = set(full_dict.keys())
    split_ids = set(split_dict.keys())

    only_full = full_ids - split_ids
    only_split = split_ids - full_ids

    # Permitir CNVs solo en 'full' cuando no tienen genes asociados
    allowed_only_full = set()
    for cnv_id in only_full:
        gene_count = full_dict.get(cnv_id, {}).get('Gene_count', None)
        try:
            if float(gene_count) == 0:
                allowed_only_full.add(cnv_id)
        except (TypeError, ValueError):
            pass

    unexpected_only_full = only_full - allowed_only_full

    if unexpected_only_full or only_split:
        only_full_preview = ', '.join(sorted(unexpected_only_full)) if unexpected_only_full else '-'
        only_split_preview = ', '.join(sorted(only_split)) if only_split else '-'
        logging.error(
            f'Muestra {sample_id}: CNVs full/split inconsistentes. '
            f'Solo en "full" ({len(unexpected_only_full)}): {only_full_preview} | '
            f'Solo en "split" ({len(only_split)}): {only_split_preview}'
        )
        raise SystemExit(1)
    
    else:
        logging.info(f'{sample_id}: {len(full_dict)} CNVs "full", {len(split_dict)} CNVs "split". '
                     f'CNVs "full" y "split" consistentes. ')


def save_outputs(bed_dict, full_dict, split_dict, header, folder_outputs, sample_id, add_id_column=False):
    """
    Guarda CNVs 'full' y 'split' en formato TSV, BED y pickle.
    
    Args:
        bed_dict (dict): Coordenadas CNVs {cnv_id: {chrom, start, end, type}}
        full_dict (dict): CNVs 'full' {cnv_id: {valor, ...}}
        split_dict (dict): CNVs 'split' {cnv_id: {valor, ...}}
        header (list): Nombres de columnas para el archivo TSV.
        folder_outputs (str): Ruta a la carpeta de salida.
        sample_id (str): ID completo de la muestra.
        add_id_column (bool): Si True, añade CNV ID en la primera columna (default: False).
    """
    # Verificar si las CNV 'full' y 'split' coinciden
    assert_full_split_id_match(full_dict, split_dict, sample_id)

    # Nombre de los archivos para guardar las CNVs 'full'
    output_tsv = os.path.join(folder_outputs, f"{sample_id}.CNVs.annotated_parsed.tsv")
    output_bed = os.path.join(folder_outputs, f"{sample_id}.CNVs.annotated_parsed.bed")

    # Ordenar CNVs por cromosoma y posición de inicio
    sorted_ids = natsorted(
        bed_dict.keys(), 
        key=lambda cnv_id: (bed_dict[cnv_id]['chrom'], int(bed_dict[cnv_id]['start']))
    )

    # Guardar coordenadas de las CNVs en un archivo BED
    with open(output_bed, 'w', newline='') as bedfile:
        bedfile.write('#chrom\tstart\tend\ttype\n')
        for cnv_id in sorted_ids:
            data = bed_dict[cnv_id]
            bedfile.write(f"{data['chrom']}\t{data['start']}\t{data['end']}\t{data['type']}\n")

    # Encabezado archivo TSV: añadir 'AnnotSV_ID' al inicio si no está incluido
    header_tsv = ['AnnotSV_ID'] + [f for f in header if f != 'AnnotSV_ID'] if add_id_column else header

    # Guardar CNVs en un archivo TSV
    with open(output_tsv, 'w', newline='') as tsvfile:
        writer = csv.DictWriter(tsvfile, fieldnames=header_tsv, delimiter='\t')
        writer.writeheader()
        for cnv_id in sorted_ids:
            row = {'AnnotSV_ID': cnv_id, **full_dict[cnv_id]}
            writer.writerow(row)

    # Guardar diccionarios 'full' y 'split' en archivos pickle
    pickle_file = os.path.join(folder_outputs, f"{sample_id}.CNVs.annotated_parsed.pkl")
    save_pickle(full_dict, pickle_file)

    pickle_split_file = pickle_file.replace('.pkl', '_split.pkl')
    save_pickle(split_dict, pickle_split_file)


def build_coord_from_row(row):
    """Construye coordenadas BED con datos de una CNV."""
    return {
        'chrom': row['SV_chrom'],
        'start': row['SV_start'],
        'end': row['SV_end'],
        'type': row['SV_type']
    }


# -------------------------------------------------------
# Funciones de solapamiento de variantes
# -------------------------------------------------------

# def ensure_sv_length(entry, start_key='SV_start', end_key='SV_end', length_key='SV_length'):
#     """
#     Comprueba o calcula que una variante tenga valor de longitud.

#     Args:
#         entry (dict): Fila/diccionario con los datos de la variante.
#         start_key (str): Clave donde buscar el inicio (default: 'SV_start').
#         end_key (str): Clave donde buscar el fin (default: 'SV_end').
#         length_key (str): Clave con la longitud (default: 'SV_length').

#     Returns:
#         dict: La fila/diccionario con la longitud correcta.
#     """
#     try:
#         length = int(float(entry.get(length_key, '')))
#     except (TypeError, ValueError):
#         start = int(entry[start_key])
#         end = int(entry[end_key])
#         length = max(0, end - start)
#         entry[length_key] = str(length)

#     return entry


def get_start_end(data):
    """Extrae las posiciones de inicio y fin de una variante."""
    try:
        start = int(data['SV_start'])
        end = int(data['SV_end'])
    except KeyError:
        start = int(data['start'])
        end = int(data['end'])
    return start, end


def check_reciprocal_overlap(data1, data2, threshold=OVERLAP_THRESHOLD):
    """
    Verifica si hay superposición recíproca ≥ umbral entre dos variantes.
    
    Args:
        data1 (dict): Datos variante 1.
        data2 (dict): Datos variante 2.
        threshold (float): Umbral mínimo de superposición recíproca.
    
    Returns:
        bool: True si hay superposición recíproca ≥ umbral, False si no.
    """
    # Obtener información de las variantes
    start1, end1 = get_start_end(data1)
    start2, end2 = get_start_end(data2)
    
    # Calcular la longitud de la superposición
    overlap = max(0, min(end1, end2) - max(start1, start2) + 1)

    # Obtener la longitud de las variantes
    len1 = max(abs(float(data1.get('SV_length', end1 - start1 + 1))), 1)
    len2 = max(abs(float(data2.get('SV_length', end2 - start2 + 1))), 1)

    # Verificar si la superposición recíproca es al menos el umbral
    return (overlap / len1 >= threshold and overlap / len2 >= threshold)


def cnv_overlaps_with_coord(row, coord_dict):
    """Comprueba solapamiento de una CNV con un diccionario de coordenadas."""
    return any(
        row['SV_chrom'] == coord['chrom']
        and row['SV_type'] == coord['type']
        and check_reciprocal_overlap(row, coord)
        for coord in coord_dict.values()
    )


# -------------------------------------------------
# Funciones de procesamiento de archivos de CNVs
# -------------------------------------------------

# ----------- Archivos de CNAG -----------
def read_cnag_cnv_file(cnv_file, window_label):
    """
    Lee un archivo de CNVs de CNAG y filtra por frecuencia interna sobre las
    'full'. Las 'split' se incluyen solo si su 'full' asociada pasó el filtro.
    
    Args:
        cnv_file (str): Ruta al archivo TSV de CNVs.
        window_label (str): Etiqueta de ventana CNAG (ej. '20k', '5k').
    
    Returns:
        tuple: (coordenadas, CNVs full, CNVs split, header).
    """
    # Inicializar diccionarios
    coord_dict = {}  # Coordenadas CNVs (archivo BED)
    full_dict = {}   # CNVs 'full'
    split_dict = {}  # CNVs 'split'

    opener = gzip.open if cnv_file.endswith('.gz') else open
    with opener(cnv_file, 'rt', newline='') as file:
        reader = csv.DictReader(file, delimiter='\t')
        header = reader.fieldnames  # Guardar encabezado
        header = header + ['CNV_window']
        
        # Para cada CNV del archivo
        for row in reader:
            cnv_id = f"{row['SV_chrom']}_{row['SV_start']}_{row['SV_end']}_{row['SV_type']}"

            # CNVs 'full'
            if row['Annotation_mode'] == 'full':
                # # Filtrar 'PASS'
                # if row['FILTER'] != 'PASS':
                #     continue

                # Filtrar por frecuencia interna
                try:
                    row['Illumina.exact.counts'] = 0 if row['Illumina.exact.counts'] == '#' else float(row['Illumina.exact.counts'])
                    if float(row['Illumina.exact.counts']) >= INTERNAL_FREQ_MAX:
                        continue
                except ValueError:
                    logging.error(f'CNV {cnv_id}: su frecuencia interna no es un valor numérico ({row["Illumina.exact.counts"]})')

                # Convertir '#' a 0 en la frecuencia interna similar
                if row['Illumina.similar.counts'] == '#':
                    row['Illumina.similar.counts'] = 0

                # Guardar las coordenadas en el diccionario BED
                coord_dict[cnv_id] = build_coord_from_row(row)

                # Anotar ventana origen (solo CNAG)
                row['CNV_window'] = window_label

                # Guardar con toda la información
                full_dict[cnv_id] = {k: v for k, v in row.items() if k != 'AnnotSV_ID'}

            # CNVs 'split'
            elif row['Annotation_mode'] == 'split':
                # Guardar CNV si su 'full' asociada pasó los filtros
                if cnv_id in full_dict:
                    # si la frecuencia interna es '#', reemplazar con la de su 'full'
                    if row.get('Illumina.exact.counts') == '#':
                        full_freq = full_dict[cnv_id].get('Illumina.exact.counts')
                        row['Illumina.exact.counts'] = full_freq
                    if row.get('Illumina.similar.counts') == '#':
                        full_freq = full_dict[cnv_id].get('Illumina.similar.counts')
                        row['Illumina.similar.counts'] = full_freq
                    split_dict.setdefault(cnv_id, []).append(row)

    return coord_dict, full_dict, split_dict, header


def merge_cnag_windows(base_coord, base_full, base_split, extra_coord, extra_full, extra_split):
    """Combina diccionarios de CNVs de dos ventanas de CNAG añadiendo
    las CNVs de la segunda no encontradas en la primera."""
    merged_coord = dict(base_coord)
    merged_full = dict(base_full)
    merged_split = dict(base_split)

    for cnv_id, cnv_data in extra_full.items():
        already_covered = cnv_overlaps_with_coord(cnv_data, merged_coord)

        if not already_covered:
            merged_coord[cnv_id] = extra_coord[cnv_id]
            merged_full[cnv_id] = cnv_data
            if cnv_id in extra_split:
                merged_split[cnv_id] = extra_split[cnv_id]

    return merged_coord, merged_full, merged_split


def process_cnag_proband(file_20k, file_5k, proband, folder_outputs):
    """
    Procesa CNVs del probando de CNAG. Combina dos ventanas (20 y 5 kbp)
    eliminando duplicados por solapamiento recíproco y aplica filtros
    ('PASS' y frecuencia interna). Guarda las CNVs en TSV, BED y pickle.

    - CNAG no proporciona QUAL - no se aplica este filtro de calidad.
    - Valores '#' en frecuencia interna se pasan a 0.
    - CNVs 'split' se filtran según las 'full', ya que no tienen frecuencia.

    Args:
        file_20k: Ruta al archivo de CNVs ventana 20,000 bp.
        file_5k: Ruta al archivo de CNVs ventana 5,000 bp.
        proband (str): ID del probando.
        folder_outputs (str): Ruta a la carpeta para guardar los resultados.

    Return:
        dict: Coordenadas CNVs {cnv_id: {chrom, start, end, type}}
    """
    # Procesar CNVs: ventana 20,000 bp
    coord_20k, full_20k, split_20k, header = read_cnag_cnv_file(file_20k, window_label='20 kbp')

    # Procesar CNVs: ventana 5,000 bp
    coord_5k, full_5k, split_5k, _ = read_cnag_cnv_file(file_5k, window_label='5 kbp')

    # Combinar ventanas eliminando duplicados por solapamiento recíproco
    prob_coord, prob_full, prob_split = merge_cnag_windows(
        coord_20k, full_20k, split_20k,
        coord_5k, full_5k, split_5k
    )

    # Guardar las CNVs (TSV y BED: ordenadas por cromosoma y posición de inicio)
    save_outputs(prob_coord, prob_full, prob_split, header, folder_outputs, proband)
    
    return prob_coord


def process_cnag_relatives(rels_files, prob_coord, cnag_to_imp, folder_outputs, output_file):
    """
    Procesa archivos TSV de CNVs de familiares de CNAG. Aplica filtros
    ('PASS' y frecuencia interna), combina ventanas y excluye solapamientos
    con el probando. Guarda las CNVs en TSV, BED y pickle.

    Args:
        rels_files (list): Lista de tuplas (cnag_id, file_20k, file_5k) de familiares.
        prob_coord (dict): Coordenadas de las CNVs del probando.
        cnag_to_imp (dict): Correspondencia de IDs de CNAG a IMPaCT.
        folder_outputs (str): Ruta a la carpeta para guardar los resultados.
        output_file (str): Archivo para guardar CNVs compartidas con el probando.
    """
    # Inicializar diccionario para guardar CNVs compartidas con el probando
    prob_common = {}

    for rel_cng, rel_file_20k, rel_file_5k in rels_files:
        rel_id = cnag_to_imp[rel_cng]

        # Procesar CNVs de ambas ventanas
        coord_20k, full_20k, split_20k, header = read_cnag_cnv_file(rel_file_20k, window_label='20 kbp')
        coord_5k, full_5k, split_5k, _ = read_cnag_cnv_file(rel_file_5k, window_label='5 kbp')

        # Combinar ventanas eliminando duplicados por solapamiento recíproco
        merged_coord, merged_full, merged_split = merge_cnag_windows(
            coord_20k, full_20k, split_20k,
            coord_5k, full_5k, split_5k
        )
        
        # Inicializar diccionarios
        rel_coord_dict = {}       # Coordenadas CNVs (archivo BED)
        rel_full_dict = {}        # CNVs 'full'
        rel_split_dict = {}       # CNVs 'split'
        prob_common[rel_id] = {}  # CNVs compartidas con el probando
        
        # Comprobar solapamientos con el probando
        for cnv_id, cnv_row in merged_full.items():
            if cnv_overlaps_with_coord(cnv_row, prob_coord):
                prob_common[rel_id][cnv_id] = cnv_row
                continue
            else:
                rel_coord_dict[cnv_id] = merged_coord[cnv_id]
                rel_full_dict[cnv_id] = cnv_row
                if cnv_id in merged_split:
                    rel_split_dict[cnv_id] = merged_split[cnv_id]
        
        # Guardar CNVs del familiar no presentes en el probando
        save_outputs(rel_coord_dict, rel_full_dict, rel_split_dict, header, folder_outputs, rel_id)
        logging.info(f'{rel_id}: {len(prob_common[rel_id])} CNVs compartidas con el probando.')

    # Guardar CNVs de los familiares compartidas con el probando
    if prob_common:
        save_pickle(prob_common, os.path.join(folder_outputs, output_file))


# ----------- Archivos de FPGMX -----------
def process_fpgmx_proband(full_file, split_file, proband, folder_outputs):
    """
    Procesa archivos de variantes (CNVs y SVs) del probando de FPGMX, aplicando
    varios filtros (CNVs, calidad y frecuencia interna). Guarda las CNVs 'full'
    y 'split' en archivos TSV, BED y pickle.

    Args:
        full_file (str): Ruta al archivo de CNVs 'full'.
        split_file (str): Ruta al archivo de CNVs 'split'.
        proband (str): ID probando.
        folder_outputs (str): Ruta a la carpeta de salida.

    Return:
        dict: Coordenadas CNVs {CNV ID: {chrom, start, end, type}}
    """
    # Inicializar diccionarios
    coord_dict = {}  # Coordenadas CNVs (archivo BED)
    full_dict = {}   # CNVs 'full'
    split_dict = {}  # CNVs 'split'
    
    #  Procesar CNVs 'full' y 'split'
    for cnv_file, is_full in [(full_file, True), (split_file, False)]:
        with gzip.open(cnv_file, 'rt') as file:
            reader = csv.DictReader(file, delimiter='\t')
            header = reader.fieldnames  # Guardar encabezado
            
            # Para cada CNV del archivo
            for row in reader:
                # Seleccionar CNVs
                if row['SV_type_original'] not in {'GAIN', 'LOSS'}:
                    continue

                # Filtrar por calidad
                if float(row['Qual']) < QUAL_MIN:
                    continue
                
                # Filtrar por frecuencia interna
                if int(row['Illumina_DRAGEN.exact.counts']) >= INTERNAL_FREQ_MAX:
                    continue
                
                # Obtener ID de la CNV
                cnv_id = f"{row['SV_chrom']}_{row['SV_start']}_{row['SV_end']}_{row['SV_type']}"

                # CNV 'full'
                if is_full and cnv_id not in full_dict:
                    # Guardar las coordenadas en el diccionario BED
                    coord_dict[cnv_id] = build_coord_from_row(row)

                    # Guardar toda la información
                    full_dict[cnv_id] = {k: v for k, v in row.items() if k != 'AnnotSV_ID'}

                # CNV 'split'
                else:
                    split_dict.setdefault(cnv_id, []).append(row)
    
    # Guardar CNVs 'full' y 'split'
    save_outputs(coord_dict, full_dict, split_dict, header, folder_outputs, proband, add_id_column=True)
    
    return coord_dict


def process_fpgmx_relatives(rels_full_files, prob_coord, folder_outputs, output_file):
    """
    Procesa archivos de CNVs y SVs de familiares de FPGMX. Aplica filtros
    (CNVs, calidad, y frecuencia interna) y elimina solapamientos con el
    probando. Guarda las CNVs 'full' y 'split' en archivos TSV, BED y pickle.

    Args:
        rels_full_files (list): Rutas a las CNVs 'full' de los familiares.
        prob_coord (dict): Coordenadas de las CNVs del probando.
        folder_outputs (str): Ruta para guardar los resultados.
        output_file (str): Archivo para guardar CNVs compartidas con el probando.
    """
    # Inicializar diccionario para CNVs compartidas con el probando
    prob_common = {}

    for rel_full_file in rels_full_files:
        # ID completo del familiar (ej. '2271-2270-4impact-01')
        rel_id = os.path.basename(rel_full_file).split('.')[0]
        
        rel_split_file = rel_full_file.replace('.full.tsv.gz', '.split.tsv.gz')
        
        # Inicializar diccionarios
        rel_coord_dict = {}       # Coordenadas CNVs (archivo BED)
        rel_full_dict = {}        # CNVs 'full'
        rel_split_dict = {}       # CNVs 'split'
        prob_common[rel_id] = {}  # CNVs compartidas con el probando
        
        # Procesar archivo del familiar
        for cnv_file, is_full in [(rel_full_file, True), (rel_split_file, False)]:
            with gzip.open(cnv_file, 'rt', encoding='iso-8859-1') as file:
                reader = csv.DictReader(file, delimiter='\t')
                header = reader.fieldnames

                for row in reader:
                    # Seleccionar CNVs
                    if row['SV_type_original'] not in {'GAIN', 'LOSS'}:
                        continue
                    # Filtrar por calidad
                    if float(row['Qual']) < QUAL_MIN:
                        continue
                    # Filtrar por frecuencia interna
                    if int(row['Illumina_DRAGEN.exact.counts']) >= INTERNAL_FREQ_MAX:
                        continue

                    cnv_id = f"{row['SV_chrom']}_{row['SV_start']}_{row['SV_end']}_{row['SV_type']}"
                    
                    # Comprobar solapamiento con el probando
                    overlaps_proband = cnv_overlaps_with_coord(row, prob_coord)

                    if overlaps_proband is False:
                        if is_full:
                            rel_coord_dict[cnv_id] = build_coord_from_row(row)
                            rel_full_dict[cnv_id] = {k: v for k, v in row.items() if k != 'AnnotSV_ID'}
                        else:
                            rel_split_dict.setdefault(cnv_id, []).append(row)
                    
                    elif overlaps_proband is True and is_full:
                        prob_common[rel_id][cnv_id] = {k: v for k, v in row.items() if k != 'AnnotSV_ID'}

        # Guardar CNVs 'full' y 'split' del familiar no presentes en el probando
        save_outputs(rel_coord_dict, rel_full_dict, rel_split_dict, header, folder_outputs,
                     rel_id, add_id_column=True)
        logging.info(f'{rel_id}: {len(prob_common[rel_id])} CNVs compartidas con el probando.')
    
    # Guardar CNVs de familiares compartidas con el probando
    if prob_common:
        save_pickle(prob_common, os.path.join(folder_outputs, output_file))


# ----------- Archivos de NASERTIC -----------
def process_nasertic_proband(input_file, proband, folder_outputs):
    """
    Procesa archivos de CNVs del probando de NASERTIC, aplicando varios filtros
    (calidad (QUAL y 'PASS') y frecuencia interna). Guarda lasCNVs 'full' y
    'split' en archivos TSV, BED y pickle.

    Args:
        input_file (str): Ruta al archivo de CNVs.
        proband (str): ID probando.
        folder_outputs (str): Ruta a la carpeta de salida.

    Return:
        dict: Coordenadas de CNVs {cnv_id: {chrom, start, end, type}}
    """
    # Inicializar diccionarios
    coord_dict = {}  # Coordenadas (archivo BED)
    full_dict = {}   # CNVs 'full'
    split_dict = {}  # CNVs 'split'

    # Leer el archivo
    with open(input_file, newline='') as file:
        reader = csv.DictReader(file, delimiter='\t')
        header = reader.fieldnames
 
        for row in reader:
            # Filtrar por calidad
            if float(row['QUAL']) < QUAL_MIN:
                continue
            # Filtrar 'PASS'
            if row['FILTER'] != 'PASS':
                continue
            # Filtrar por frecuencia interna
            try:
                if int(row['Illumina_DRAGEN.exact.counts']) >= INTERNAL_FREQ_MAX:
                    continue
            except ValueError:
                pass  # conservar valores no numéricos ('#')

            # Obtener ID de la CNV
            cnv_id = f"{row['SV_chrom']}_{row['SV_start']}_{row['SV_end']}_{row['SV_type']}"
            
            # Añadir CNV, si no está
            if row['Annotation_mode'] == 'full' and cnv_id not in full_dict:
                # Guardar coordenadas en el diccionario BED
                coord_dict[cnv_id] = build_coord_from_row(row)
                
                # Guardar toda la información
                full_dict[cnv_id] = {k: v for k, v in row.items() if k != 'AnnotSV_ID'}

            # CNV 'split'
            if row['Annotation_mode'] == 'split':
                split_dict.setdefault(cnv_id, []).append(row)

    # Guardar CNVs 'full' y 'split'
    save_outputs(coord_dict, full_dict, split_dict, header, folder_outputs, proband)
    
    return coord_dict


def process_nasertic_relatives(rels_files, prob_coord, folder_outputs, output_file):
    """
    Procesa archivos de CNVs de familiares de NASERTIC. Aplica filtros y
    elimina CNVs solapantes con el probando. Guarda las CNVs 'full' y 'split'
    en archivos TSV, BED y pickle.

    Args:
        rels_files (list): Rutas a las CNVs de los familiares.
        prob_coord (dict): Coordenadas de las CNVs del probando.
        folder_outputs (str): Ruta donde guardar los resultados.
        output_file (str): Archivo para guardar CNVs compartidas con el probando.
    """
    # Inicializar diccionario para guardar CNVs compartidas con el probando
    prob_common = {}

    for rel_file in rels_files:
        # ID completo del familiar
        rel_id = os.path.basename(rel_file).split('.')[0]
        
        # Inicializar diccionarios
        rel_coord_dict = {}       # Coordenadas CNVs (archivo BED)
        rel_full_dict = {}        # CNVs 'full'
        rel_split_dict = {}       # CNVs 'split'
        prob_common[rel_id] = {}  # CNVs compartidas con el probando
        
        # Procesar archivo del familiar
        with open(rel_file, newline='') as file:
            reader = csv.DictReader(file, delimiter='\t')
            header = reader.fieldnames

            for row in reader:
                # Aplicar filtros
                if float(row['QUAL']) < QUAL_MIN:
                    continue
                if row['FILTER'] != 'PASS':
                    continue
                if int(row['Illumina_DRAGEN.exact.counts']) >= INTERNAL_FREQ_MAX:
                    continue

                cnv_id = f"{row['SV_chrom']}_{row['SV_start']}_{row['SV_end']}_{row['SV_type']}"

                # Comprobar solapamiento con el probando
                overlaps_proband = cnv_overlaps_with_coord(row, prob_coord)

                if overlaps_proband is False:
                    # CNV 'full'
                    if row['Annotation_mode'] == 'full':
                        rel_coord_dict[cnv_id] = build_coord_from_row(row)
                        rel_full_dict[cnv_id] = {k: v for k, v in row.items() if k != 'AnnotSV_ID'}
                    # CNV 'split'
                    elif row['Annotation_mode'] == 'split':
                        rel_split_dict.setdefault(cnv_id, []).append(row)
                
                elif overlaps_proband is True and row['Annotation_mode'] == 'full':
                    prob_common[rel_id][cnv_id] = {k: v for k, v in row.items() if k != 'AnnotSV_ID'}

        # Guardar CNVs del familiar no presentes en el probando
        save_outputs(rel_coord_dict, rel_full_dict, rel_split_dict, header, folder_outputs, rel_id)
        logging.info(f'{rel_id}: {len(prob_common[rel_id])} CNVs compartidas con el probando')

    # Guardar CNVs de familiares compartidas con el probando
    if prob_common:
        save_pickle(prob_common, os.path.join(folder_outputs, output_file))


# -------------------------------
# Función principal
# -------------------------------

def validate_germinal_files(pheno_dict, input_files, node, imp_to_cnag=None):
    """
    Valida que existan los archivos de las muestras germinales esperadas.
    
    Args:
        pheno_dict (dict): Diccionario de fenotipos.
        input_files (list): Lista de archivos encontrados en entrada.
        node (str): Nodo de secuenciación ('CNAG', 'FPGMX', 'NASERTIC').
        imp_to_cnag (dict): Mapa IMPaCT→CNAG (solo para CNAG).
    
    Raises:
        SystemExit: Si faltan archivos de muestras germinales esperadas.
    """
    # Obtener las muestras germinales esperadas
    expected_imp_ids = {
        key for key, val in pheno_dict.items()
        if 'ADN germinal' in val.get('Muestra', '')
    }
    
    # Obtener los IDs encontrados
    if node == 'CNAG':
        found_ids = {f.split('_')[0].split('.')[0] for f in input_files
                     if f.endswith(CNAG_FILE_PATTERNS['20k']) or f.endswith(CNAG_FILE_PATTERNS['5k'])}
        expected_ids = {imp_to_cnag[imp_id] for imp_id in expected_imp_ids}
    elif node == 'FPGMX':
        found_ids = {f.split('.')[0] for f in input_files if f.endswith('full.tab.gz')}
        expected_ids = expected_imp_ids
    elif node == 'NASERTIC':
        found_ids = {f.split('.')[0] for f in input_files}
        expected_ids = expected_imp_ids
    
    # Detectar si faltan archivos
    missing_ids = sorted(expected_ids - found_ids)
    if missing_ids:
        raise SystemExit(
            f'[ ERROR ] - Faltan archivos de CNVs germinales de las muestras: '
            f'{", ".join(missing_ids)}'
        )


def main(folder_path):
    """
    Función principal para procesar CNVs de probando y familiares. Según el
    nodo de secuenciación, aplica filtros, evita duplicados por solapamiento
    y guarda los archivos necesarios para los análisis posteriores.

    Args:
        folder_path (str): Ruta a la carpeta de análisis de la familia.
    """
    # ID de la familia y nodo de secuenciación
    family_id = folder_path.split('/')[-2]
    node = folder_path.split('/')[-3]

    logging.info(f"Analizando familia: {family_id} | Nodo de secuenciación: {node}")
    logging.info("INICIANDO EL PREPROCESAMIENTO DE CNVs (PASO 1)")

    # Rutas a las carpetas de entrada y salida
    folder_inputs = (
        os.path.join(folder_path, 'INPUTS/') if node == 'FPGMX'
        else os.path.join(folder_path, 'CNV/INPUTS/')
    )
    folder_outputs = os.path.join(folder_path, 'CNV/INTERMEDIATES/')

    # Crear la carpeta para guardar los resultados si no existe
    if not os.path.exists(folder_outputs):
        os.makedirs(folder_outputs, exist_ok=True)
        logging.info(f'Creada carpeta "{folder_outputs}"')
    else:
        logging.warning(f'La carpeta "{folder_outputs}" ya existe')

    # Archivo correspondencia muestras-fenotipo
    pheno_dict = load_phenotype_data(SAMPLE_PHENO_FILE, family_id)

    # Nombre del archivo para guardar CNVs compartidas entre familiares y probando
    output_file = f'{family_id}.CNVs.relatives_overlap_proband.pkl'
    
    # Obtener ID del probando (ej. '1943-1943-4impact-01')
    matches = [
        key for key, val in pheno_dict.items()
        if key.startswith(f'{family_id}-{family_id}-4impact') and 'ADN germinal' in val.get('Muestra', '')
    ]
    proband = matches[0]
    
    # Archivos con las CNVs de probando y familiares
    input_files = [f for f in os.listdir(folder_inputs) if f.endswith(('.tsv', '.tab', '.tab.gz'))]
    if not input_files:
        logging.warning('No hay archivos esperados en el directorio')
        return

    # Procesar CNVs según nodo
    if node == 'CNAG':
        # Correspondencia entre IDs de CNAG
        imp_to_cnag, cnag_to_imp = build_sample_id_maps(pheno_dict)
        prob_cnag = imp_to_cnag[proband]
        
        # Validar que hay archivos para las muestras germinales esperadas
        validate_germinal_files(pheno_dict, input_files, node, imp_to_cnag)

        # Ruta a todos los archivos por muestra y ventana
        cnag_files = {}
        for f in input_files:
            if f.endswith(CNAG_FILE_PATTERNS['20k']):
                window = '20k'
            elif f.endswith(CNAG_FILE_PATTERNS['5k']):
                window = '5k'
            else:
                continue
            cnag_id = f.split('_')[0].split('.')[0]
            cnag_files.setdefault(cnag_id, {})[window] = os.path.join(folder_inputs, f)

        # Procesar CNVs probando
        prob_file_20k = cnag_files.get(prob_cnag, {}).get('20k')
        prob_file_5k = cnag_files.get(prob_cnag, {}).get('5k')
        if not prob_file_20k or not prob_file_5k:
            raise SystemExit(
                f'[ ERROR ] - No se encontraron ambas ventanas (20k/5k) para el probando {prob_cnag}'
            )
        prob_coord = process_cnag_proband(prob_file_20k, prob_file_5k, proband, folder_outputs)
        
        # Procesar CNVs familiares
        rels_files = []
        for cnag_id, files in cnag_files.items():
            if cnag_id == prob_cnag:
                continue

            file_20k = files.get('20k')
            file_5k = files.get('5k')
            if not file_20k or not file_5k:
                logging.warning(f'{cnag_id}: faltan archivos de una o ambas ventanas (20k/5k). Se omite esta muestra.')
                continue
            rels_files.append((cnag_id, file_20k, file_5k))

        process_cnag_relatives(rels_files, prob_coord, cnag_to_imp, folder_outputs, output_file)
    
    elif node == 'FPGMX':
        # Aumentar el límite de tamaño de campo CSV
        set_csv_field_size()

        # Validar que hay archivos para las muestras germinales esperadas
        validate_germinal_files(pheno_dict, input_files, node)
        
        # Procesar archivos del probando
        for f in input_files:
            if f.startswith(proband) and f.endswith('full.tab.gz'):
                prob_full_file = os.path.join(folder_inputs, f)
            elif f.startswith(proband) and f.endswith('split.tab.gz'):
                prob_split_file = os.path.join(folder_inputs, f)
        prob_coord = process_fpgmx_proband(prob_full_file, prob_split_file, proband, folder_outputs)
        
        # Procesar archivos de los familiares
        rels_full_files = [
            os.path.join(folder_inputs, f) for f in input_files
            if not f.startswith(family_id) and f.endswith('full.tab.gz')
        ]
        process_fpgmx_relatives(rels_full_files, prob_coord, folder_outputs, output_file)

    elif node == 'NASERTIC':
        # Validar que hay archivos para las muestras germinales esperadas
        validate_germinal_files(pheno_dict, input_files, node)
        
        # Procesar archivo del probando
        prob_file = next(
            os.path.join(folder_inputs, f) for f in input_files
            if f.startswith(proband)
        )
        prob_coord = process_nasertic_proband(prob_file, proband, folder_outputs)
        
        # Procesar archivos de familiares
        rels_files = [
            os.path.join(folder_inputs, f) for f in input_files
            if not f.startswith(family_id)
        ]
        process_nasertic_relatives(rels_files, prob_coord, folder_outputs, output_file)

    logging.info("PREPROCESAMIENTO DE CNVs COMPLETADO")



if __name__ == '__main__':
    # folder_path = 'C:/Users/edurne.urrutia/Desktop/VUSCAN_PIPELINES/DATA/SAMPLES/CNAG/1282/'
    WORKDIR = sys.argv[1]
    folder_path = WORKDIR

    if not os.path.isdir(folder_path):
        logging.error('La ruta introducida no es válida o no existe.')
    else:
        main(folder_path)
