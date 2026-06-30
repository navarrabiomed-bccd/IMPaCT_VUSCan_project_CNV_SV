# -*- coding: utf-8 -*-
"""
CNV and SV pipeline: Preprocessing CNVs and SVs files to generate annotated
and parsed TSV files for the SNV pipeline. For use with the VUSCan families.

Last modification: April 2026
"""

import os
import csv
import logging
import pandas as pd
import sys
from natsort import natsorted
import gzip
import shutil



# -----------------------------------------------
# PARAMETERS CONFIGURATION
# -----------------------------------------------

QUAL_MIN = 30             # Minimum quality QUAL
INTERNAL_FREQ_MAX = 10    # Maximum internal frequency count
SAMPLE_PHENO_FILE = '/mnt/c/Users/bioinfor/Desktop/VUSCAN_PIPELINES/DATA/muestras_fenotipo.txt'


# -----------------------------------------------
# LOGGING SETUP
# -----------------------------------------------

logging.basicConfig(
    level=logging.DEBUG,  # Minimum level to display
    format='%(asctime)s - [ %(levelname)s ] - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)


# -----------------------------------------------
# Funciones de lectura y guardado de archivos
# -----------------------------------------------

def build_sample_id_maps(family_id, file_path=SAMPLE_PHENO_FILE):
    """Construye mapas de conversión entre IDs CNAG e IMPaCT."""
    try:
        df = pd.read_csv(SAMPLE_PHENO_FILE, sep='\t', encoding='utf-8')
    except UnicodeDecodeError:
        df = pd.read_csv(SAMPLE_PHENO_FILE, sep='\t', encoding='ISO-8859-1')
    
    pattern = f'-{family_id}-'
    df = df[df['Muestra ID'].astype(str).str.contains(pattern)]
    pheno_dct = df.set_index('Muestra ID').T.to_dict()

    imp_to_cnag, cnag_to_imp = {}, {}

    for imp_id, data in pheno_dct.items():
        cnag_id = data.get('ID_CNAG')
        if cnag_id:
            imp_to_cnag[imp_id] = cnag_id
            cnag_to_imp[cnag_id] = imp_id

    return imp_to_cnag, cnag_to_imp


def build_variant_id(row):
    """Genera un ID estandar de SV: chrom_start_end_type."""
    chrom = row['SV_chrom']
    start = int(row['SV_start'])
    end = int(row['SV_end'])
    sv_type = row['SV_type']
    return f'{chrom}_{start}_{end}_{sv_type}'


def copy_somatic_folder(folder_inputs, folder_outputs, node, variant_type):
    # Carpeta con los archivos somáticos
    if node == 'FPGMX' and variant_type == 'CNV':
        somatic_input = os.path.join(folder_inputs, 'SOMATIC_CNV')
    elif node == 'FPGMX' and variant_type == 'SV':
        somatic_input = os.path.join(folder_inputs, 'SOMATIC_SV')
    else:
        somatic_input = os.path.join(folder_inputs, 'SOMATIC')
    somatic_output = os.path.join(folder_outputs, 'SOMATIC')

    if not os.path.isdir(somatic_input):
        logging.warning(
            f'No existe la carpeta SOMATIC de entrada ({somatic_input}). '
            'Se continúa sin copiar archivos somáticos.'
        )
        return

    # Copiar carpeta completa en estos casos
    copy_entire_folder = (node == 'CNAG' or (node == 'NASERTIC' and variant_type == 'CNV'))

    extension_files = None
    if not copy_entire_folder:
        if node == 'NASERTIC' and variant_type == 'SV':
            extension_files = ('.tsv', '.tsv.gz', '.tab')

        elif node == 'FPGMX' and variant_type in {'CNV', 'SV'}:
            extension_files = ('.tsv', '.tsv.gz')

    # Copiar carpeta/archivos
    if copy_entire_folder:
        shutil.copytree(somatic_input, somatic_output, dirs_exist_ok=True)
        logging.info(f'Carpeta SOMATIC completa copiada para {node} y {variant_type}.')
        return
    
    if not extension_files:
        logging.warning(f'No hay reglas para copiar SOMATIC para {node} y {variant_type}.')
        return
    
    os.makedirs(somatic_output, exist_ok=True)

    copied = 0
    for fname in os.listdir(somatic_input):
        if fname.endswith(extension_files):
            src = os.path.join(somatic_input, fname)
            dst = os.path.join(somatic_output, fname)
            shutil.copy2(src, dst)
            copied += 1

    logging.info(f'Carpeta SOMATIC copiada ({copied} archivos) para {node} y {variant_type}.')


# -----------------------------------------------
# CNV - Funciones de procesamiento de archivos
# -----------------------------------------------

def process_cnv_cnag(file_path, output_tsv):
    cnv_dct = {}

    with gzip.open(file_path, 'rt', encoding='iso-8859-1') as f:
        reader = csv.DictReader(f, delimiter='\t')
        
        for row in reader:
            if row['Annotation_mode'] == 'full':
                # # Filtrar 'PASS' y frecuencia interna
                # if row['FILTER'] != 'PASS':
                #     continue

                # Filtar por frecuencia interna (convertir '#' a 0)
                try:
                    row['Illumina.exact.counts'] = 0 if row['Illumina.exact.counts'] == '#' else float(row['Illumina.exact.counts'])
                    if float(row['Illumina.exact.counts']) >= INTERNAL_FREQ_MAX:
                        continue
                except ValueError:
                    logging.error(f'CNV {row["AnnotSV_ID"]}: frecuencia interna no es un valor numérico ({row["Illumina.exact.counts"]})')
                    continue

                if row['Illumina.similar.counts'] == '#':
                    row['Illumina.similar.counts'] = 0

                cnv_id = build_variant_id(row)

                if cnv_id not in cnv_dct:
                    cnv_dct[cnv_id] = {k: v for k, v in row.items() if k not in {'AnnotSV_ID', None}}
    
    # Ordenar por cromosoma y posición de inicio
    sorted_cnvs = natsorted(cnv_dct.keys(), key=lambda id: (cnv_dct[id]['SV_chrom'], int(cnv_dct[id]['SV_start'])))
    header_id = ['AnnotSV_ID'] + [k for k in next(iter(cnv_dct.values())).keys() if k is not None]
    
    # Guardar CNVs en archivo TSV
    with open(output_tsv, 'w', newline='') as tsvfile:
        writer = csv.DictWriter(tsvfile, fieldnames=header_id, delimiter='\t')
        writer.writeheader()
        for id in sorted_cnvs:
            row = {'AnnotSV_ID': id, **cnv_dct[id]}
            writer.writerow(row)


def process_cnv_fpgmx(file_path, output_tsv):
    cnv_dct = {}

    # Aumentar el límite de tamaño de campo CSV
    max_size = sys.maxsize
    while True:
        try:
            csv.field_size_limit(max_size)
            break
        except OverflowError:
            max_size = int(max_size / 10)
    
    # Procesar CNVs 'full'
    with gzip.open(file_path, 'rt', encoding='iso-8859-1') as file:
        reader = csv.DictReader(file, delimiter='\t')
        header = reader.fieldnames

        for row in reader:
            if float(row['Qual']) < QUAL_MIN:
                continue
            if row['SV_type_original'] not in {'GAIN', 'LOSS'}:
                continue
            if int(row['Illumina_DRAGEN.exact.counts']) >= INTERNAL_FREQ_MAX:
                continue

            cnv_id = build_variant_id(row)
            if cnv_id not in cnv_dct:
                cnv_dct[cnv_id] = {k: v for k, v in row.items() if k != 'AnnotSV_ID'}
    
    # Ordenar por cromosoma y posición de inicio
    sorted_cnvs = natsorted(cnv_dct.keys(), key=lambda id: (cnv_dct[id]['SV_chrom'], int(cnv_dct[id]['SV_start'])))

    # Guardar CNVs en un archivo TSV
    header_id = ['AnnotSV_ID'] + [field for field in header if field != 'AnnotSV_ID']
    with open(output_tsv, 'w', newline='') as tsvfile:
        writer = csv.DictWriter(tsvfile, fieldnames=header_id, delimiter='\t')
        writer.writeheader()
        for id in sorted_cnvs:
            row = {'AnnotSV_ID': id, **cnv_dct[id]}
            writer.writerow(row)


def process_cnv_nasertic(input_file, output_tsv):
    cnv_dct = {}

    with open(input_file, newline='') as file:
        reader = csv.DictReader(file, delimiter='\t')
        header = reader.fieldnames

        for row in reader:
            if row['Annotation_mode'] == 'full':
                if row['FILTER'] != 'PASS':
                    continue
                if float(row['QUAL']) < QUAL_MIN:
                    continue
                if int(row['Illumina_DRAGEN.exact.counts']) >= INTERNAL_FREQ_MAX:
                    continue

                cnv_id = row['AnnotSV_ID']
                cnv_dct[cnv_id] = {k: v for k, v in row.items() if k != 'AnnotSV_ID'}

    # Ordenar por cromosoma y posición de inicio
    sorted_cnvs = natsorted(cnv_dct.keys(), key=lambda id: (cnv_dct[id]['SV_chrom'], int(cnv_dct[id]['SV_start'])))
    
    # Guardar CNVs en un archivo TSV
    with open(output_tsv, 'w', newline='') as tsvfile:
        writer = csv.DictWriter(tsvfile, fieldnames=header, delimiter='\t')
        writer.writeheader()
        for id in sorted_cnvs:
            row = {'AnnotSV_ID': id, **cnv_dct[id]}
            writer.writerow(row)


# -------------------------------------------
# SV - Funciones de procesamiento de archivos
# -------------------------------------------

def process_sv_cnag(file_path):
    sv_dict = {}

    # Abrir archivo de entrada según el formato
    opener = gzip.open if file_path.endswith('.gz') else open
    with opener(file_path, 'rt', encoding='iso-8859-1') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row['Annotation_mode'] == 'full':
                if float(row['QUAL']) < QUAL_MIN:
                    continue
                if row['FILTER'] != 'PASS':
                    continue
                try:
                    row['Illumina.exact.counts'] = 0 if row['Illumina.exact.counts'] == '#' else float(row['Illumina.exact.counts'])
                    if float(row['Illumina.exact.counts']) >= INTERNAL_FREQ_MAX:
                        continue
                except ValueError:
                    logging.error(f'SV {row["AnnotSV_ID"]}: frecuencia interna no es un valor numérico ({row["Illumina.exact.counts"]})')
                    continue

                if row['Illumina.similar.counts'] == '#':
                    row['Illumina.similar.counts'] = 0

                sv_id = build_variant_id(row)
                if sv_id not in sv_dict:
                    sv_dict[sv_id] = {k: v for k, v in row.items() if k not in {'AnnotSV_ID', None}}
                    
    return sv_dict


def process_sv_fpgmx(file_path):
    sv_dict = {}
   
    with gzip.open(file_path, 'rt', encoding='iso-8859-1') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if int(row['Qual']) < QUAL_MIN:
                continue
            if row['SV_type_original'] in {'GAIN', 'LOSS'}:
                continue
            if int(row['Illumina_DRAGEN.exact.counts']) >= INTERNAL_FREQ_MAX:
                continue

            sv_id = build_variant_id(row)
            sv_dict[sv_id] = {k: v for k, v in row.items() if k not in {'AnnotSV_ID', None}}

    return sv_dict


def process_sv_nasertic(file_path):
    sv_dict = {}

    with open(file_path, 'rt') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row['Annotation_mode'] == 'full':
                if float(row['QUAL']) < QUAL_MIN:
                    continue
                if row['FILTER'] != 'PASS':
                    continue
                if int(row['Illumina_DRAGEN.exact.counts']) >= INTERNAL_FREQ_MAX:
                    continue

                sv_id = build_variant_id(row)
                if sv_id not in sv_dict:
                    sv_dict[sv_id] = {k: v for k, v in row.items() if k not in {'AnnotSV_ID', None}}

    return sv_dict


def process_family_sv(folder_path, folder_results, node, family_id):
    input_files = os.listdir(folder_path)
    if not input_files:
        logging.error('No se encontraron archivos de SVs en la carpeta')
        sys.exit(1)

    _, cnag_to_imp = build_sample_id_maps(family_id) if node == 'CNAG' else ({}, {})
    processed_count = 0
    
    try:
        for file_name in input_files:
            if node == 'CNAG' and file_name.endswith(('.tab', '.tab.gz')):
                file_path = os.path.join(folder_path, file_name)

                cnag_id = file_name.split('_')[0].split('.')[0]
                sample_id = cnag_to_imp[cnag_id]

                logging.info(f'Procesando SVs de muestra: {cnag_id} / {sample_id}')

                sv_dct = process_sv_cnag(file_path)

            elif node == 'FPGMX':
                # Seleccionar archivo con CNVs 'full'
                if file_name.endswith('.tab.gz') and ('full' in file_name):
                    file_path = os.path.join(folder_path, file_name)
                    sample_id = file_name.split('.')[0]
                    logging.info(f'Procesando SVs de muestra: {sample_id}')

                    sv_dct = process_sv_fpgmx(file_path)
        
            elif node == 'NASERTIC' and file_name.endswith(('.tsv', '.tab')):
                file_path = os.path.join(folder_path, file_name)
                sample_id = file_name.split('.')[0]
                logging.info(f'Procesando SVs de muestra: {sample_id}')

                sv_dct = process_sv_nasertic(file_path)              

            # Guardar SVs procesadas en un archivo TSV
            sorted_ids = natsorted(sv_dct.keys(), key=lambda id: (sv_dct[id].get('SV_chrom', ''), int(sv_dct[id].get('SV_start', 0))))
            header = ['AnnotSV_ID'] + [k for k in next(iter(sv_dct.values())).keys() if k is not None]
            output_tsv = os.path.join(folder_results, f'{sample_id}.SVs.annotated_parsed.tsv')

            with open(output_tsv, 'w', newline='') as tsvfile:
                writer = csv.DictWriter(tsvfile, fieldnames=header, delimiter='\t', extrasaction='ignore')
                writer.writeheader()
                for id in sorted_ids:
                    row = {'AnnotSV_ID': id, **sv_dct[id]}
                    writer.writerow(row)

            processed_count += 1

        if processed_count == 0:
            logging.warning('No se encontraron archivos de SVs con los criterios esperados para procesar')

    except FileNotFoundError:
        logging.error(f'Carpeta no encontrada: {folder_path}')
        sys.exit(1)
    except Exception as e:
        logging.error(f'Error inesperado al procesar la familia {family_id}: {e}')
        sys.exit(1)


# -----------------------------------------------
# FUNCIÓN PRINCIPAL
# -----------------------------------------------
def main(folder_path):
    # Carpeta para guardar los resultados
    folder_outputs = '/mnt/c/Users/bioinfor/Desktop/VUSCAN_PIPELINES/archivos_Dido/'

    # ID familia y nodo de secuenciación
    family_id = folder_path.split('/')[-2]
    node = folder_path.split('/')[-3]

    # ---------- CNVs ----------
    # Rutas a las carpetas con los ficheros de entrada y salida
    if node == 'FPGMX':
        cnv_inputs = os.path.join(folder_path, 'INPUTS/')
    else:
        cnv_inputs = os.path.join(folder_path, 'CNV/INPUTS/')
    cnv_outputs = os.path.join(folder_outputs, family_id,'CNV/')

    # Crear la carpeta de salida si no existe
    if not os.path.exists(cnv_outputs):
        os.makedirs(cnv_outputs)
        logging.info(f'Creada carpeta de CNVs "{cnv_outputs}"')
    else:
        logging.warning(f'Carpeta de CNVs "{cnv_outputs}" ya existe')

    # Copiar carpeta 'SOMATIC' con los archivos de somático
    if node == 'FPGMX':
        copy_somatic_folder(os.path.join(folder_path, 'INPUTS/'), cnv_outputs, node, 'CNV')
    else:
        copy_somatic_folder(cnv_inputs, cnv_outputs, node, 'CNV')

    # Archivos con las CNVs de probando y familiares
    input_files = [f for f in os.listdir(cnv_inputs) if f.endswith(('.tsv', '.tab', '.tab.gz'))]
    if not input_files:
        logging.error('No se encontraron archivos de CNVs en la carpeta')  
        sys.exit(1)
    
    is_file_processed = False

    # Procesar CNVs
    if node == 'CNAG':
        # Obtener correspondencia entre IDs de CNAG
        _, cnag_to_imp = build_sample_id_maps(family_id)

        # Seleccionar CNVs - ventana 20 kbp
        for file_name in input_files:
            if file_name.endswith(('fixed.tab.gz', '_20000.CNVs.p.value.annotated.IntFreq.tab.gz')):
                file_path = os.path.join(cnv_inputs, file_name)
                cnag_id = file_name.split('_')[0].split('.')[0]
                sample_id = cnag_to_imp[cnag_id]
                logging.info(f'Procesando CNVs de muestra: {cnag_id} / {sample_id}')
                
                output_tsv = os.path.join(cnv_outputs, f'{sample_id}.CNVs.annotated_parsed.tsv')
                process_cnv_cnag(file_path, output_tsv)
                is_file_processed = True

    elif node == 'FPGMX':
        for file_name in input_files:
            if file_name.endswith('full.tab.gz'):
                file_path = os.path.join(cnv_inputs, file_name)
                sample_id = file_name.split('.')[0]
                logging.info(f'Procesando CNVs de muestra: {sample_id}')
                
                output_tsv = os.path.join(cnv_outputs, f'{sample_id}.CNVs.annotated_parsed.tsv')
                process_cnv_fpgmx(file_path, output_tsv)
                is_file_processed = True

    elif node == 'NASERTIC':
        for file_name in input_files:
            if file_name.endswith(('.tsv', '.tab')):
                file_path = os.path.join(cnv_inputs, file_name)
                sample_id = file_name.split('.')[0]
                logging.info(f'Procesando CNVs de muestra: {sample_id}')

                output_tsv = os.path.join(cnv_outputs, f'{sample_id}.CNVs.annotated_parsed.tsv')
                process_cnv_nasertic(file_path, output_tsv)
                is_file_processed = True
    
    if not is_file_processed:
        logging.error('No se encontraron archivos de CNVs con los criterios esperados para procesar')
        sys.exit(1)

    # ---------- SVs ----------
    # Rutas a las carpetas con los ficheros de entrada y salida
    if node == 'FPGMX':
        sv_inputs = os.path.join(folder_path, 'INPUTS/')
    else:
        sv_inputs = os.path.join(folder_path, 'SV/INPUTS/')
    sv_outputs = os.path.join(folder_outputs, family_id, 'SV/')

    # Crear la carpeta de salida si no existe
    if not os.path.exists(sv_outputs):
        os.makedirs(sv_outputs)
        logging.info(f'Creada carpeta de SVs "{sv_outputs}"')
    else:
        logging.warning(f'Carpeta de SVs "{sv_outputs}" ya existe')

    # Copiar carpeta 'SOMATIC' con los archivos de somático
    if node == 'FPGMX':
        copy_somatic_folder(os.path.join(folder_path, 'INPUTS/'), sv_outputs, node, 'SV')
    else:
        copy_somatic_folder(sv_inputs, sv_outputs, node, 'SV')

    # Procesar SVs de probando y familiares
    process_family_sv(sv_inputs, sv_outputs, node, family_id)


    # ---------- Comprimir carpeta con los archivos ----------
    family_path = os.path.join(folder_outputs, family_id)
    if not os.path.exists(family_path):
        logging.error(f'No existe la carpeta de la familia {family_path}')
        sys.exit(1)

    shutil.make_archive(family_path, 'zip', root_dir=folder_outputs, base_dir=family_id)
    shutil.rmtree(family_path)
    logging.info(f'Carpeta comprimida: {family_path}.zip')



if __name__ == '__main__':
    # folder_path = 'C:/Users/edurne.urrutia/Desktop/VUSCAN_PIPELINES/DATA/SAMPLES/FPGMX/3839/'
    WORKDIR = sys.argv[1]
    folder_path = f'{WORKDIR}'
    
    if not os.path.isdir(folder_path):
        logging.error('La ruta introducida no es válida o no existe')
    else:
        main(folder_path)