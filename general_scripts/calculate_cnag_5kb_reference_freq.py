# -*- coding: utf-8 -*-
"""
Compute an internal frequency of CNVs detected in the 5 kbp window
of CNAG using all samples from the IMPaCT cohort to use as a reference.

Last updated: 2026-03-24
"""

import csv
import glob
import logging
import os
from collections import defaultdict


BASE_DIR = '/mnt/c/Users/bioinfor/Desktop/VUSCAN_PIPELINES/DATA/SAMPLES/CNAG/IMPaCT_cohort/Analizadas_v3'
OUTPUT_FILE = '/mnt/c/Users/bioinfor/Desktop/VUSCAN_PIPELINES/DATA/SAMPLES/CNAG/Frequency_VUSCan/cnag_5kb_ref_freq.tsv'


logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - [ %(levelname)s ] - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
)


def count_family_folders(base_dir):
    """Count family folders in the main directory."""
    family_count = 0

    for entry in os.scandir(base_dir):
        if entry.is_dir() and entry.name.isdigit():
            family_count += 1

    return family_count


def get_cnag_5kb_files(base_dir):
    """Search for 5 kbp CNV files in the base directory."""
    # Patterns for 5 kbp CNV files
    patterns = [os.path.join(base_dir, '[0-9]*/CNV/INPUTS/*.CNV.tsv')]

    # Search for matches
    files = set()
    for pattern in patterns:
        files.update(glob.glob(pattern))

    logging.info(f'CNV files (5 kbp) found: {len(files)}')
    return sorted(files)


def compute_reference_counts(cnv_files):
    """"Load CNV files and calculate the count of samples carrying
    each CNV and the total number of samples."""
    # CNV_ID -> set of samples where it appears
    cnv_to_samples = defaultdict(set)
    # Global set of unique samples evaluated
    total_samples = set()

    # Iterate over all detected files
    for cnv_file in cnv_files:
        sample_id = os.path.basename(cnv_file).replace('.CNV.tsv', '')
        total_samples.add(sample_id)

        with open(cnv_file, newline='') as file:
            reader = csv.DictReader(file, delimiter='\t')

            for row in reader:
                # Select 'full' CNVs
                if row.get('Annotation_mode') != 'full':
                    continue
                # CNV ID
                cnv_id = f"{row['SV_chrom']}_{row['SV_start']}_{row['SV_end']}_{row['SV_type']}"

                # Add the sample to the set of carriers for this CNV
                cnv_to_samples[cnv_id].add(sample_id)

    return cnv_to_samples, total_samples


def save_output(cnv_to_samples, total_samples, output_file):
    """"Save the count and frequency of each CNV in a TSV file."""
    if not os.path.exists(os.path.dirname(output_file)):
        logging.error(f'Output folder does not exist: {os.path.dirname(output_file)}')
        return

    # Total number of samples used
    total_n = len(total_samples)

    with open(output_file, 'w', newline='') as tsvfile:
        fieldnames = [
            'CNV_ID',
            'Reference_exact_count',
            'Reference_frequency',
            'Reference_total_samples',
            'Reference_sample_ids',
        ]
        writer = csv.DictWriter(tsvfile, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()

        for cnv_id in sorted(cnv_to_samples.keys()):
            # Get samples carrying the CNV
            sample_ids = sorted(cnv_to_samples[cnv_id])
            # Exact count of samples with the CNV
            exact_count = len(sample_ids)
            # Relative frequency
            freq = (exact_count / total_n)

            # Save CNV data
            writer.writerow({
                'CNV_ID': cnv_id,
                'Reference_exact_count': exact_count,
                'Reference_frequency': f'{freq:.6f}',
                'Reference_total_samples': total_n,
                'Reference_sample_ids': ';'.join(sample_ids),
            })

    logging.info(f'Output file generated: {output_file}')


def main():
    logging.info('Starting CNAG 5 kbp reference frequency calculation')

    # Count family folders in the main directory
    total_families = count_family_folders(BASE_DIR)
    logging.info(f'Families: {total_families}')

    # Search for 5 kbp CNV files
    cnv_files = get_cnag_5kb_files(BASE_DIR)
    if not cnv_files:
        raise FileNotFoundError(f'No *.CNV.tsv files found in: {BASE_DIR}')

    # Calculate CNV presence per sample and total samples
    cnv_to_samples, total_samples = compute_reference_counts(cnv_files)
    # Save output with reference counts and frequencies
    save_output(cnv_to_samples, total_samples, OUTPUT_FILE)

    logging.info(f'Samples: {len(total_samples)}')
    logging.info(f'Reference CNVs: {len(cnv_to_samples)}')
    logging.info('Frequency calculation completed successfully.')


if __name__ == '__main__':
    main()
