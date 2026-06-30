# -*- coding: utf-8 -*-
"""
CNV pipeline - Step 3 (Excel report generation). Pipeline step for
regenerating the final Excel report from existing intermediate files.

This script is for reruns where the intermediate files are already
available and only the final Excel report must be regenerated. The
report suffix is configurable.

NOTE: It was used in the CNAG families to add internal frequency to
the 5 kbp window computed using the IMPaCT cohort.The final output
was stored as `*_prioritized_CNVs_<YYYYMMDD>_5kbp_freq.xlsx`.
"""

import logging
import os
import shutil
import sys
from datetime import datetime


# Import the existing utilities from the main CNV reporting pipeline
from cnv_annotation_reporting import (
    SAMPLE_PHENO_FILE,
    get_proband_id,
    load_phenotype_data,
    parse_ped_file,
    run_generate_excel_report,
    split_analysis_folder,
)


logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - [ %(levelname)s ] - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
)


DEFAULT_REPORT_SUFFIX = '_5kbp_freq'


def main(folder_path, report_suffix=DEFAULT_REPORT_SUFFIX):
    """
    Generate the final CNV Excel report for one family.

    Args:
        folder_path (str): Path to the folder with intermediate files.
        report_suffix (str): Optional suffix for the report name.
    """
    _, family_id = split_analysis_folder(folder_path)
    output_folder = os.path.join(folder_path, 'OUTPUTS')
    logging.info(f'Generating final Excel report for family {family_id}')

    # Load phenotypes, PED, and proband ID
    pheno_dict = load_phenotype_data(folder_path, SAMPLE_PHENO_FILE)
    ped_dict = parse_ped_file(folder_path)
    proband = get_proband_id(folder_path, pheno_dict)

    report_date = datetime.now().strftime('%Y%m%d')
    default_report = os.path.join(
        output_folder,
        f'{family_id}_prioritized_CNVs_{report_date}.xlsx'
    )
    report_existed_before = os.path.exists(default_report)

    # Generate the report using the imported CNV function
    run_generate_excel_report(folder_path, output_folder, pheno_dict, proband, ped_dict)

    # If a previous report exists, keep a suffixed copy. Otherwise,
    # rename the generated file to avoid duplicates.
    final_report = os.path.join(
        output_folder,
        f'{family_id}_prioritized_CNVs_{report_date}{report_suffix}.xlsx'
    )

    if report_existed_before:
        shutil.copy2(default_report, final_report)
        logging.info(f'Excel report copied as: {final_report}')
    else:
        os.replace(default_report, final_report)
        logging.info(f'Excel report renamed as: {final_report}')

    logging.info('Final Excel report generated successfully.')


if __name__ == '__main__':
    WORKDIR = sys.argv[1]
    folder_path = f'{WORKDIR}CNV/'
    
    if not os.path.isdir(folder_path):
        logging.error('The provided path is not valid or does not exist.')
    else:
        main(folder_path)
