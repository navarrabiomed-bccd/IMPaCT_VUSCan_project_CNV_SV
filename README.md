# IMPaCT-VUSCan - CNV and SV pipelines

Repository with CNV (copy number variant) and SV (structural variant) prioritization and annotation pipelines for the IMPaCT-VUSCan project. Outputs include Excel reports and input files for the SNV pipeline.

## Repository structure

```
VUSCAN/
├── general_scripts/                         # Orchestrator and shared scripts
│   ├── run_family_cnv_sv.sh                 # Full family-level execution
│   ├── run_annotsv.sh                       # Step 0: somatic annotation with AnnotSV
│   ├── preprocess_cnv_sv_impact.py          # CNV/SV input preparation for SNV (IMPaCT cohort)
│   ├── preprocess_cnv_sv_vuscan.py          # CNV/SV input preparation for SNV (VUSCan cohort)
│   ├── execution_summary_excel.py           # Execution summary in Excel format
│   └── calculate_cnag_5kb_reference_freq.py # CNAG 5 kbp reference frequency calculation
├── CNV_PIPELINE/                            # CNV pipeline
│   ├── cnv_preprocess.py                    # Step 1: preprocessing
│   ├── run_classifycnv.sh                   # Step 2.0: ClassifyCNV scoring
│   ├── cnv_annotation_reporting.py          # Steps 2-3: annotation and report
│   ├── generate_cnv_excel.py                # Excel report generation
│   └── run_generate_cnv_excel.sh            # Excel generation launcher
├── SV_PIPELINE/                             # SV pipeline
│   └── sv_pipeline.py                       # Steps 1-3: full pipeline
└── envs/                                    # Conda environments
    ├── pipeline.yml
    ├── annotsv.yml
    └── classifycnv.yml
```

## Requirements

- Linux OS (tested on WSL/Ubuntu).
- Conda or Mamba for environment management.
- AnnotSV installed (path configurable in your environment).
- ClassifyCNV installed (path configurable in your environment).

## Environment setup

```bash
conda env create -f envs/pipeline.yml
conda env create -f envs/annotsv.yml
conda env create -f envs/classifycnv.yml
```

## Orchestrator usage

Detailed docs:

- [general_scripts/README_run_family_cnv_sv.md](general_scripts/README_run_family_cnv_sv.md)
- [general_scripts/INFORME_FINAL_WORD_run_family_cnv_sv.md](general_scripts/INFORME_FINAL_WORD_run_family_cnv_sv.md)

Run from the general_scripts directory:

```bash
cd general_scripts
bash run_family_cnv_sv.sh <NODE> <FAMILY> [FAMILY2 ...]
bash run_family_cnv_sv.sh <NODE> --all
```

Examples:

```bash
bash run_family_cnv_sv.sh FPGMX 1042
bash run_family_cnv_sv.sh NASERTIC 3022 3023
bash run_family_cnv_sv.sh CNAG --all
```

Supported nodes: FPGMX, NASERTIC, CNAG.

## Pipeline overview

```
run_family_cnv_sv.sh
│
├── run_annotsv.sh                          Step 0: somatic annotation (AnnotSV)
├── CNV_PIPELINE/cnv_preprocess.py          CNV Step 1: parsing and filtering
├── CNV_PIPELINE/run_classifycnv.sh         CNV Step 2.0: ClassifyCNV
├── CNV_PIPELINE/cnv_annotation_reporting.py CNV Steps 2-3: annotation and report
├── SV_PIPELINE/sv_pipeline.py              SV Steps 1-3: full pipeline
└── preprocess_cnv_sv_vuscan.py             Input preparation for SNV
```

## Expected input structure

```
DATA/SAMPLES/VUSCan_families/<NODE>/<FAMILY>/
├── CNV/
│   └── INPUTS/
└── SV/
    └── INPUTS/
        └── SOMATIC/    # depends on node and sample type
```

## Main outputs

- Family log: <WORKDIR>/<FAMILY>_family.log
- Prioritized CNV and SV Excel reports
- TSV files for SNV pipeline input
