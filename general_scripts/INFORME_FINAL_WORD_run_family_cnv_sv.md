---
title: "DOCUMENTACIÓN DE LOS ALGORITMOS DE PRIORIZACIÓN DE CNVs Y SVs"
subtitle: "IMPaCT-VUSCan project"
author: "Edurne Urrutia y Estefania Huergo (Navarrabiomed)"
date: \today

# quitar indice y numeración de secciones
toc: false
numbersections: false

mainfont: "Times New Roman"
fontsize: 11pt
geometry: margin=2.5cm

header-includes:
  # bullets
  - \usepackage{enumitem}
  - \setlist[itemize]{label=\textbullet}
  # word breaks
  - \usepackage{microtype}
  - \usepackage{xurl}
  - \usepackage{fvextra}
  - \DefineVerbatimEnvironment{Highlighting}{Verbatim}{breaklines,commandchars=\\\{\}}
  - \sloppy
  # tables
  - \usepackage{booktabs}
  - \usepackage{longtable}
  - \usepackage{array}
  - \usepackage{xltabular}
  - \usepackage{multirow}
  - \usepackage{caption}
---

\thispagestyle{empty}

\newpage

\pagenumbering{gobble}
\setcounter{tocdepth}{2}
\tableofcontents

\newpage

\pagenumbering{arabic}


# 1. Introducción

En este informe se describen de forma completa los algoritmos de anotación y priorización de Copy Number Variations (CNVs) y Structural Variants (SVs) del proyecto IMPaCT-VUSCan. La priorización de variantes se basa en una serie de filtros y anotaciones funcionales desarrolladas específicamente para datos familiares. Se anotan con información de frecuencia poblacional, genes asociados a cáncer y funciones proteicas relevantes, entre otras. También se evalúa si las variantes están en el probando y en otros miembros de la familia para determinar segregación o herencia compartida.

Los pasos principales de los algoritmos son: (1) preprocesamiento, (2) priorización y (3) generación del informe. Sin embargo, cada uno tiene sus particularidades y cada nodo proporciona archivos diferentes, por lo que requieren un preprocesamiento distinto. También se incorpora la información sobre el solapamiento entre CNVs y SVs en el informe final.

El script `run_family_cnv_sv.sh` es el orquestador principal que ejecuta de forma automatizada estos algoritmos en familias. El objetivo del orquestador es estandarizar el procesamiento por familia en una sola ejecución reproducible, evitando la ejecución separada de CNVs y SVs.

De manera que su función es ejecutar para cada familia con trazabilidad por log y control de estado global:

1. Anotación de archivos somáticos con AnnotSV cuando es necesario.
2. Algoritmo de CNV completo.
3. Algoritmo de SV completo.
4. Preparación de archivos de CNV y SV para el algoritmo de SNVs.

Adicionalmente:

- Procesa familias en paralelo.
- Mueve la carpeta de análisis de cada familia completada a otra carpeta (`PENDIENTES_SUBIR`).
- Genera un log con el análisis de cada familia.


\newpage


# 2. Archivos de entrada

## 2.1. Archivos base

- Archivo PED con el árbol genealógico para cada familia: `<FAMILY_ID>.ped` (carpeta de la familia correspondiente)
- Archivo con información de muestras-fenotipo de todas las familias, que incluye la correspondencia entre IDs de IMPaCT y CNAG: `muestras_fenotipo.txt` (carpeta del proyecto)
- Archivos para anotar variantes:
    - Lista de 367 genes asociados con cáncer hereditario: `lista_CH_genes.txt` (carpeta cancer_genes/)
    - Genes asociados con cáncer de la base de datos COSMIC: `tsg_og.tsv` (carpeta cancer_genes/)
    - gnomAD: `gnomad.v4.1.cnv.non_neuro_controls_filtered.tsv` (CNVs) y `gnomad.v4.1.sv.non_neuro_controls_final.tsv` (SVs) (carpeta gnomAD/)
    - Funciones de Uniprot: `uniprot.txt` (carpeta pathway_annotation_db/)
    - Rutas de KEGG: `kegg_pathways.tsv` (carpeta pathway_annotation_db/)
    - Rutas de Reactome: `Gene2Reactome.txt` (carpeta pathway_annotation_db/)
- Glosarios para los informes Excel: `glosario_cnv.xlsx` (CNVs) y `glosario_sv.xlsx` (SVs) (carpeta glosarios/)


## 2.2. Archivos de los nodos de secuenciación

Estos son los archivos que entrega cada nodo de secuenciación para cada muestra, con las variantes anotadas con AnnotSV salvo que se indique lo contrario (información adicional en Anexo A):

- **CNVs**
    - **Muestras de ADN germinal**
        - CNAG: Dos archivos que contienen CNVs detectadas en las ventanas de 20,000 bp y 5,000 bp.
            - IMPaCT-Genómica: `*CNV.fixed.tsv` (ventana 20 kbp) y `*CNV.tsv` (ventana de 5 kbp).
            - VUSCan: `*fixed.tab.gz` y `*20000.CNVs.p.value.annotated.IntFreq.tab.gz` (ventana de 20 kbp), y `*5000.CNVs.p.value.annotated.tsv` y `*5000.CNVs.p.value.annotated.IntFreq.tab.gz` (ventana de 5 kbp).
        - FPGMX: Dos archivos que contienen CNVs y SVs, `*split_sv.full.tab.gz` ('full') y `*split_sv.split.tab.gz` ('split').
        - NASERTIC: Un único archivo `*.tab` (IMPaCT-Genómica) o `*CNVs.annotated.updated.tsv` (VUSCan).
    - **Muestras de ADN somático**
        - CNAG:
            - Tumor Only: `*TumorOnly.ann.tag.tsv`
            - Tumor vs Normal: `*Tumor_vs_Normal.ann.tag.tsv`
            - Tumor vs PON: `*Tumor_vs_Baseline.ann.tag.tsv`
        - FPGMX:
            - Tumor Only: `*_DTO_.nonref.cnv.vcf.gz` (sin anotar; el *pipeline* genera el archivo `*_DTO_.nonref.cnv.filtered_annotated.tsv.gz` filtrado por 'PASS' y anotado).
            - Tumor vs Normal: con patrón `_DTN_` y puede estar anotado (TSV) o no (VCF).
            - Tumor vs PON: con patrón `_D_` y puede estar anotado (TSV) o no (VCF).
            Para VUSCan algunas muestras de tumor pueden tener una 'L' adicional en el patrón (LDTO, LDTN o LD, respectivamente).
        - NASERTIC:
            - Tumor Only: `*CNVs.annotated.tsv` (IMPaCT-Genómica) o `*4impact.cnv.vcf.gz` (VUSCan; modificado a `*4impact-{xx}.cnv.to.vcf`; sin anotar; el *pipeline* genera el archivo `*CNVs.TO.filtered_annotated.tsv.gz` filtrado por 'PASS' y anotado).
            - Tumor vs Normal: no hay para IMPaCT-Genómica; `*CNVs.annotated.tsv` para VUSCan (modificado a `*.CNVs.TN.annotated.tsv` y comprimido a `.gz` al ejecutar el *pipeline*).
            - Tumor vs PON: `*.CNVs.PON.annotated.tsv` (IMPaCT-Genómica) o `*.cnv.vcf` (VUSCan; modificado a `*.cnv.pon.vcf`; sin anotar; el *pipeline* genera el archivo `*CNVs.PON.filtered_annotated.tsv.gz` filtrado por 'PASS' y anotado)

- **SVs**
    - **Muestras de ADN germinal**
        - CNAG: `*.SV.fixed.tsv`
        - FPGMX: Dos archivos que contienen CNVs y SVs, `*split_sv.full.tab.gz` ('full') y `*split_sv.split.tab.gz` ('split').
        - NASERTIC: `*.SVs.annotated.tab` para IMPaCT-Genómica y `*.SVs.annotated.updated.tsv` para VUSCan.
    - **Muestras de ADN somático**
        - CNAG:
            - Tumor Only: `*.tumor_sv.annotated.CandidateGenes.tsv`
            - Tumor vs Normal: `*.somatic_sv.annotated.CandidateGenes.tsv`
        - FPGMX:
            - Tumor vs PON: `*_D_sv.annotated.tsv.gz`
            - Tumor vs Normal: `*_DTN_sv.annotated.tsv.gz`
            - Tumor Only: `*_DTO_.sv.vcf.gz` (sin anotar; el *pipeline* genera `*_DTO_.sv.filtered_annotated.tsv.gz`, filtrado por 'PASS' y anotado).
        - NASERTIC:
            - Tumor Only: `*.SVs.TO.annotated.PASSfiltered.tsv.gz` (IMPaCT-Genómica) o `*.sv.vcf` (VUSCan; modificado a `*.sv.to.vcf`; sin anotar; el *pipeline* genera el archivo `*SVs.TO.filtered_annotated.tsv.gz` filtrado por 'PASS' y anotado).
            - Tumor vs Normal: `*.TN.annotated.tsv` (IMPaCT-Genómica) o `*.SVs.annotated.tsv` (VUSCan; modificado a `*.SVs.TN.annotated.tsv`).


Según el nodo, el software utilizado para generar los archivos es diferente. En la siguiente tabla se recoge un resumen de los archivos proporcionados por cada nodo.

\captionof{table}{Resumen de los archivos proporcionados por cada nodo para CNVs y SVs.}

| **Nodo** | **Variante** | **Tipo** | **Software**    | **Archivos**        |
|----------|--------------|----------|-----------------|---------------------|
| CNAG     | CNV          | WGS      | Control-FREEC 7 | Dos archivos (5 k y 20 kbp) |
|          | SCNA         | WES      | CNVkit          | Tumor Only; Tumor vs Normal; Tumor vs PON |
|          | SV           | WGS      | MANTA           | Un archivo |
|          | SV somático  | WES      | MANTA           | Tumor Only; Tumor vs Normal |
| FPGMX    | CNV / SV     | WGS      | DRAGEN          | Archivos con CNVs y SVs ('full'/'split' separado) |
|          | SCNA         | WES      | DRAGEN          | Tumor Only; Tumor vs Normal; Tumor vs PON |
|          | SV somático  | WES      | DRAGEN          | Tumor Only; Tumor vs Normal; Tumor vs PON |
| NASERTIC | CNV          | WGS      | DRAGEN          | Un archivo |
|          | SCNA         | WES      | DRAGEN          | Tumor Only*; Tumor vs Normal (VUSCan); Tumor vs PON |
|          | SV           | WGS      | DRAGEN          | Un archivo |
|          | SV somático  | WES      | DRAGEN          | Tumor Only*; Tumor vs Normal |

*Autonormalización de DRAGEN para Tumor Only. Para Tumor vs Normal y Tumor vs PON, se utiliza la muestra pareada de sangre y un pool de normales no parafinadas, respectivamente.*



Retos identificados:
- En general, no hay comparación entre nodos.
- La información de somático debe ser interpretada con precaución.
- Tumor vs Normal: Se compara WES frente a WGS. Además, para las muestras de IMPaCT-Genómica de NASERTIC no hay este tipo de archivos de CNVs somáticas.
- Tumor vs PON: Falta de un panel de normales adecuado, ya que se utiliza un pool de muestras no parafinadas.


\newpage


# 3. Priorización de CNVs

## 3.1. Descripción del proceso

El proceso de priorización de Copy Number Variations (CNVs) consiste en tres pasos principales:
1. Preprocesamiento
2. Priorización
3. Generación del informe en formato Excel

#### Paso 1. Preprocesamiento

En el primer paso, se realiza el parseo y filtrado de los archivos de entrada proporcionados por los nodos de secuenciación. Dependiendo del nodo, los archivos de entrada varían. Las CNVs germinales, anotadas con AnnotSV por los nodos, se filtran de acuerdo a los siguientes criterios:
- Calidad mínima: se eliminan variantes con una calidad QUAL inferior a 30. Para CNAG no se aplica este filtro porque no se proporciona esta información.
- Calidad por FILTER: se conservan variantes con valor 'PASS'. Para FPGMX, los archivos ya vienen filtrados por 'PASS' y no se incluyen este campo.
- Frecuencia interna exacta (campo 'Exact counts'): se eliminan variantes con un valor superior a 10.
- Filtro adicional para FPGMX: como el mismo archivo contiene CNVs y SVs, se seleccionan las CNVs que son las que tienen 'GAIN' o 'LOSS' en el campo 'SV_type_original'.

Los filtros se aplican sobre las CNVs 'full' y, cuando una pasa los filtros, se seleccionan también sus CNVs 'split' correspondientes.

Posteriormente, para CNAG se combinan las CNVs únicas de los dos archivos proporcionados (ventanas de 5,000 y 20,000 bp). Para ello, se añaden las CNVs de la ventana de 5 kbp a las de 20 kbp, eliminando solapes recíprocos (≥ 70%) para evitar duplicados. En el campo 'CNV_window' se indica la ventana de detección.

En CNAG, cuando una variante no se encuentra en la referencia usada para calcular frecuencias internas (exact, similar y other), se codifica como '#', por lo que se transforma a 0. Además, dado que CNAG calcula la frecuencia interna usando muestras de 20 kbp como referencia, se añade una frecuencia interna para la ventana de 5 kbp calculada con las 186 muestras de la cohorte de IMPaCT-Genómica (campo 'IMPaCT exact counts (5 kbp)').

En los archivos de somático, las CNVs se filtran por calidad (FILTER = 'PASS' y QUAL ≥ 30) y se anotan con AnnotSV si no lo están. Excepto para CNAG que ya aplica el filtro 'PASS' y no proporciona el campo 'QUAL', por lo que no se aplica ningún filtro adicional.

### Paso 2. Priorización
Las CNVs germinales filtradas, se priorizan mediante dos herramientas basadas en los criterios del ACMG (American College of Medical Genetics and Genomics) y ClinGen (Clinical Genome Resource) para CNVs (Riggs et al., 2019): AnnotSV y ClassifyCNV. Los criterios implementados se recogen en las Tablas \ref{table:cnv-loss} y \ref{table:cnv-gain} para CNVs de pérdida y de ganancia, respectivamente.

Dado que AnnotSV y ClassifyCNV no emplean exactamente los mismos criterios ni las mismas bases de datos, se mantienen ambas puntuaciones por separado. El campo 'Discrepancies' recoge los criterios que se aplican de forma diferente entre ambas herramientas.

\newpage


\footnotesize

\captionof{table}{Implementación de los criterios de ACMG/ClinGen para CNVs de pérdida (*abreviaturas: ACMG = American College of Medical Genetics and Genomics; CNV = Copy Number Variation; HI = haploinsuficiencia; TS = triplosensibilidad; LOF = pérdida de función*).}
\label{table:cnv-loss}

\setlength{\extrarowheight}{2pt}
\renewcommand{\arraystretch}{1.15}

\begin{xltabular}{\textwidth}{>{\raggedright\arraybackslash}X>{\centering\arraybackslash}p{1.8cm}>{\centering\arraybackslash}p{1.6cm}>{\centering\arraybackslash}p{1.5cm}>{\raggedright\arraybackslash}p{1.6cm}>{\raggedright\arraybackslash}p{1.6cm}>{\raggedright\arraybackslash}p{1.4cm}}
\hline

\multirow{2}{*}{Criterio} & 
\multicolumn{3}{c}{Puntos} & \multicolumn{2}{c}{Fuentes} & 
\multirow{2}{*}{Evaluación} \\

\cmidrule(r){2-4} \cmidrule(l){5-6}

& ACMG/ClinGen & ClassifyCNV & AnnotSV & ClassifyCNV & AnnotSV & \\
\hline

\multicolumn{7}{l}{\parbox{\textwidth}{\textbf{Sección 1: Evaluación inicial del contenido genómico}}} \\
\hline

1A. Contiene elementos codificadores de proteínas u otros elementos conocidos importantes funcionalmente &
0 & 0 & 0 & 
RefGene, FANTOM5 Enhancers, Ensembl regulatory build &
RefSeq, ENSEMBL, EnhancerAtlas 2.0 (Consultar documentación) &
Sí \\
\hline

1B. No contiene elementos codificadores ni elementos funcionales conocidos &
-0.6 & -0.6 & -0.6 & 
RefGene, FANTOM5 Enhancers, Ensembl regulatory build &
RefSeq, ENSEMBL, EnhancerAtlas 2.0 (Consultar documentación) &
Sí \\
\hline

\multicolumn{7}{l}{\parbox{\textwidth}{\textbf{Sección 2: Solapamiento con genes/regiones de HI o establecida/predicha o benignos establecidos}}} \\
\hline

2A. Solapamiento completo de gen/región genómica HI establecido &
1 & 1 & 1 & 
ClinGen, HI, and TS genes and curated regions &
 &
Sí \\
\hline

2B. Solapamiento parcial de una región genómica HI establecida\newline
- La CNV NO contiene el gen causante o la región crítica para esta región HI O\newline
- No está claro si el gen causante o la región crítica está afectada O\newline
- No se ha establecido ningún gen causal específico o región crítica para esta región HI.
 &
0 & 0 & 0 &
ClinGen HI and TS genes and curated regions &
 &
Sí \\
\hline

2C. Solapamiento parcial con el extremo 5' de un gen HI / mórbido establecido (extremo 3' no implicado)... &
 &  &  & 
refGene gene features, ClinGen HI and TS genes and curated regions &
 &  \\
\hline

2C-1. ... y la secuencia codificante está implicada &
0.90 (0.45-1.00) & 0.90 & 0.90 &  &  & Sí \\
\hline

2C-2. ... y sólo está implicada la 5' UTR &
0 (0–0.45) & 0 & 0 &  &  & Sí \\
\hline

2D. Solapamiento parcial con el extremo 3' de un gen HI mórbido establecido (extremo no implicado) ... &
 &  &  & refGene gene features, ClinGen HI and TS genes and curated regions &  & \\
\hline

2D-1. ... y solo está implicada la UTR 3' &
0 & 0 & 0 &  &  & Sí \\
\hline

2D-2. ... y solo está implicado el último exón. Se han descrito variantes patogénicas en este exón. &
0.90 (0.45-0.90) & 0.30 & 0.90 &  &  & 
AnnotSV: Sí\newline
ClassifyCNV: Parcialmente (0.30 puntos existan o no variantes patogénicas descritas) \\
\hline

2D-3. ... y sólo el último exón está implicado. No se han descrito otras variantes patogénicas en este exón. &
0.30 (0-0.45) & 0.30 & 0.45 &  &  & 
AnnotSV: Sí\newline
ClassifyCNV: Parcialmente (0.30 puntos existan o no variantes patogénicas descritas) \\
\hline

2D-4. ... e incluye otros exones además del último exón. Se espera NMD. &
0.90 (0.45-1.00) & 0.90 & 0.90 &  &  & 
Sí \\
\hline

2E. Ambos breakpoints están dentro del mismo gen (CNV intragénica) ...\newline
2E-1. ... e interrumpe el marco de lectura\newline
2E-2. ... y >=1 exón eliminado Y otros SNV/indel patogénicos en la CNV Y la variante elimina >= 10% de la proteína\newline
2E-3. ... y >=1 exón eliminado Y otros SNV/indel patogénicos en la CNV Y la variante elimina < 10% de la proteína\newline
2E-4. ... y >=1 exón eliminado Y NO SNV/indel patogénicos en la CNV Y la variante elimina > 10% de la proteína
 &
2E-1\newline
2E-2\newline
2E-3\newline
2E-4\newline & 
- PVS1 = 0.90\newline
- N/A = Sin puntos, pero continuar la evaluación\newline & 
2E-1: 0.9\newline
2E-2: 0.45\newline
2E-3: 0.3\newline
2E-4: 0.2 &  &  & 
AnnotSV: Sí\newline
ClassifyCNV: Parcialmente. Se asignan puntos a las variantes intragénicas en transcritos biológicamente relevantes que alteran el marco de lectura y se prevé que provoquen una desintegración sin sentido. No se asignan puntos a otros tipos de deleciones intragénicas. \\
\hline

2F. Completamente contenida dentro de una región CNV benigna. &
-1 & -1 & -1 & ClinGen HI and TS genes and curated regions &  & Sí \\
\hline

2G. Se solapa con una CNV benigna, pero incluye material genómico adicional. &
0 & 0 & 0 & ClinGen HI and TS genes and curated regions &  & Sí \\
\hline

2H. Dos o más predictores HI sugieren que AL MENOS UN gen del intervalo es HI. &
0.15 & 0.15 & 0.15 & DECIPHER, ExAC pLI scores &  & Sí \\
\hline

\multicolumn{7}{l}{\parbox{\textwidth}{\textbf{Sección 3: Evaluación del número de genes}}} \\
\hline

3A. 0-24 genes. &
0 & 0 & 0 & RefGene &  & Sí \\
\hline

3B. 25-34 genes. &
0.45 & 0.45 & 0.45 & RefGene &  & Sí \\
\hline

3C. 35+ genes. &
0.90 & 0.90 & 0.90 & RefGene &  & Sí \\
\hline

\multicolumn{7}{>{\raggedright\arraybackslash}p{\dimexpr\textwidth-2\tabcolsep\relax}}{\textbf{Sección 4: Evaluación detallada del contenido genómico utilizando casos de la literatura publicada, bases de datos públicas y/o datos internos del laboratorio}} \\
\hline

4A-4N. &
- & - & - & - & - & - \\
\hline

4O. Solapamiento con variación común en la población. &
-1 (0 - -1) & -1 & -1 & DGV Gold\newline
Standard Variants & gnomAD, ClinVar, ClinGen, DDD, dbVar & Sí \\
\hline

\end{xltabular}

\renewcommand{\arraystretch}{1.0}
\setlength{\extrarowheight}{0pt}


\newpage

\footnotesize

\captionof{table}{Implementación de los criterios de ACMG/ClinGen para CNVs de ganancia (*abreviaturas: ACMG = American College of Medical Genetics and Genomics; CNV = Copy Number Variation; HI = haploinsuficiencia; TS = triplosensibilidad; LOF = pérdida de función*).}
\label{table:cnv-gain}

\setlength{\extrarowheight}{2pt}
\renewcommand{\arraystretch}{1.15}

\begin{xltabular}{\textwidth}{>{\raggedright\arraybackslash}X>{\centering\arraybackslash}p{1.8cm}>{\centering\arraybackslash}p{1.6cm}>{\centering\arraybackslash}p{1.5cm}>{\raggedright\arraybackslash}p{1.6cm}>{\raggedright\arraybackslash}p{1.6cm}>{\raggedright\arraybackslash}p{1.4cm}}
\hline

\multirow{2}{*}{Criterio} &
\multicolumn{3}{c}{Puntos} & \multicolumn{2}{c}{Fuentes} &
\multirow{2}{*}{Evaluación} \\

\cmidrule(r){2-4} \cmidrule(l){5-6}

& ACMG/ClinGen & ClassifyCNV & AnnotSV & ClassifyCNV & AnnotSV & \\
\hline

\multicolumn{7}{l}{\parbox{\textwidth}{\textbf{Sección 1: Evaluación inicial del contenido genómico}}} \\
\hline

1A. Contiene elementos codificadores de proteínas u otros elementos funcionalmente importantes. &
0 & 0 & 0 &
RefGene, FANTOM5 Enhancers, Ensembl regulatory build &
RefSeq, ENSEMBL, EnhancerAtlas 2.0. Consultar documentación. &
Sí \\
\hline

1B. NO contiene elementos codificadores de proteínas ni otros elementos funcionales conocidos. &
-0.6 & -0.6 & -0.6 &
RefGene, FANTOM5 Enhancers, Ensembl regulatory build &
 &
Sí \\
\hline

\multicolumn{7}{l}{\parbox{\textwidth}{\textbf{Sección 2: Solapamiento con genes/regiones genómicas de HI establecida/predicha o benignos establecidos}}} \\
\hline

2A. Superposición completa; el gen TS o la región crítica mínima está totalmente contenida dentro de la CNV. &
1 & 1 & 1 &
ClinGen HI and TS genes and curated regions &
 &
Sí \\
\hline

2B. Solapamiento parcial de una región genómica TS:\newline
- La CNV NO contiene el gen causante o la región crítica para esta región TS O\newline
- No está claro si el gen causante o la región crítica está afectada O\newline
- No se ha establecido ningún gen causal específico o región crítica para esta región TS. &
0 & 0 & 0 &
ClinGen HI and TS genes and curated regions &
 &
Sí \\
\hline

2C. Idéntico en contenido génico a CNV de ganancia benigna. &
-1 & -1 & -1 &
ClinGen HI and TS genes and curated regions &
 &
Sí \\
\hline

2D. Menor que la CNV de ganancia benigna; los breakpoints no interrumpen genes codificantes. &
-1 & -1 & -1 &
RefGene, ClinGen HI and TS genes and curated regions &
 &
Sí \\
\hline

2E. Menor que la CNV de ganancia benigna; el/los breakpoints interrumpe(n) potencialmente un gen codificante. &
0 & 0 & 0 &
RefGene, refGene gene features, ClinGen HI and TS genes and curated regions &
 &
Sí \\
\hline

2F. Mayor que la CNV de ganancia benigna, no incluye genes codificantes adicionales. &
-1 (0 - -1.00) & -1 & -1 &
RefGene, ClinGen HI and TS genes and curated regions &
 &
Sí \\
\hline

2G. Se solapa con una CNV benigna, pero incluye material genómico adicional. &
0 & 0 & 0 &
RefGene, ClinGen HI and TS genes and curated regions &
 &
Sí \\
\hline

2H. Gen HI totalmente contenido en la CNV. &
0 & 0 & 0 &
ClinGen HI and TS genes and curated regions &
 &
Sí \\
\hline

2I. Ambos breakpoints se encuentran en el mismo gen (posible LOF):\newline
2I-1. ... e interrumpe el marco de lectura.\newline
2I-2. ... y el fenotipo del paciente es altamente específico y consistente con lo descrito para este gen HI/mórbido.\newline
2I-3. ... y el fenotipo del paciente no concuerda con lo descrito para este gen HI/mórbido. &
2I-1\newline
2I-2\newline
2I-3 &
2I-1: 0.45\newline
2I-2: 0.45\newline
2I-3: 0 &
 &
 &
AnnotSV: Sí.\newline
ClassifyCNV: No. \\
\hline

2J. Un breakpoint está dentro de un gen HI, y el fenotipo del paciente es inconsistente o desconocido. &
0 & 0 & 0 &
 &
 &
Sí \\
\hline

2K. Un breakpoint está dentro de un gen HI, y el fenotipo del paciente es altamente específico y consistente. &
0.45 & - & 0.45 &
 &
 &
AnnotSV: Sí.\newline
ClassifyCNV: No. \\
\hline

2L. Uno o ambos puntos de rotura están dentro de gen(es) sin importancia clínica. &
0 & 0 & 0 &
 &
 &
Sí \\
\hline

\multicolumn{7}{l}{\parbox{\textwidth}{\textbf{Sección 3: Evaluación del número de genes}}} \\
\hline

3A. 0-34 genes &
0 & 0 & 0 & RefGene &  & Sí \\
\hline

3B. 35-49 genes &
0.45 & 0.45 & 0.45 & RefGene &  & Sí \\
\hline

3C. 50+ genes &
0.9 & 0.9 & 0.9 & RefGene &  & Sí \\
\hline

\multicolumn{7}{>{\raggedright\arraybackslash}p{\dimexpr\textwidth-2\tabcolsep\relax}}{\textbf{Sección 4: Evaluación detallada del contenido genómico utilizando casos de la literatura publicada, bases de datos públicas y/o datos internos del laboratorio}} \\
\hline

4A-4N &
- & - & - & - & - & No \\
\hline

4O. Solapamiento con variación común en la población. &
-1 (0 - -1) & -1 &  & DGV Gold\newline
Standard Variants &  & Sí \\
\hline

\end{xltabular}

\renewcommand{\arraystretch}{1.0}
\setlength{\extrarowheight}{0pt}


Además de la priorización, se anotan las variantes con la siguiente información:
- CNVs del probando compartidas con familiares (solapamiento recíproco ≥ 70%).
- Genes asociados: relacionados con cáncer (lista de cáncer hereditario y base de datos COSMIC), dosis génica y rutas/funciones (Uniprot, KEGG y Reactome).
- Frecuencia poblacional usando la base de datos gnomAD (solapamiento recíproco ≥ 70%).
- Presencia en somático de las CNVs del probando (solapamiento recíproco ≥ 70%).
- Solapamiento con SVs para el probando (solapamiento recíproco ≥ 70%).

Finalmente, se combinan las CNVs de los familiares que no están presentes en el probando.


### Paso 3. Generación del informe

Finalmente, se genera un informe en formato Excel con el nombre del ID de la familia seguido por 'prioritized_CNVs.xlsx' que incluye las siguientes hojas:
- Proband: CNVs del probando ordenadas por puntuación de AnnotSV, después por puntuación de ClassifyCNV y, por último, por longitud. La tabla recoge coordenadas y tipo de CNV, clasificación de ClassifyCNV y AnnotSV con discrepancias entre ellas, anotación de los genes asociados (incluyendo asociación a cáncer hereditario, COSMIC, sensibles a la dosis génica, OMIM y Uniprot, KEGG y Reactome), información familiar (genotipo, estado, sexo y breakpoints), presencia en tumor, frecuencia interna, gnomAD y solapamiento con SVs.
- Proband_split: Desglose por genes de las CNVs del probando incluidas en la hoja 'Proband' con información detallada. Si una CNV no tiene ningún gen asociado, no aparece en este listado. Incluye anotación completa proporcionada por AnnotSV, como elementos reguladores cercanos, overlap, frameshift, …
- Relatives: CNVs de familiares no presentes en el probando, con estructura similar a la hoja 'Proband'.
- Relatives_split: Una hoja para cada familiar con sus CNVs desglosadas por gen (estructura similar a 'Proband_split').
- Hojas con las variantes somáticas: Tumor Only, Tumor vs PON y Tumor vs Normal para cada muestra de tumor.

Las columnas de AnnotSV y ClassifyCNV (hojas Proband y Relatives) están coloreadas atendiendo a los criterios del ACMG:
- Patogénica (P) (puntuación ≥ 0.99): rojo
- Probablemente patogénica (LP) (puntuación 0.90 a 0.98): naranja
- Variante de significado incierto (VUS) (puntuación 0.89 a -0.89): amarillo
- Probablemente benigna (LB) (puntuación -0.90 a -0.98) y benigna (B) (puntuación ≤ -0.99): sin colorear.

Nota sobre CNVs somáticas: La información de CNVs somáticas debe interpretarse con precaución. Desde NASERTIC se informó de posibles limitaciones por el panel de normales empleado. El uso de un panel de normales no parafinadas puede no ser la mejor aproximación, ya que para estudios de CNVs somáticas debería emplearse un panel de muestras tumorales embebidas en parafina, idealmente del mismo tipo de cáncer que la muestra analizada. Por lo tanto, aunque los archivos somáticos filtrados se incluyen en el informe, su interpretación queda a criterio de cada comité.


## 3.2. Implementación de los scripts

El paso 1 (preprocesamiento) se ejecuta desde el script `cnv_preprocess.py`. La primera parte del paso 2 (priorización) es la priorización con la herramienta ClassifyCNV que se ejecuta desde el script de bash `run_classifycnv.sh`. El resto del paso 2, junto con el paso 3 (generación del informe) se ejecutan desde el script `cnv_annotation_reporting.py`. Todos necesitan como argumento de entrada la ruta a la carpeta de análisis de la familia y se pueden ejecutar de forma automática desde el orquestador `run_family_cnv_sv.sh`.

Con el script `cnv_preprocess.py` se comienza creando la carpeta `INTERMEDIATES` dentro de `CNV/` para guardar los archivos intermedios generados en este paso. Después, se buscan los archivos de CNVs germinales en la carpeta `CNV/INPUTS/` (información detallada sobre los archivos de entrada en la sección 2.2) y se valida que están todos los archivos de las muestras esperadas con la información de cada familia.

A continuación, según el nodo, se procesan estos archivos filtrando por calidad (FILTER = 'PASS' y QUAL ≥ 30, excepto CNAG que no proporciona QUAL) y frecuencia interna (Illumina.exact.counts (CNAG) o Illumina_DRAGEN.exact.counts (FPGMX y NASERTIC) < 10). Los archivos de FPGMX ya están filtrados por 'PASS' y no se dispone de este campo. Además, para FPGMX hay un único archivo con CNVs y SVs, por lo que se aplica un filtro adicional para seleccionar solo las CNVs, conservando aquellas que tengan el valor 'GAIN' o 'LOSS' en el campo 'SV_type_original'.
Cuando una CNV 'full' pasa los filtros, se selecciona también su 'split'.
En CNAG, las CNVs que no se encuentran en la referencia usada para la frecuencia interna, tienen el valor '#', por lo que se cambia a 0. Además, dado que CNAG calcula la frecuencia interna usando muestras de referencia de la ventana de 20 kbp, se añade una frecuencia interna para la ventana de 5 kbp calculada con las 186 muestras de la cohorte de IMPaCT-Genómica (campo 'IMPaCT exact counts (5 kbp)').
Para el probando se guardan archivos en formato BED, TSV y pickle con las CNVs 'full' parseadas y filtradas como `<sample>.CNVs.annotated_parsed.*`. Además, las CNVs 'split' se guardan en un archivo pickle `<sample>.CNVs.annotated_parsed_split.pkl`. Para los familiares, se guardan los mismos archivos con las CNVs pero excluyendo las que solapen recíprocamente con el probando (≥ 70%), mientras que estas CNVs compartidas con el probando se guardan en un archivo pickle adicional `<family_id>.CNVs.relatives_overlap_proband.pkl`.
Los archivos BED contienen las posiciones cromosómicas de las CNVs filtradas y se emplean como entrada para la herramienta ClassifyCNV que se ejecuta en el siguiente paso.
En cuanto a los archivos de somático, se procesan en el paso 2.7.

La primera parte del paso 2 se ejecuta con el script `run_classifycnv.sh` que aplica la herramienta ClassifyCNV para priorizar las CNVs filtradas del probando y familiares usando los archivos BED generados en el paso anterior. Se generan archivos de salida con la priorización de cada muestra (`<sample>_Scoresheet.txt`) y una carpeta con la salida completa de ClassifyCNV para cada una (`<sample>_ClassifyCNV_result`) dentro de `INTERMEDIATES/ClassifyCNV`.

En el script `cnv_annotation_reporting.py` se completa el paso de priorización y se genera el informe en formato Excel, y se comienza creando la carpeta `OUTPUTS` dentro de `CNV/` para guardar los archivos generados y el informe final. Después, se ejecutan secuencialmente los siguientes subpasos:

- Paso 2.1. Comparación de la priorización entre AnnotSV y ClassifyCNV. Como archivos de entrada se utilizan las CNVs 'full' preprocesadas en el paso anterior `<sample>.CNVs.annotated_parsed.pkl` y el resultado de la priorización con ClassifyCNV `<sample>_Scoresheet.txt`. Para cada muestra, los archivos se combinan y se guarda un archivo con las CNVs del probando `<family_id>_proband_priorization.pkl` y otro con las únicas de los familiares `<family_id>_relatives_priorization.pkl`.

- Paso 2.2. Anotación de las CNVs del probando compartidas con familiares (solapamiento recíproco ≥ 70%) añadiendo información de genotipo y estado de cada familiar. Para ello, se utilizan las CNVs del probando guardadas en el paso 2.1 y se combina con las CNVs de familiares preprocesadas en el paso 1 (`<family_id>.CNVs.relatives_overlap_proband.pkl`). El resultado se guarda en el archivo `<family_id>_proband_sharedcnvs.pkl`.

- Paso 2.3. Combinación de las CNVs de los familiares que no están en el probando añadiendo información fenotípica. Se utiliza el archivo de los familiares generado en el paso 2.1 y se guarda el archivo `<family_id>_relatives_mergedcnvs.pkl`.

- Paso 2.4. Anotación de genes relacionados con cáncer (lista de cáncer hereditario y base de datos COSMIC). No se guarda ningún archivo.

- Paso 2.5. Anotación de rutas y funciones de los genes afectados (Uniprot, KEGG y Reactome). No se guarda ningún archivo.

- Paso 2.6. Anotación con frecuencia poblacional usando la base de datos gnomAD (solapamiento recíproco ≥ 70%). Se añaden estos campos con la información de gnomAD para probando y familiares:
    - gnomAD (int): Número de CNVs coincidentes
    - gnomAD_ID (str | None): ID de la mejor coincidencia (mayor frecuencia poblacional y mayor longitud)
    - gnomAD_freq (float | None): Frecuencia poblacional
    - gnomAD_length (int | None): Longitud de la CNV
Las CNVs únicas de familiares se guardan en el archivo final `<family_id>_relatives_cnvs_annotated.pkl`.

- Paso 2.7. Procesamiento de CNVs somáticas y anotación de las CNVs germinales del probando (solapamiento recíproco ≥ 70%). 
En los archivos de somático, las CNVs se filtran por calidad (FILTER = 'PASS' y QUAL ≥ 30) y se anotan con AnnotSV si no lo están. Excepto para CNAG que ya aplica el filtro 'PASS' y no proporciona el campo 'QUAL', por lo que no se aplica ningún filtro adicional.
Se guarda un archivo con las CNVs germinales del probando anotadas `<family_id>_proband_tumour.pkl` y otro con las CNVs de somático de todas las muestras de tumor `<family_id>_family_tumours.pkl`.

- Paso 2.8. Anotación con SVs en el probando (solapamiento recíproco ≥ 70%). Se guarda un archivo final con las CNVs del probando anotadas `<family_id>_proband_cnvs_annotated.pkl`.

Nota: Todos los solapamientos se calculan teniendo en cuenta que sean variantes del mismo tipo y que solapen recíprocamente al menos un 70%.


Como último paso del algoritmo, se genera el informe final en formato Excel `<family_id>_prioritized_CNVs.xlsx` en la carpeta `CNV/OUTPUTS/`. Se utilizan las CNVs 'full' anotadas de probando (hoja 'Proband') y familiares (hoja 'Relatives') con los archivos generados en los pasos 2.8 y 2.6, respectivamente. Para las hojas de SVs 'split' se utilizan los archivos generados en el paso 1 (`<sample>.CNVs.annotated_parsed_split.pkl`) anotando los genes con la lista de genes asociados a cáncer hereditario (columna solicitada por el comité de mama). Para las hojas de somático se utiliza el archivo generado en el paso 2.7 (`<family_id>_family_tumours.pkl`). Si la familia no tiene muestras tumorales, el informe se genera sin los campos de somático en la hoja 'Proband' y sin las hojas de SVs somáticas.

Al informe se añade un glosario con la descripción de cada hoja y campo para facilitar su interpretación, el cual se encuentra en la carpeta `glosarios/`. Además, se incluye un resumen de la familia y las variantes detectadas en el análisis (hoja 'Summary'), que se genera con la función `add_summary_sheet` importada del script `execution_summary_excel.py`.


\newpage


# 4. Priorización de SVs

## 4.1. Descripción del proceso

El proceso de priorización de variantes estructurales (SVs) consiste en tres pasos principales:
1. Preprocesamiento
2. Priorización
3. Generación del informe en formato Excel

### Paso 1. Preprocesamiento
En el primer paso, se realiza el parseo y filtrado de los archivos de entrada proporcionados por los nodos de secuenciación. Dependiendo del nodo, los archivos de entrada varían. Las SVs germinales, anotadas con AnnotSV por los nodos, se filtran de acuerdo a los siguientes criterios:
- Calidad mínima: se eliminan variantes con una calidad QUAL inferior a 30.
- Calidad por 'FILTER': se conservan variantes con valor 'PASS'. Para FPGMX, los archivos ya vienen filtrados por 'PASS' y no se incluyen este campo.
- Frecuencia interna exacta (campo 'Exact counts'): se eliminan variantes con un valor superior a 10.
- Filtro adicional para FPGMX: como el mismo archivo contiene CNVs y SVs, se seleccionan las SVs que son las que no tienen 'GAIN' o 'LOSS' en el campo 'SV_type_original'.

Los filtros se aplican sobre las SVs 'full' y 'split' comprobando que las mismas variantes pasan los filtros.

En CNAG, cuando una variante no se encuentra en la referencia usada para calcular frecuencias internas (exact, similar y other), se codifica como '#', por lo que se transforma a 0.

Tras el filtrado, para cada muestra se agrupan variantes que solapan recíprocamente (≥ 70%), seleccionando como variante representativa de cada grupo la que tiene mayor puntuación de AnnotSV y, en caso de empate, la de mayor longitud. Este paso es necesario ya que hay SVs que difieren en unas pocas posiciones.

En los archivos de somático, las SVs se filtran por calidad (FILTER = 'PASS' y QUAL ≥ 30) y se anotan con AnnotSV si no lo están.

### Paso 2. Priorización
La priorización de SVs se basa en las recomendaciones de la ACMG para SVs germinales (Raca et al., 2023) usando la herramienta AnnotSV. Los criterios implementados se recogen en la Tabla 3.

<INSERTAR TABLA EN WORD>

Además de la priorización, se anotan las variantes con la siguiente información:
- SVs del probando compartidas con familiares (solapamiento recíproco ≥ 70%).
- Genes asociados: relacionados con cáncer (lista de cáncer hereditario y base de datos COSMIC), rutas/funciones (Uniprot, KEGG y Reactome).
- Frecuencia poblacional usando la base de datos gnomAD (solapamiento recíproco ≥ 70%).
- Presencia en somático de las SVs del probando (solapamiento recíproco ≥ 70%).
- Solapamiento con CNVs para el probando (solapamiento recíproco ≥ 70%).

Finalmente, se combinan las SVs de los familiares que no están en el probando.


### Paso 3. Generación del informe

Finalmente, se genera un informe en formato Excel con el nombre del ID de la familia seguido por 'prioritized_SVs.xlsx' que incluye las siguientes hojas:
- Proband: SVs del probando ordenadas por puntuación de AnnotSV y, después, por longitud. La tabla recoge coordenadas y tipo de SV, clasificación AnnotSV, anotación de los genes asociados (incluyendo asociación a cáncer hereditario, COSMIC, OMIM y Uniprot, KEGG y Reactome), información familiar (genotipo, estado, sexo y breakpoints), presencia en tumor, frecuencia interna, gnomAD, coincidencia con otras SVs y solapamiento con CNVs.
- Proband_split: Desglose por genes de las SVs del probando (hoja 'Proband') con información detallada. Si una SV no tiene ningún gen asociado, no aparece en este listado. Incluye anotación completa proporcionada por AnnotSV, como elementos reguladores cercanos, overlap, frameshift, …
- Relatives: SVs de familiares no presentes en el probando, con estructura similar a la hoja 'Proband'.
- Relatives_split: Una hoja para cada familiar con sus SVs desglosadas por gen (estructura similar a 'Proband_split').
- Hojas con las variantes somáticas: Tumor Only, Tumor vs PON y Tumor vs Normal para cada muestra de tumor.

Las columnas de AnnotSV (hojas Proband y Relatives) están coloreadas atendiendo a los criterios del ACMG:
    - Patogénica (P) (puntuación ≥ 0.99): rojo
    - Probablemente patogénica (LP) (puntuación 0.90 a 0.98): naranja
    - Variante de significado incierto (VUS) (puntuación 0.89 a -0.89): amarillo
    - Probablemente benigna (LB) (puntuación -0.90 a -0.98) y benigna (B) (puntuación ≤ -0.99): sin colorear.


## 4.2. Implementación de los scripts

Todos los pasos del algoritmo (preprocesamiento, priorización y generación del informe) se ejecutan desde el script `sv_pipeline.py`, excepto la anotación de los archivos de somático usando AnnotSV cuando no lo están, que se hace desde el script `run_annotsv.sh`. Este script recibe como entrada la ruta a la carpeta de análisis de la familia y se ejecuta de forma automática desde el orquestador `run_family_cnv_sv.sh`. Si hay archivos somáticos pendientes de anotar, se comienza con este paso.

Al iniciar `sv_pipeline.py` se crean las carpetas `INTERMEDIATES` y `OUTPUTS` dentro de `SV/` para guardar los archivos generados.

En el paso 1, busca y valida los archivos germinales de entrada en `SV/INPUTS/` según el nodo (información detallada sobre los archivos de entrada en la sección 2.2). 

A continuación, según el nodo se filtran estos archivos por calidad (FILTER = 'PASS' y QUAL ≥ 30) y frecuencia interna (Illumina.exact.counts (CNAG) o Illumina_DRAGEN.exact.counts (FPGMX y NASERTIC) < 10). Los archivos de FPGMX ya están filtrados por 'PASS' y no se dispone de este campo. Además, para FPGMX hay un único archivo con CNVs y SVs, por lo que se aplica un filtro adicional para seleccionar solo las SVs, conservando aquellas que no tengan el valor 'GAIN' o 'LOSS' en el campo 'SV_type_original'.
En CNAG, las CNVs que no se encuentran en la referencia usada para la frecuencia interna, tienen el valor '#', por lo que se cambia a 0. 

Los filtros se aplican a todas las SVs y se valida la consistencia entre anotaciones full/split. Los archivos generados con las SVs 'full' procesadas se guardan en formato TSV (`<sample>.SVs.annotated_parsed.tsv`; SVs ordenadas por cromosoma e inicio) y pickle (`<sample>.SVs.annotated_parsed.pkl`). Además, las SVs 'split' se guardan en un archivo pickle (`<sample>.SVs.annotated_parsed_split.pkl`). Este proceso se realiza para cada integrante de una familia.


En el paso 2, el script ejecuta secuencialmente los subpasos de priorización y anotación:

- Paso 2.1. Filtrado de SVs solapantes del probando a partir del archivo con las SVs 'full' procesadas en el paso 1. Para ello, se selecciona la variante representativa (mayor puntuación AnnotSV y mayor longitud) y se guarda el ID del resto de SVs. Con el resultado se genera el archivo `<family_id>_proband_svs_filtered.pkl`.

- Paso 2.2. Anotación de SVs del probando compartidas con familiares (solapamiento recíproco ≥ 70%) y añade información de segregación/fenotipo. Por un lado, se guardan estas variantes anotadas en el archivo `<family_id>_proband_sharedsvs.pkl`, pero también se guardan las SVs de familiares no compartidas con el probando en el archivo `<family_id>_relatives_unique.pkl` para su posterior combinación en el paso 2.3.

- Paso 2.3. Combinación de SVs de los familiares no presentes en el probando añadiendo información de fenotipo. Se guarda el archivo `<family_id>_relatives_mergedsvs.pkl`.

- Paso 2.4. Filtrado de SVs solapantes en los familiares, seleccionando la variante representativa (mayor puntuación AnnotSV y mayor longitud). Guarda el ID del resto de SVs. Guarda el archivo `<family_id>_relatives_svs_filtered.pkl`.

- Paso 2.5. Anotación de genes relacionados con cáncer (lista de genes de cáncer hereditario y COSMIC) en probando y familiares. No se guarda ningún archivo.

- Paso 2.6. Anotación de rutas y funciones de los genes afectados (Uniprot, KEGG y Reactome). No se guarda ningún archivo.

- Paso 2.7. Anotación con frecuencia poblacional usando la base de datos gnomAD (solapamiento recíproco ≥ 70%). Se añaden estos campos con la información de gnomAD para probando y familiares:
    - gnomAD (int): Número de CNVs coincidentes
    - gnomAD_ID (str | None): ID de la mejor coincidencia (mayor frecuencia poblacional y mayor longitud)
    - gnomAD_AC (int | None): Número de alelos
    - gnomAD_AF (float | None): Frecuencia poblacional
    - gnomAD_FREQ_HOMALT (int | None): Frecuencia de homocigotos alterados
Las SVs de familiares se guardan en el archivo final `<family_id>_relatives_sv_annotated.pkl`.

- Paso 2.8. Procesamiento de SVs somáticas y anotación de las SVs germinales del probando (solapamiento recíproco ≥ 70%).
En los archivos de somático, las SVs se filtran por calidad (FILTER = 'PASS' y QUAL ≥ 30) y se anotan con AnnotSV si no lo están. Se guarda un archivo con las SVs germinales del probando anotadas `<family_id>_<family_id>_proband_tumour.pkl` y otro con las SVs de somático de todas las muestras de tumor `<family_id>_family_tumours.pkl`.

- Paso 2.9. Anotación con CNVs en el probando (solapamiento recíproco ≥ 70%). Se guarda un archivo final con las SVs del probando anotadas `<family_id>_<family_id>_proband_svs_annotated.pkl`.

Nota: Todos los solapamientos se calculan teniendo en cuenta que sean variantes del mismo tipo y que solapen recíprocamente al menos un 70%.

En el paso 3, se genera el informe final en formato Excel `<family_id>_prioritized_SVs.xlsx` en la carpeta `SV/OUTPUTS/`. Se utilizan las SVs 'full' anotadas de probando (hoja 'Proband') y familiares (hoja 'Relatives') con los archivos generados en los pasos 2.9 y 2.7, respectivamente. Para las hojas de SVs 'split' se utilizan los archivos generados en el paso 1 (`<sample>.SVs.annotated_parsed_split.pkl`) anotando los genes con la lista de genes asociados a cáncer hereditario (columna solicitada por el comité de mama). Para las hojas de somático se utiliza el archivo generado en el paso 2.8 (`<family_id>_family_tumours.pkl`). Si la familia no tiene muestras tumorales, el informe se genera sin los campos de somático en la hoja 'Proband' y sin las hojas de SVs somáticas.

Al informe se añade un glosario con la descripción de cada hoja y campo para facilitar su interpretación, el cual se encuentra en la carpeta `glosarios/`. Además, se incluye un resumen de la familia y las variantes detectadas en el análisis (hoja 'Summary'), que se genera con la función `add_summary_sheet` importada del script `execution_summary_excel.py`.


\newpage


# Nota para CNVs y SVs

Para las muestras de germinal, la información de genotipo (GT) y número de copias (CN) depende del nodo. 
En NASERTIC se proporcionan con un formato típico de genotipo en archivos generados por DRAGEN. Para CNVs, en la columna 'sample_id' se dan los valores separados por ':' y en la columna 'FORMAT' se indica el orden. En este caso el formato es GT:SM:CN:BC:GC:CT:AC:PE, donde GT es el genotipo y CN el número de copias estimado en esa región.
En SVs está en la columna con el ID de la muestra (por ejemplo, '3832-3832-4impact-01') y tiene la información GT:FT:GQ:PL:PR:SR:SB:FS:VF, donde GT también es el genotipo.
En FPGMX, el genotipo para CNVs y SVs está en la columna 'Zygosity' y el número de copias para CNVs se indica en la columna 'CopyNumber'.
En CNAG, para CNVs el número de copias se encuentra en la columna 'ControlFreeC_CopyNumber' y el genotipo en la que tiene de nombre el ID de la muestra para las de IMPaCT-Genómica, mientras que las de VUSCan no tienen información de genotipo. Para SVs el genotipo está en la columna con el ID de la muestra.


\newpage


# 5. Arquitectura del flujo

## 5.1. Script principal (orquestador)

- [run_family_cnv_sv.sh](run_family_cnv_sv.sh)

## 5.2. Scripts llamados por el orquestador

- [run_annotsv.sh](run_annotsv.sh)
- [../CNV_PIPELINE/cnv_preprocess.py](../CNV_PIPELINE/cnv_preprocess.py)
- [../CNV_PIPELINE/run_classifycnv.sh](../CNV_PIPELINE/run_classifycnv.sh)
- [../CNV_PIPELINE/cnv_annotation_reporting.py](../CNV_PIPELINE/cnv_annotation_reporting.py)
- [../SV_PIPELINE/sv_pipeline.py](../SV_PIPELINE/sv_pipeline.py)
- [preprocess_cnv_sv_vuscan.py](preprocess_cnv_sv_vuscan.py)


# 6. Requisitos y dependencias

## 6.1 Sistema

- Linux/WSL
- Bash

## 6.2 Entornos Conda

- `annotsv_env`: Preprocesar y filtrar archivos somáticos VCF con bcftools, y después anotar con AnnotSV.
- `classifycnv_env`: Priorizar CNVs con ClassifyCNV.
- `pipeline`: Ejecutar scripts de Python de CNV/SV y preparación de archivos para SNV.

Todos los entornos se generan con los YAML del repositorio en la ruta `VUSCAN/envs/`.

## 6.3 Dependencias de rutas

- El orquestador usa rutas relativas entre carpetas del repositorio.
- Debe ejecutarse desde la carpeta [VUSCAN/general_scripts](.) para evitar errores de ruta relativa.


\newpage


# 7. Manual de ejecución

## 7.1. Preparación del entorno de trabajo (solo la primera vez)
1. Se requiere sistema Linux o una máquina virtual (como WSL).

2. Instalar Conda y los entornos necesarios (ver sección 6.2).
```bash
# Instalar Conda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
conda install mamba -n base -c conda-forge

# Instalar entorno `annotsv_env`
conda env create -f ../envs/annotsv.yml -n annotsv_env
conda activate annotsv_env
git clone https://github.com/lgmgeo/AnnotSV.git
cd AnnotSV
chmod +x AnnotSV
make PREFIX=. install
make PREFIX=. install-human-annotation
conda deactivate

# Instalar entorno `classifycnv_env`
conda env create -f ../envs/classifycnv.yml -n classifycnv_env
conda activate classifycnv_env
mamba install git -c conda-forge
git clone http://github.com/Genotek/ClassifyCNV.git
cd ClassifyCNV/
chmod +x update_clingen.sh
./update_clingen.sh
# Ejemplo de uso para comprobar la instalación
python3 ClassifyCNV.py --infile YourCNVFile.bed --GenomeBuild {hg19,hg38}
ll ClassifyCNV_results/
rm -r ClassifyCNV_results/

# Instalar entorno `pipeline`
conda env create -f ../envs/pipeline.yml -n pipeline
conda activate pipeline
```

3. Clonar el repositorio de VUSCan desde GitHub y acceder a la carpeta `general_scripts`.
```bash
git clone <repository_url>  ## https://github.com/estefaniahuergo/VUSCAN.git
cd VUSCAN/general_scripts
```

4. Crear carpeta del proyecto.
```bash
mkdir <BASEDIR>
# Donde <BASEDIR> es la ruta a la carpeta del proyecto, por ejemplo: `/mnt/c/Users/bioinfor/Desktop/VUSCAN_PIPELINES/DATA/SAMPLES/VUSCan_families/`.
```

5. Guardar los archivos necesarios (ver sección 2.1):
    - Archivo correspondencia muestras-fenotipo
    ```bash
    cp <path_to_muestras_fenotipo.txt> <BASEDIR>/muestras_fenotipo.txt
    ```

    - Archivos adicionales (COSMIC, lista genes cáncer hereditario, gnomAD, funciones Uniprot, glosario)
```bash
mkdir <BASEDIR>/cancer_genes <BASEDIR>/glosarios <BASEDIR>/gnomad <BASEDIR>/pathway_annotation_db
cp <path_to_cancer_hereditary_genes_file> <BASEDIR>/cancer_genes/lista_CH_genes.txt  # Genes asociados con cáncer hereditario
cp <path_to_cosmic_file> <BASEDIR>/cancer_genes/tsg_og.tsv  # COSMIC
cp <path_to_gnomad_file> <BASEDIR>/gnomad/gnomad.v4.1.cnv.non_neuro_controls_filtered.tsv  # Base de datos gnomAD para CNVs
cp <path_to_gnomad_file> <BASEDIR>/gnomad/gnomad.v4.1.sv.non_neuro_controls_final.tsv  # Base de datos gnomAD para SVs
cp <path_to_uniprot_file> <BASEDIR>/pathway_annotation_db/uniprot.txt  # Funciones de Uniprot
cp <path_to_kegg_file> <BASEDIR>/pathway_annotation_db/kegg_pathways.tsv  # Anotación de rutas KEGG
cp <path_to_reactome_file> <BASEDIR>/pathway_annotation_db/Gene2Reactome.txt  # Anotación de rutas Reactome
cp <path_to_glossary_cnv_file> <BASEDIR>/glosarios/glosario_cnv.xlsx  # Glosario del informe de CNVs
cp <path_to_glossary_sv_file> <BASEDIR>/glosarios/glosario_sv.xlsx  # Glosario del informe de SVs
```

6. Modificar los scripts para que apunten a la carpeta de análisis y a los archivos necesarios (punto 5).


## 7.2. Preparación de los archivos de una familia

1. Crear una carpeta para guardar los archivos de las familias y otra carpeta para cada nodo si no existen.

```bash
mkdir <BASEDIR>/SAMPLES
mkdir <BASEDIR>/SAMPLES/<NODO>
```

2. Crear carpeta con el ID de la familia en la carpeta del nodo.

```bash
mkdir <WORKDIR>
# Donde <WORKDIR> es <BASEDIR>/SAMPLES/<NODO>/<FAMILY_ID>
```

3. Guardar dentro el archivo PED (revisar que la información del archivo es correcta).
```bash
cp <path_to_ped_file> <WORKDIR>
```

4. Crear los directorios para guardar los archivos de entrada según el nodo. Para CNAG y NASERTIC, las carpetas son `CNV/INPUTS` y `SV/INPUTS`, mientras que para FPGMX es una única carpeta `INPUTS` porque CNVs y SVs están en el mismo archivo.
```bash
mkdir <WORKDIR>/CNV/INPUTS <WORKDIR>/SV/INPUTS # CNAG y NASERTIC
mkdir <WORKDIR>/INPUTS # FPGMX
```

5. Guardar los archivos de germinal según el nodo.
```bash
# Para CNAG y NASERTIC:
cp <path_to_cnv_files> <WORKDIR>/CNV/INPUTS
cp <path_to_sv_files> <WORKDIR>/SV/INPUTS
# Para FPGMX
cp <path_to_cnv_files> <WORKDIR>/INPUTS
```

6. Crear las carpetas de somático y guardar los archivos de somático.
```bash
# CNAG y NASERTIC
mkdir <WORKDIR>/CNV/INPUTS/SOMATIC <WORKDIR>/SV/INPUTS/SOMATIC
cp <path_to_somatic_cnv_files> <WORKDIR>/CNV/INPUTS/SOMATIC
cp <path_to_somatic_sv_files> <WORKDIR>/SV/INPUTS/SOMATIC

# FPGMX
mkdir <WORKDIR>/INPUTS/SOMATIC_CNV <WORKDIR>/INPUTS/SOMATIC_SV
cp <path_to_somatic_cnv_sv_files> <WORKDIR>/INPUTS/SOMATIC_CNV
cp <path_to_somatic_cnv_sv_files> <WORKDIR>/INPUTS/SOMATIC_SV
```

7. Añadir la información de la familia en el archivo `muestras_fenotipo.txt` si no está. Esta información se puede encontrar en el Sharepoint de VUSCan.


## 7.3. Uso del orquestador

Uso del orquestador desde [VUSCAN/general_scripts](.) mediante la línea de comandos:

```bash
run_family_cnv_sv.sh <NODO> <FAMILIA> [FAMILIA2 ...]
run_family_cnv_sv.sh <NODO> --all
```

Parámetros:

- `NODO`: nodo de secuenciación (`FPGMX`, `NASERTIC`, `CNAG`).
- `FAMILIA`: ID de la familia a analizar (identificador numérico). Se pueden indicar varias familias separadas por espacio.
- `--all`: analizar todas las familias en la carpeta del nodo.

Ejemplos de ejecución:

- Una familia:
```bash
run_family_cnv_sv.sh FPGMX 1042
```

- Varias familias con el listado de IDs:

```bash
run_family_cnv_sv.sh FPGMX 1042 1043 1044
```

- Todas las familias en la carpeta de un nodo:

```bash
run_family_cnv_sv.sh NASERTIC --all
```

Hay establecido un límite para ejecución en paralelo de 2 familias (`MAX_PARALLEL=2`). Para ejecutar más familias a la vez, se puede modificar este parámetro en el script `run_family_cnv_sv.sh`.


## 7.4. Rutas de trabajo

- Ruta base dada al orquestador para la ejecución (`<BASEDIR>`):

```bash
/mnt/c/Users/bioinfor/Desktop/VUSCAN_PIPELINES/DATA/SAMPLES/VUSCan_families/<NODO>
```

- Carpeta de destino final de cada familia procesada:

```bash
<BASEDIR>/SAMPLES/<NODO>/PENDIENTES_SUBIR
```

- Para cada familia:
    - Carpeta de trabajo (<WORKDIR>): `<BASEDIR>/SAMPLES/<NODO>/<FAMILIA>/`
    - Ruta al archivo log: `<WORKDIR>/<FAMILIA>_family.log`

Dado que el script usa rutas relativas para varios scripts, debe ejecutarse desde [VUSCAN/general_scripts](.) para evitar errores de ruta relativa.


## 7.5. Flujo completo por pasos

### Anotación somática (AnnotSV)

Comando:

```bash
conda run --no-capture-output -n annotsv_env bash run_annotsv.sh "$NODO" "$WORKDIR"
```

Ejecución según nodo:

- CNAG: No requiere anotación somática.
- FPGMX:
    - Comprime archivos TSV de CNVs y SVs somáticos anotados, si no lo están.
    - Procesa archivos VCF de CNVs y SVs sin anotar (`<WORKDIR>/SOMATIC/*.vcf.gz`). Genera archivos filtrados por 'PASS', anotados con AnnotSV y comprimidos:
        - CNV TO: `*.filtered_annotated.tsv.gz`
        - SV TO: `*.sv.filtered_annotated.tsv.gz`
- NASERTIC: Procesa archivos de CNVs SVs somáticas sin anotar (`<WORKDIR>/CNV/INPUTS/SOMATIC/*.vcf.gz` y `<WORKDIR>/SV/INPUTS/SOMATIC/*.vcf.gz`). Genera archivos filtrados por 'PASS', anotados con AnnotSV y comprimidos:
        - CNV TO: `<SAMPLE>.CNVs.TO.annotated.PASSfiltered.tsv.gz`
        - CNV PON: `<SAMPLE>.CNVs.PON.annotated.PASSfiltered.tsv.gz`
        - SV TO: `<SAMPLE>.SVs.TO.annotated.PASSfiltered.tsv.gz`


### Algoritmo CNVs - Paso 1: Preprocesado de CNVs

Comando:

```bash
conda run --no-capture-output -n pipeline python3 ../CNV_PIPELINE/cnv_preprocess.py "$WORKDIR"
```

Como resultado, genera archivos TSV, BED y pickle en `<WORKDIR>/CNV/INTERMEDIATES/`:

- Archivos generados por muestra (probando y familiares; para familiares solo se incluyen variantes que no están en el probando):

    - <sample>.CNVs.annotated_parsed.tsv
    - <sample>.CNVs.annotated_parsed.bed
    - <sample>.CNVs.annotated_parsed.pkl
    - <sample>.CNVs.annotated_parsed_split.pkl


- Archivo adicional para los familiares:

    - `<family>.CNVs.relatives_overlap_proband.pkl` (CNVs de familiares compartidas con el probando)


### Algoritmo CNVs - Paso 2.0: Priorización con ClassifyCNV

Comando:

```bash
conda run --no-capture-output -n classifycnv_env bash ../CNV_PIPELINE/run_classifycnv.sh "$WORKDIR"
```

Entrada esperada:

- Archivos BED para probando y familiares en la ruta `<WORKDIR>/CNV/INTERMEDIATES/{sample_id}.bed`

Como resultado, genera en `<WORKDIR>/CNV/INTERMEDIATES/ClassifyCNV/`:

- `<ID>_Scoresheet.txt`: para cada muestra genera un archivo con la priorización de ClassifyCNV
- `<ID>_ClassifyCNV_result/`: carpeta con la salida completa de ClassifyCNV para la familia


### Algoritmo CNVs - Pasos 2 y 3: Anotación de CNVs y generación del informe

Comando:

```bash
conda run --no-capture-output -n pipeline python3 ../CNV_PIPELINE/cnv_annotation_reporting.py "$WORKDIR"
```

Como resultado, genera archivos pickles de priorización/anotación y un informe final Excel en `<WORKDIR>/CNV/OUTPUTS/`:

- Archivos intermedios en formato pickle:

- `<family>_proband_priorization.pkl`
- `<family>_relatives_priorization.pkl` (si hay familiares)
- `<family>_proband_sharedcnvs.pkl`
- `<family>_relatives_mergedcnvs.pkl` (si hay familiares)
- `<family>_proband_tumour.pkl`
- `<family>_family_tumours.pkl` (si hay muestras tumorales)
- `<family>_proband_cnvs_annotated.pkl`
- `<family>_relatives_cnvs_annotated.pkl` (si hay familiares)

- Informe final:

- `<family>_prioritized_CNVs.xlsx`


### Algoritmo SVs - Pasos 1-3: Pipeline completo

Comando:

```bash
conda run --no-capture-output -n pipeline python3 ../SV_PIPELINE/sv_pipeline.py "$WORKDIR"
```

Resultado:

- Genera archivos intermedios por muestra en `<WORKDIR>/SV/INTERMEDIATES/`.
- Genera archivos de priorización/anotación e informe final Excel en `<WORKDIR>/SV/OUTPUTS/`.


Archivos intermedios en `SV/INTERMEDIATES` (por muestra):

- `<sample>.SVs.annotated_parsed.tsv`
- `<sample>.SVs.annotated_parsed.pkl`
- `<sample>.SVs.annotated_parsed_split.pkl`

Pickles de salida en `SV/OUTPUTS` (familia):

- `<family>_proband_svs_filtered.pkl`
- `<family>_proband_sharedsvs.pkl`
- `<family>_relatives_unique.pkl` (si hay familiares)
- `<family>_relatives_mergedsvs.pkl` (si hay familiares)
- `<family>_relatives_svs_filtered.pkl` (si hay familiares)
- `<family>_<family>_proband_tumour.pkl`
- `<family>_family_tumours.pkl` (si hay muestras tumorales)
- `<family>_<family>_proband_svs_annotated.pkl`
- `<family>_relatives_svs_annotated.pkl` (si hay familiares)

Reporte final SV:

- `<family>_prioritized_SVs.xlsx`

### Preparación de archivos para algoritmo de SNV (opcional)

Comando:

```bash
conda run --no-capture-output -n pipeline python3 preprocess_cnv_sv_vuscan.py "$WORKDIR"
```

Este script realiza estos pasos:

- Preprocesa archivos somáticos de CNVs y SVs.
- Copia archivos somáticos de la carpeta `SOMATIC` según nodo.
- Genera una carpeta en ruta externa.
- Comprime la carpeta final de la familia a ZIP
(`/mnt/c/Users/bioinfor/Desktop/VUSCAN_PIPELINES/archivos_Dido/<family>.zip`).

Ruta de salida externa usada por el script:

- `/mnt/c/Users/bioinfor/Desktop/VUSCAN_PIPELINES/archivos_Dido/`

La lista completa de archivos generados y la estructura de salida se detallan en el Anexo A.


### Paso adicional - Comprobación de la ejecución y ruta de salida

La ejecución del orquestador genera un archivo de log por familia en `<WORKDIR>/<FAMILIA>_family.log` con salida duplicada en terminal. Se recomienda revisar este archivo para comprobar que la ejecución se ha realizado correctamente y que no hay errores o resultados inesperados.

La carpeta de análisis de cada familia procesada se mueve a `<BASEDIR>/SAMPLES/<NODO>/PENDIENTES_SUBIR`. Se aconseja revisar los informes de CNVs y SVs generados en la carpeta `OUTPUTS` antes de subirlos al SharePoint del proyecto. Una vez subidos también se aconseja mover la carpeta de la familia a otra carpeta, por ejemplo `<BASEDIR>/SAMPLES/<NODO>/PROCESADAS` para tener un control de las familias pendientes de analizar, las procesadas y las entregadas.


\newpage


# Anexo A. Archivos generados y estructura de salida

## A.1. Estructura de carpetas

```
<BASEDIR>/
├── cancer_genes/
├── glosarios/
├── gnomad/
├── pathway_annotation_db/
└── SAMPLES/
    └── <NODO>/
        ├── <FAMILIA>/  # Carpeta de trabajo de la familia (<WORKDIR>)
        │   ├── CNV/
        │   │   ├── INPUTS/
        │   │   ├── INTERMEDIATES/
        │   │   │   └── ClassifyCNV/
        │   │   └── OUTPUTS/
        │   └── SV/
        │       ├── INPUTS/
        │       ├── INTERMEDIATES/
        │       └── OUTPUTS/
        └── PENDIENTES_SUBIR/
```

Además, se genera una carpeta externa para la entrega de archivos a SNV:

```
/mnt/c/Users/bioinfor/Desktop/VUSCAN_PIPELINES/archivos_Dido/
└── <FAMILIA>.zip
```


## A.2. Archivos de entrada

Los archivos de NASERTIC se descargan desde el servidor Aspera, aunque cada uno está en una carpeta diferente. Para VUSCan están dentro de la carpeta 'VUSCAN' en estas rutas:
- Los archivos de CNVs y SVs de ADN germinal están en la carpeta de análisis de WGS para cada batch de muestras dentro de AnnotSV. Terminan en `.updated.tsv` que son los que incluyen la frecuencia interna.
- Los archivos de CNVs y SVs de Tumor vs Normal anotados están en las carpetas de análisis de WES para cada batch de muestras dentro de AnnotSV.
- El archivo de CNVs de Tumor vs PON está en la carpeta 'Somatic_analysis_WES'. En cambio el archivo de SVs de Tumor vs PON está dentro de 'tumor_only_analysis' > '{sample_id}_tumor_only'. Ninguno de los dos está anotado, por lo que se filtran por 'PASS' y se anotan con AnnotSV, y se guardan con el nombre ...
- El archivo de CNVs de Tumor vs PON está en 'tumor_only_analysis' > 'CNV_PON_tumor_VUSCAN'.


## A.3. Archivos generados

### A.3.1. Logging

Log por familia con salida duplicada en terminal y archivo log en `<WORKDIR>/<FAMILIA>_family.log`.


### A.3.2. CNV/INTERMEDIATES

Esta tabla consolida los pasos y archivos generados en todo el flujo, unificando la información que antes aparecía en la sección 7.5 con las tablas específicas de salida de los apartados A.3.2 a A.3.8.

| Carpeta de salida | Archivos generados | Contenido |
|-------------------|-------------------|-----------|
| Anotación somática | | |
| CNAG: Los archivos se entregan anotados | | |
| FPGMX (SVs): `INPUTS/SOMATIC/` | `<SAMPLE>.filtered_annotated.tsv.gz` y `<SAMPLE>.sv.filtered_annotated.tsv.gz` | SVs y CNVs filtradas por `PASS`, anotadas con AnnotSV y comprimidas |
| NASERTIC: `SV/INPUTS/SOMATIC/` | `<SAMPLE>.SVs.TO.annotated.PASSfiltered.tsv.gz` en NASERTIC; `<SAMPLE>.filtered_annotated.tsv.gz` y `<SAMPLE>.sv.filtered_annotated.tsv.gz` en FPGMX | SVs y CNVs somáticas filtradas por `PASS`, anotadas con AnnotSV y comprimidas según el nodo |

| CNV - Preprocesado (paso 1) | | |
| `CNV/INTERMEDIATES/` | `<sample>.CNVs.annotated_parsed.tsv`, `<sample>.CNVs.annotated_parsed.bed`, `<sample>.CNVs.annotated_parsed.pkl`, `<sample>.CNVs.annotated_parsed_split.pkl`, `<family>.CNVs.relatives_overlap_proband.pkl` | CNVs germinales parseadas, filtradas y preparadas para priorización |
| CNV - ClassifyCNV | `CNV/INTERMEDIATES/ClassifyCNV/` | `<ID>_Scoresheet.txt`, `<ID>_ClassifyCNV_result/` | Resultado de priorización con ClassifyCNV |
| CNV - Anotación e informe | `CNV/OUTPUTS/` | `<family>_proband_priorization.pkl`, `<family>_relatives_priorization.pkl`, `<family>_proband_sharedcnvs.pkl`, `<family>_relatives_mergedcnvs.pkl`, `<family>_proband_tumour.pkl`, `<family>_family_tumours.pkl`, `<family>_proband_cnvs_annotated.pkl`, `<family>_relatives_cnvs_annotated.pkl`, `<family>_prioritized_CNVs.xlsx` | Priorización final y reporte Excel de CNVs |
| SV - Pipeline completo | `SV/INTERMEDIATES/` y `SV/OUTPUTS/` | `<sample>.SVs.annotated_parsed.tsv`, `<sample>.SVs.annotated_parsed.pkl`, `<sample>.SVs.annotated_parsed_split.pkl`, `<family>_proband_svs_filtered.pkl`, `<family>_proband_sharedsvs.pkl`, `<family>_relatives_unique.pkl`, `<family>_relatives_mergedsvs.pkl`, `<family>_relatives_svs_filtered.pkl`, `<family>_<family>_proband_tumour.pkl`, `<family>_family_tumours.pkl`, `<family>_<family>_proband_svs_annotated.pkl`, `<family>_relatives_svs_annotated.pkl`, `<family>_prioritized_SVs.xlsx` | SVs germinales y somáticas parseadas, filtradas, anotadas y priorizadas |
| Preparación para SNV | `/mnt/c/Users/bioinfor/Desktop/VUSCAN_PIPELINES/archivos_Dido/` | `<family>.zip` | Carpeta final comprimida para entrega |


Archivos preprocesados de CNVs germinales para probando y familiares en `<WORKDIR>/CNV/INTERMEDIATES/`:

| Nombre archivo | Descripción | Relación |
|---|---|---|
| <sample>.CNVs.annotated_parsed.tsv | Archivo TSV con CNVs 'full' | Probando (todas las variantes)\newline
Familiares (no compartidas con el probando) |
| <sample>.CNVs.annotated_parsed.bed | Archivo BED con CNVs 'full' para ClassifyCNV | Probandos |
| <sample>.CNVs.annotated_parsed.pkl | Archivo pickle con CNVs 'full' | Probandos |
| <sample>.CNVs.annotated_parsed_split.pkl | Archivo pickle con CNV 'split' | Probandos |
| <family>.CNVs.relatives_overlap_proband.pkl | Archivo pickle con CNVs 'full' de familiares compartidas con el probando |

### A.3.3. CNV/INTERMEDIATES/ClassifyCNV

| Nombre archivo | Descripción |
|---|---|
| <ID>_Scoresheet.txt | Resultado de priorización con ClassifyCNV |
| <ID>_ClassifyCNV_result/ | Carpeta de resultados completos de ClassifyCNV |

### A.3.4. CNV/OUTPUTS

| Nombre archivo | Descripción |
|---|---|
| <family>_proband_priorization.pkl | Priorización del probando |
| <family>_relatives_priorization.pkl | Priorización de familiares |
| <family>_proband_sharedcnvs.pkl | CNV compartidas en probando |
| <family>_relatives_mergedcnvs.pkl | CNV familiares fusionadas |
| <family>_proband_tumour.pkl | Probando anotado con somático |
| <family>_family_tumours.pkl | Variantes de tumores familiares |
| <family>_proband_cnvs_annotated.pkl | CNV probando final anotado |
| <family>_relatives_cnvs_annotated.pkl | CNV familiares anotadas |
| <family>_prioritized_CNVs.xlsx | Informe final de CNV |

### A.3.5. SV/INTERMEDIATES

| Patrón de archivo | Descripción |
|---|---|
| <sample>.SVs.annotated_parsed.tsv | Tabla SV parseada |
| <sample>.SVs.annotated_parsed.pkl | SV full en pickle |
| <sample>.SVs.annotated_parsed_split.pkl | SV split en pickle |

### A.3.6. SV/OUTPUTS

| Patrón de archivo | Descripción |
|---|---|
| <family>_proband_svs_filtered.pkl | SV probando filtradas |
| <family>_proband_sharedsvs.pkl | SV probando compartidas |
| <family>_relatives_unique.pkl | SV únicas de familiares |
| <family>_relatives_mergedsvs.pkl | SV fusionadas de familiares |
| <family>_relatives_svs_filtered.pkl | SV familiares filtradas |
| <family>_<family>_proband_tumour.pkl | Probando anotado con SV somáticas |
| <family>_family_tumours.pkl | Somático familiar |
| <family>_<family>_proband_svs_annotated.pkl | SV probando anotadas finales |
| <family>_relatives_svs_annotated.pkl | SV familiares anotadas |
| <family>_prioritized_SVs.xlsx | Informe final de SV |


### A.3.7. Salida externa para SNV

Ruta base:

- /mnt/c/Users/bioinfor/Desktop/VUSCAN_PIPELINES/archivos_Dido/

Resultado:

| Ruta | Descripción |
|---|---|
| .../archivos_Dido/<family>.zip | Entrega final comprimida |


### A.3.8. Movimiento final de carpeta familiar

- Origen: <BASEDIR>/SAMPLES/<NODO>/<FAMILIA>
- Destino: <BASEDIR>/SAMPLES/<NODO>/PENDIENTES_SUBIR/<FAMILIA>
- Si existe destino: <FAMILIA>_YYYYmmdd_HHMMSS




\newpage

<!-- NOTA: La siguiente sección no se incluye en el PDF generado -->

<!-- 
# Conversión a PDF

Convertir este documento a PDF usando `pandoc` con `xelatex` desde `VUSCAN/general_scripts`:

```bash
pandoc INFORME_FINAL_WORD_run_family_cnv_sv.md \
    --from markdown \
    --pdf-engine=xelatex \
    --output /mnt/c/Users/bioinfor/Desktop/Documentacion_VUSCAN_CNVs_SVs.pdf
```
-->

