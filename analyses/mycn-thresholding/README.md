## Assessing the status of the MYCN gene amplification for thesholding
**Module author:** Eric Wafula

## Purpose
The purpose of this analysis module is to access the status of the MYCN gene amplification status between the OpenPedCan Neuroblastoma`GMKF` and `TARGET` consensus copy number variants (CNV) calls from `WSX/WGS` data and clinical information in the master histologies file. The `mycn-thresholding` module performs the following analyses to produce results required for consensus CNV thresholding: 
1. Identifies discordant samples in MYCN amplification status between the consensus CNV calls and the clinical information.
2. Extracts segments in the MYCN cytogenic region (`2p24.3`) reported in the `CNVkit`, `Control-FREEC`, and `MantaSV` variant callers analysis results.
3. Computes copy number classification metrics using the clinical information in the v8 histologies file as the ground truth compared to the inferred consensus CNV calls.
4. Estimates overlap coverage between variant callers for samples that are discordant in MYCN amplification status between consensus CNV calls and clinical information.
5. Plots chromosome 2 segment mean for the MYCN cytogenic region (`2p24.3`) reported in the CNVkit analysis results for discordant samples to illustrate the MYCN copy number status (`deletion`, `amplification`, or `neutral`). 


## Input files
- **Histologies file**
`../../data/histologies.tsv`
- **Consensus CNV file**
`../../data/consensus_wgs_plus_cnvkit_wxs.tsv.gz`
- **CNVkit SEG file**
`../../data/cnv-cnvkit.seg.gz`
- **Control-FREEC CNV file**
`../../data/cnv-controlfreec.tsv.gz`
- **MantaSV SV file**
`../../data/sv-manta.tsv.gz`

## Analysis scripts

#### `run-mycn-thresholding.sh`
This is a wrapper bash script for setting input file paths for the main anlysis script. All file paths set in this script are based on root directory of module. Therefore, the script should always be run from the root directory, `OpenPedCan-analysis/analyses/mycn-thresholding/`


Usage:
```bash
bash run-mycn-thresholding.sh
```

#### `01-mycn-thresholding.py`
Python script to create summary results to enable consensus copy number variation (CNV) thresholding assessment for Neuroblastoma MYCN gene

Usage:
```
python3 01-mycn-thresholding.py HISTOLOGY_FILE CONSENSUS_CNV_FILE CNVKIT_CNV_FILE
                                    CONTROLFREEC_CNV_FILE MANTASV_SV_FILE
```

Parameter Options:
```
usage: 01-mycn-thresholding.py [-h] [-v] HISTOLOGY_FILE CONSENSUS_CNV_FILE CNVKIT_CNV_FILE
                                            CONTROLFREEC_CNV_FILE MANTASV_SV_FILE

The 01-mycn-thresholding.py script creates summary results to enable consensus copy number
variation (CNV) thresholding assessment for Neuroblastoma MYCN gene.

positional arguments:
  HISTOLOGY_FILE        OPenPedCan histology file (histologies.tsv)
                        
  CONSENSUS_CNV_FILE    OPenPedCan consensus CNV file (consensus_wgs_plus_cnvkit_wxs.tsv.gz)
                        
  CNVKIT_CNV_FILE       OPenPedCan CNVKit CNV file (cnv-cnvkit.seg.gz)
                        
  CONTROLFREEC_CNV_FILE
                        OPenPedCan Control-FREEC CNV file (cnv-controlfreec.tsv.gz)
                        
  MANTASV_SV_FILE       OPenPedCan mantaSV SV file (sv-manta.tsv.gz)
                        

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         Print the current 01-mycn-thresholding.py version and exit
```

#### `02-mycn-cnv-clinical-discordant.Rmd`
CNVkit amplification status notebook plots for the MYCN region `2p24.3`

Usage:
```
Rscript -e "rmarkdown::render('02-mycn-cnv-clinical-non-matching.Rmd', clean = TRUE)"
```

## Results

- **Summary table Neuroblastoma TARGET and GMKF samples for both consensus CNV and clinical MYCN amplification status**
`results/mycn_cnv_vs_clinical_status.tsv`
- **Neuroblastoma TARGET and GMKF samples with discordant MYCN amplification status between the consensus CNV calls and the clinical information, including percentage overlap coverage  between the variant callers**
`results/mycn_cnv_vs_clinical_status.tsv`
- **Copy number classification metrics of MYCN gene amplification status, including true/false positives, true/false negative, accuracy, true negative rate (sensitivity), false postive rate (specificity), and the sum of true negative and false postive rates**
`results/mycn_cnv_vs_clinical_status_classification_metrics.tsv`
- **CNVkit analysis MYCN cytogenic region (`2p24.3`) base alignment genomic segments**
`results/cnvkit_mycn_cnv_vs_clinical_diff_status.tsv`
- **Control-FREEC analysis MYCN cytogenic region (`2p24.3`) base alignment genomic segments**
`results/ccontrolfreec_mycn_cnv_vs_clinical_diff_status.tsv
- **MantaSV analysis MYCN cytogenic region (`2p24.3`) base alignment genomic segments** 
`results/mantaSV_mycn_cnv_vs_clinical_diff_status.tsv`
- **CNVkit analysisc copy number notebook plots for the MYCN cytogenic region (`2p24.3`)** 
`02-mycn-cnv-clinical-discordant.html`
