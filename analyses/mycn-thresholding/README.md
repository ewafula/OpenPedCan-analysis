## Assessing the status of the MYCN gene amplification for thesholding
**Module author:** Eric Wafula

## Purpose
The purpose of this analysis module is to access the status of the MYCN gene amplification status between the OpenPedCan Neuroblastoma`GMKF` and `TARGET` consensus copy number variants (CNV) calls from `WSX/WGS` data and `TARGET` SNP array (focal and segments) with clinical information in the master histologies file. The `mycn-thresholding` module performs the following analyses to produce results required for consensus CNV thresholding:
#### Consensus vs clinical
1. Identifies discordant samples in MYCN amplification status between the consensus CNV calls and the clinical information.
2. Extracts segments in the MYCN cytogenic region (`2p24.3`) reported in the `CNVkit`, `Control-FREEC`, and `MantaSV` variant callers analysis results.
3. Computes copy number classification metrics using the Neuroblastoma clinical MYCN amplification status information in the histologies file as the ground truth compared to the inferred consensus CNV calls.
4. Estimates overlap coverage between variant callers for samples that are discordant in MYCN amplification status between consensus CNV calls and clinical information.
5. Plots chromosome 2 segment mean for the MYCN cytogenic region (`2p24.3`) reported in the CNVkit analysis results for discordant samples to illustrate the MYCN copy number status (`deletion`, `amplification`, or `neutral`). Adapted from the AlexsLemonade/OpenPBTA-analysis [molecular-subtyping-embryonal module](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/analyses/molecular-subtyping-embryonal/03-clean-c19mc-data.Rmd). 
#### SNP arrays vs clinical
1. Extracted Neuroblastoma tumor MYCN gene copy number calls (that intersect with `hg19 genecode v26lift37` annotations) from hundreds of TARGET SNP arrays (`focal` and `segments`) sample files in the Maris Lab folder on the `Respublica HPC server` using the scripts available in the utils folder of this module. Excluded sample without any calls.
2. Using the `copy number ratio` cutoffs of `0.1 - 2.0` to score MYCN gene amplification, computed classification metrics using the clinical Neuroblastoma MYCN amplification status information in the histologies file as the ground truth compared to the inferred SNP arrays calls. TARGET WXS samples from the histologies file included those that intersect with samples in SNP arrays and have known MYCN amplification status.

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

- **MYCN focal SNP array file**

  `input/mycn_tumor_focal_SNParray.tsv`

- **MYCN paired segment SNP array file**

  `.input/mycn_tumor_paired_segment_SNParray.tsv`
  
- **MYCN single segment SNP array file**

  `input/mycn_tumor_single_segment_SNParray.tsv`


## Analysis scripts

#### `run-mycn-thresholding.sh`
This is a wrapper bash script for setting input file paths for the main anlysis script. All file paths set in this script are based on root directory of module. Therefore, the script should always be run from the root directory, `OpenPedCan-analysis/analyses/mycn-thresholding/`


Usage:
```bash
bash run-mycn-thresholding.sh
```

#### `01-mycn-consensus-thresholding.py`
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

The 01-mycn-consensus-thresholding.py script creates summary results to enable consensus copy number
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
  -v, --version         Print the current 01-mycn-consensus-thresholding.py version and exit
```

#### `02-mycn-consensus-clinical-discordant.Rmd`
CNVkit amplification status notebook plots for the MYCN region `2p24.3`

Usage:
```
Rscript -e "rmarkdown::render('02-mycn-consensus-clinical-discordant.Rmd', clean = TRUE)"
```

#### `03-mycn-SNParray-thresholding.py`
Python script to create ssummary results to enable focal and segment SNP arrays copy number variation (CNV) amplification thresholding assessment for Neuroblastoma MYCN gene.

Usage:
```
python3 03-mycn-SNParray-thresholding.py FOCAL_SNPARRAY_FILE PAIRED_SEGMENT_SNPARRAY_FILE
                                    SINGLE_SEGMENT_SNPARRAY_FILE HISTOLOGY_FILE
```

```
usage: 03-mycn-SNParray-thresholding.py [-h] [-v]
                                        FOCAL_SNPARRAY_FILE PAIRED_SEGMENT_SNPARRAY_FILE SINGLE_SEGMENT_SNPARRAY_FILE HISTOLOGY_FILE

The 03-mycn-SNParray-thresholding.py script creates summary results to enable focal and segment SNP arrays copy number variation (CNV) amplification thresholding assessment for Neuroblastoma MYCN gene.

positional arguments:
  FOCAL_SNPARRAY_FILE             TARGET MYCN gene focal SNP array file
                        
  PAIRED_SEGMENT_SNPARRAY_FILE    TARGET MYCN gene paired segmnent array file
                        
  SINGLE_SEGMENT_SNPARRAY_FILE    TARGET MYCN gene single segmnent array file
                        
  HISTOLOGY_FILE                  OPenPedCan histology file (histologies.tsv)
                        

optional arguments:
  -h, --help                     show this help message and exit
  -v, --version                  Print the current 03-mycn-SNParray-thresholding.py version and exit
```

## Utility scripts
Python scripts to extract Neuroblastoma tumor MYCN gene copy number calls  from hundreds of TARGET `focal` and `segments`SNP arrays sample files available in the Maris Lab folder on the `Respublica HPC server`. The scripts should be executed in working directory that contains 1) the sub-folder with sample focal SNP array files (`hg19-files`), 2) the paired segment SNP array file (`nexus_tumor_paired_824_geo.seg`) , 3) the paired segment SNP array file (`nexus_tumor_single_924_geo.seg`), 4) the segment SNP array metadata file (`metadata.txt`; requires cleaning to remove embedded carriage return white space characters using this command: `perl -p -e 's/\r//g' metadata.txt > metadata.txt.clean`), and 5) hg19 genecode GTF file (gmkf_paired_seg_array_intersect_histologies.tsv) for referencing MYCN gene coordinates. 

Usage:
```
python3 parse_mycn_focal_SNParray.py
python3 parse_mycn_segment_SNParray.py
```


## Results

- **Summary table Neuroblastoma TARGET and GMKF samples for both consensus CNV and clinical MYCN amplification status**

  `results/mycn_consensus_vs_clinical_status.tsv`

- **Neuroblastoma TARGET and GMKF samples with discordant MYCN amplification status between the consensus CNV calls and the clinical information, including percentage overlap coverage  between the variant callers**

  `results/mycn_consensus_vs_clinical_diff_status.tsv`

- **Consensus copy number classification metrics of MYCN gene amplification status, including true/false positives, true/false negative, accuracy, true negative rate (sensitivity), false postive rate (specificity), the sum of true negative and false postive rates, and ROC AUC score**

  `results/mycn_consensus_vs_clinical_metrics.tsv`

- **CNVkit, Control-FREEC, and MantaSV analysis MYCN cytogenic region (`2p24.3`) base alignment genomic segments**

  - `results/cnvkit_mycn_consensus_vs_clinical_diff_status.tsv`
  - `results/ccontrolfreec_mycn_consensus_vs_clinical_diff_status.tsv`
  - `results/mantaSV_mycn_consensus_vs_clinical_diff_status.tsv`

- **CNVkit analysis copy number notebook plots for the MYCN cytogenic region (`2p24.3`)**

  `02-mycn-consensus-clinical-discordant.html`

- **Neuroblastoma TARGET amplification status for both the SNP arrays calls and the clinical information**

  - `mycn_focal_SNParray_vs_clinical_status.tsv`
  - `mycn_paired_segment_SNParray_vs_clinical_status.tsv`
  - `mycn_single_segment_SNParray_vs_clinical_status.tsv`

- **SNP arrays copy number classification metrics of MYCN gene amplification status, including true/false positives, true/false negative, accuracy, true negative rate (sensitivity), false postive rate (specificity), the sum of true negative and false postive rates, and ROC AUC score**

  - `mycn_focal_SNParray_vs_clinical_metrics.tsv`
  - `mycn_paired_segment_SNParray_vs_clinical_metrics.tsv`
  - `mycn_single_segment_SNParray_vs_clinical_metrics.tsv`

