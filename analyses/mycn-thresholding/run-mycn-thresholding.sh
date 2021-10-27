#!/bin/bash
#
# This shell script runs the analysis module mycn-thresholding for assessing 
# the MYCN amplification status in OPenPedCan Neubroblastoma GMFK and TARGET 
# WGS/WXS data to determine the copy number and consensus overlap coverage 
# threshold
#
# Eric Wafula

set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
# copied from the run_in_ci.sh file at
# <https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/scripts/>
script_directory="$(perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# create the results directory
mkdir -p results

# Histologies file path 
histology_file=../../data/histologies.tsv

# CNV consensus file path
cnv_file=../../data/consensus_wgs_plus_cnvkit_wxs.tsv.gz

# CNVkit SEG file path
cnvkit_file=../../data/cnv-cnvkit.seg.gz

# Control-FREEC CNV file path
controlfreec_file=../../data/cnv-controlfreec.tsv.gz

# MantaSV SV file path
mantasv_file=../../data/sv-manta.tsv.gz

# MYCN focal SNParray file
mycn_focal=input/mycn_tumor_focal_SNParray.tsv

# MYCN paired segment SNParray file
mycn_paired_segment=input/mycn_tumor_paired_segment_SNParray.tsv

# MYCN single segment SNParray file
mycn_single_segment=input/mycn_tumor_single_segment_SNParray.tsv

# Run the first script in this module that identifies discordant samples between
# consensus CNV call and clinical data, identifies MYCN cytogenetic regions in
# CNVkit, Control-FREEC, and MantaSV analysis results, and computes classification
# metrics and overlap coverage between callers for discordant samples for determining
# the threshold for MYCN amplification status in the consensus CNV analysis.
echo 'Run 01-mycn-consensus-thresholding.py...'
python3 01-mycn-consensus-thresholding.py $histology_file $cnv_file $cnvkit_file $controlfreec_file $mantasv_file

echo 'Done 01-mycn-consensus-thresholding.py...'

echo 'Run 02-mycn-consensus-clinical-discordant.Rmd...'
# Run notebook that plots chromosome 2 segment mean from CNVkit anaysis
# for the MYCN cytogenetic region to illustrate the MYCN CNV status
# samples with discordance bewteen consensus call and clinical data
Rscript -e "rmarkdown::render('02-mycn-consensus-clinical-discordant.Rmd', clean = TRUE)"

echo 'Done 02-mycn-consensus-clinical-discordant.Rmd...'

# Run the third script in this module that computes classification metrics at mutiple copy number 
# ratio cutoffs between focal/paired segment/single segment SNParray calls and the clinical MYCN 
# amplification status in OpenPedCan histologies file to aid in determining a suitable threshold 
# for MYCN amplification calls
echo 'Run 03-mycn-SNParray-thresholding.py...'
python3 03-mycn-SNParray-thresholding.py $mycn_focal $mycn_paired_segment $mycn_single_segment $histology_file

echo 'Done 03-mycn-SNParray-thresholding.py...'


