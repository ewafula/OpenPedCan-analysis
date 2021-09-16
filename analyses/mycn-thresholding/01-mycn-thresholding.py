#!/usr/bin/env python


"""
01-mycn-thresholding.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Functions to assess Neuroblastoma MYCN gene consensus copy number variation (CNV) thresholding for GMFK and TARGET cohorts for OPenPedCan analyses modules 
"""


__author__ = ('Eric Wafula (wafulae@chop.edu)')
__version__ = '1.0'
__date__ = '31 August 2021'


import os
import sys
import math
import argparse
import subprocess
import numpy as np
import pandas as pd
from collections import OrderedDict
from sklearn.metrics import roc_auc_score
from sklearn.metrics import confusion_matrix

def read_parameters():
     p = argparse.ArgumentParser(description=("The 01-mycn-thresholding.py script creates summary results to enable consensus copy number variation (CNV) thresholding assessment for Neuroblastoma MYCN gene."), formatter_class=argparse.RawTextHelpFormatter)
     p.add_argument('HISTOLOGY_FILE', type=str, default=None, help="OPenPedCan histology file (histologies.tsv)\n\n")
     p.add_argument('CONSENSUS_CNV_FILE', type=str, default=None, help="OPenPedCan consensus CNV file (consensus_wgs_plus_cnvkit_wxs.tsv.gz)\n\n")
     p.add_argument('CNVKIT_CNV_FILE', type=str, default=None, help="OPenPedCan CNVKit CNV file (cnv-cnvkit.seg.gz)\n\n")
     p.add_argument('CONTROLFREEC_CNV_FILE', type=str, default=None, help="OPenPedCan Control-FREEC CNV file (cnv-controlfreec.tsv.gz)\n\n")
     p.add_argument('MANTASV_SV_FILE', type=str, default=None, help="OPenPedCan mantaSV SV file (sv-manta.tsv.gz)\n\n")
     p.add_argument('-v', '--version', action='version', version="01-mycn-thresholding.py version {} ({})".format(__version__, __date__), help="Print the current 01-mycn-thresholding.py version and exit\n\n")
     return p.parse_args()


def get_mycn_data(consensus_cnv_file, histologies_file):
    """
    Get CNV calls for MYCN genes for each Neuroblastoma sample in both the consensus CNV files and the clinical histologies files
    Parameter: OpendeCan consensus file with CNV calls
    Parameter: OpendeCan clinical histologies file
    Returns: a Pandas dataframe with Neuroblastoma samples and their associated MYCN gene status 
    """

    # Load consensus CNV file into a pandas dataframe and filter to retain only rows with MYCN genes and the relevant columns
    cnv_df = pd.read_csv(consensus_cnv_file, sep="\t", dtype=str)
    cnv_df = cnv_df[cnv_df.gene_symbol == "MYCN"]
    cnv_df = cnv_df[["biospecimen_id", "copy_number", "status"]]
    
    # Rename columns and update MYCN status values to either amplified or not amplified to match the `molecular_subtype` column in the histologies files
    cnv_df = cnv_df.rename(columns={"biospecimen_id": "Kids_First_Biospecimen_ID", "status": "cnv_status"}, errors="raise")
    cnv_df.loc[cnv_df.cnv_status != "amplification", "cnv_status"] = "MYCN non-amp"
    cnv_df.loc[cnv_df.cnv_status == "amplification", "cnv_status"] = "MYCN amp"
    
    # Add the MYCN clinical  and cohort columns to the dataframe and update with MYCN molecular subtype values from corresponding sample ID values in the histologies files
    cnv_df = pd.concat([cnv_df.iloc[:, :3], pd.DataFrame(columns=["clinical_status", "cohort", "experimental_strategy"], index=cnv_df.index).fillna("NA"), cnv_df.iloc[:, 3:]], axis=1).fillna("NA")
    histology_df = pd.read_csv(histologies_file, sep="\t", na_filter=False, dtype=str)
    for row in histology_df.itertuples(index=False):
        if row.Kids_First_Biospecimen_ID in list(cnv_df.Kids_First_Biospecimen_ID):
            cnv_df.loc[cnv_df.Kids_First_Biospecimen_ID == row.Kids_First_Biospecimen_ID, "clinical_status"] = row.molecular_subtype
            cnv_df.loc[cnv_df.Kids_First_Biospecimen_ID == row.Kids_First_Biospecimen_ID, "cohort"] = row.cohort
            cnv_df.loc[cnv_df.Kids_First_Biospecimen_ID == row.Kids_First_Biospecimen_ID, "experimental_strategy"] = row.experimental_strategy
    
    # Filter dataframe to exclude samples with non-MYCN molecular subtype values in the `clinical_status` columns. These are mostly PBTA cohort samples and samples with `Unknown` molecular subtypes
    cnv_df = cnv_df[(cnv_df.clinical_status == "MYCN amp") | (cnv_df.clinical_status == "MYCN non-amp")] 
    return(cnv_df)  


def compute_classification_metrics(cnv_df):
    """
    Compute confusion matrix to evaluate the accuracy of CNV calls for MYCN amplification status. CNV MYCN status is classified relative to the clinical MYCN status to determine a suitable copy number cutoff for calling a MYCN gene in a sample as either amplified or not amplified.
    Parameter: a Pandas dataframe with Neuroblastoma samples and their associated MYCN gene status 
    Returns: a Pandas dataframe with classification metrics, including accuracy, sensitivity (TPR), specificity (TNR), and ROC AUC score to determine a suitable cutoff for copy number MYCN amplification
    """
    metrics_dfs = []
    copy_number_cutoff = [x for x in range(int(cnv_df.copy_number.min()) + 1, 16)]
    cnv_df = cnv_df.replace({"MYCN amp": 1, "MYCN non-amp": 0})
    for cutoff in copy_number_cutoff:
        cnv_df = cnv_df[cnv_df.copy_number.astype("int64") >= cutoff]
        CM = confusion_matrix(list(cnv_df.clinical_status), list(cnv_df.cnv_status))
        TN = CM[0][0]
        FN = CM[1][0]
        TP = CM[1][1]
        FP = CM[0][1]
        if (TN + FP) == 0:
            TNR = 1.0
        else:
            TNR = round(TN / (TN + FP), 6)
        TPR = round(TP / (TP + FN), 6)
        Accuracy = round((TP + TN)/(TN + FN + TP + FP), 6)
        TNRTPR = round((TNR + TPR), 6)
        ROCAUC = roc_auc_score((list(cnv_df.clinical_status), list(cnv_df.cnv_status)), 6)
        data = OrderedDict()
        data["Copy_number"] = cutoff
        data["TNs"] = TN
        data["FNs"] = FN
        data["TPs"] = TP
        data["TNR"] = TNR
        data["TPR"] = TPR
        data["Accuracy"] = Accuracy
        data["TNR + TPR"] = TNRTPR
        data["ROCAUC"] = ROCAUC
        df = pd.DataFrame(data, index=[copy_number_cutoff.index(cutoff)])
        metrics_dfs.append(df)
    metrics_df = pd.concat(metrics_dfs)
    return(metrics_df)


def get_mycn_genomic_segments(variant_file, chr_col, seg_start_col, seg_end_col, sample_col, cnv_df):
    """
    Get MYCN genomic segments for MYCN samples with amplification status that does not match between the CNV calls and the clinical data for CNVkit, Control-FREEC, and mantaSV callers
    Parameter: a Pandas dataframe with Neuroblastoma samples and their associated MYCN gene status 
    Returns: a Pandas dataframe with genomic segments for MYCN samples and their corresponding segment mean value
    """
    # Load CNV segments file into a pandas dataframe and filter to retain only sample rows with the MYCN region listed in consensus CNV file (2p24.3 regions with gene ENSG00000134323)
    seg_df = pd.read_csv(variant_file, sep="\t", dtype=str)
    seg_df["build"] = "GRCh37/hg19"
    # From the UCSC Genome Brower, these are the coordinates corresponding to the region 2p24.3, a GRCh37/hg19 multi-region: chr2:12000,001-16700000, includinhg MYCN gene ENSG00000134323.11 (chr2:16080686-16087129)
    seg_start = 12200001
    seg_end = 16700000
    mycn_seg_df = seg_df[((seg_df[chr_col] == "2") | (seg_df[chr_col] == "chr2")) & (seg_df[seg_start_col].astype("float").astype("int64") >= seg_start) & (seg_df[seg_end_col].astype("float").astype("int64") <= seg_end)]

    # Filter for only the relevant MYCN sample IDs with amplification status that do not match between the CNV calls and the clinical
    cnv_diff_df = cnv_df[cnv_df.cnv_status != cnv_df.clinical_status]
    caller_mycn_seg_df = mycn_seg_df[mycn_seg_df[sample_col].isin(cnv_diff_df.Kids_First_Biospecimen_ID)]
    return(caller_mycn_seg_df)


def compute_caller_overlaps(cnvkit_mycn_df, controlfreec_mycn_df, mantaSV_mycn_df, cnv_diff_df):
    """
    Compute coverage overlaps between genomic segments for MYCN samples with amplification status that does not match between the CNV calls and the clinical data for CNVkit, Control-FREEC, and mantaSV callers
    Parameter: CNV caller Pandas dataframes for discordant segments between the consesus calls and clinical data 
    Returns: a Pandas dataframe with overlap coverage ratio values between callers within the MYCN cyogenetic region for samples where there is discordance between the consesus calls and clinical data
    """
    # create caller1 and caller2 bed files for MYCN segments
    cnvkit_bed_file = "results/cnvkit.bed"
    cnvkit_bed_df = cnvkit_mycn_df[["ID", "loc.start", "loc.end",]]
    cnvkit_bed_df.to_csv(cnvkit_bed_file, sep="\t", index=False, header=None, encoding="utf-8")
    controlfreec_bed_file = "results/controlfreec.bed"
    controlfreec_bed_df = controlfreec_mycn_df[["Kids_First_Biospecimen_ID", "start", "end",]]
    controlfreec_bed_df.to_csv(controlfreec_bed_file, sep="\t", index=False, header=None, encoding="utf-8")
    mantaSV_bed_file = "results/mantaSV.bed"
    mantaSV_bed_df = mantaSV_mycn_df[["Kids.First.Biospecimen.ID.Tumor", "SV.start", "SV.end",]]
    mantaSV_bed_df.to_csv(mantaSV_bed_file, sep="\t", index=False, header=None, encoding="utf-8")

    # summarize bedtools coverage between callers
    cnvkit_controlfreec_df = summarize_bedtools_caller_overlaps(cnvkit_bed_file, controlfreec_bed_file)
    cnvkit_controlfreec_df.rename(columns={"coverage": "cnvkit-controlfreec"}, inplace=True)
    cnvkit_controlfreec_df = cnvkit_controlfreec_df[["Kids_First_Biospecimen_ID", "cnvkit-controlfreec"]]
    cnvkit_mantaSV_df = summarize_bedtools_caller_overlaps(cnvkit_bed_file, mantaSV_bed_file)
    cnvkit_mantaSV_df.rename(columns={"coverage": "cnvkit-mantaSV"}, inplace=True)
    cnvkit_mantaSV_df = cnvkit_mantaSV_df[["Kids_First_Biospecimen_ID", "cnvkit-mantaSV"]]
    cnv_diffs_overlaps_df = pd.merge(cnvkit_controlfreec_df, cnvkit_mantaSV_df, how="outer", on=["Kids_First_Biospecimen_ID"])
    controlfreec_mantaSV_df = summarize_bedtools_caller_overlaps(controlfreec_bed_file, mantaSV_bed_file)
    controlfreec_mantaSV_df.rename(columns={"coverage": "controlfreec-mantaSV"}, inplace=True)
    controlfreec_mantaSV_df = controlfreec_mantaSV_df[["Kids_First_Biospecimen_ID", "controlfreec-mantaSV"]]
    cnv_diffs_overlaps_df = pd.merge(cnv_diffs_overlaps_df, controlfreec_mantaSV_df, how="outer", on=["Kids_First_Biospecimen_ID"])
    cnv_diffs_overlaps_df = pd.merge(cnv_diff_df, cnv_diffs_overlaps_df, how="left", on=["Kids_First_Biospecimen_ID"])
    cnv_diffs_overlaps_df.fillna(0.0, inplace=True)
    cnv_diffs_overlaps_df[["cnvkit-controlfreec", "cnvkit-mantaSV", "controlfreec-mantaSV"]] = cnv_diffs_overlaps_df[["cnvkit-controlfreec", "cnvkit-mantaSV", "controlfreec-mantaSV"]].apply(lambda x: round(x*100, 3))
    return(cnv_diffs_overlaps_df)


def summarize_bedtools_caller_overlaps(caller1_bed_file, caller2_bed_file):
    """
    Summarize base coverage overlaps in the MYCN cyogenetic region bewteen callers by selecting the larget bidirectional for the mean coverage of segments within the MYCN cytogenic region. 
    Parameter: Caller bed files for discordant segments between the consesus calls and clinical data
    Returns: a Pandas dataframe of with summarized bedtools caller coverage overlaps
    """
    # rename required bedtools output columns for determining caller overlap coverage
    col_list = [0,6]
    caller_1_list = ["Kids_First_Biospecimen_ID", "caller1"]
    caller_2_list = ["Kids_First_Biospecimen_ID", "caller2"]

    # summarize overlaps between CNVkit and Control-FREEC
    with open("results/caller1_coverage.bed", "w") as bed:
        subprocess.run(["bedtools", "coverage", "-a", caller1_bed_file, "-b", caller2_bed_file], stdout=bed, check=True)
    caller1_df = pd.read_csv("results/caller1_coverage.bed", sep="\t", usecols=col_list, names=caller_1_list)
    caller1_df = caller1_df.groupby("Kids_First_Biospecimen_ID").mean().reset_index()
    with open("results/caller2_coverage.bed", "w") as bed:
        subprocess.run(["bedtools", "coverage", "-a", caller2_bed_file, "-b", caller1_bed_file], stdout=bed, check=True)
    caller2_df = pd.read_csv("results/caller2_coverage.bed", sep="\t", usecols=col_list, names=caller_2_list)
    caller2_df = caller2_df.groupby("Kids_First_Biospecimen_ID").mean().reset_index()
    summarize_df = pd.merge(caller1_df, caller2_df, how="outer", on=["Kids_First_Biospecimen_ID"])
    summarize_df["coverage"] = summarize_df[["caller1", "caller2"]].max(axis=1)
    return(summarize_df)


def main():
# get input parameters
    args = read_parameters()

    # call functions
    cnv_df = get_mycn_data(args.CONSENSUS_CNV_FILE, args.HISTOLOGY_FILE)
    metrics_df = compute_classification_metrics(cnv_df)
    cnvkit_mycn_df = get_mycn_genomic_segments(args.CNVKIT_CNV_FILE, "chrom", "loc.start", "loc.end", "ID", cnv_df)
    controlfreec_mycn_df = get_mycn_genomic_segments(args.CONTROLFREEC_CNV_FILE, "chr", "start", "end", "Kids_First_Biospecimen_ID", cnv_df)
    mantaSV_mycn_df = get_mycn_genomic_segments(args.MANTASV_SV_FILE, "SV.chrom", "SV.start", "SV.end", "Kids.First.Biospecimen.ID.Tumor", cnv_df)
    cnv_diffs_df = cnv_df[cnv_df.cnv_status != cnv_df.clinical_status]
    cnv_diffs_overlaps_df = compute_caller_overlaps(cnvkit_mycn_df, controlfreec_mycn_df, mantaSV_mycn_df, cnv_diffs_df)

    # write summary files
    cnv_df.to_csv("results/mycn_cnv_vs_clinical_status.tsv", sep="\t", index=False, encoding="utf-8")
    metrics_df.to_csv("results/mycn_cnv_vs_clinical_status_classification_metrics.tsv", sep="\t", index=False, encoding="utf-8")
    cnvkit_mycn_df.to_csv("results/cnvkit_mycn_cnv_vs_clinical_diff_status.tsv", sep="\t", index=False, encoding="utf-8")
    controlfreec_mycn_df.to_csv("results/controlfreec_mycn_cnv_vs_clinical_diff_status.tsv", sep="\t", index=False, encoding="utf-8")
    mantaSV_mycn_df.to_csv("results/mantaSV_mycn_cnv_vs_clinical_diff_status.tsv", sep="\t", index=False, encoding="utf-8")
    cnv_diffs_overlaps_df.to_csv("results/mycn_cnv_vs_clinical_diff_status.tsv", sep="\t", index=False, encoding="utf-8")

    # remove temporary bed files
    results_files = os.listdir("results")
    bed_files = [file for file in results_files if file.endswith(".bed")]
    for file in bed_files:
        bed_file = os.path.join("results", file)
        os.remove(bed_file)
    sys.exit(0)


if __name__ == '__main__':
     main()