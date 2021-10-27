#!/usr/bin/env python


"""
03-mycn-SNParray-thresholding.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Functions to assess Neuroblastoma MYCN gene focal and segment SNP array copy number variation (CNV) amplification for OPenPedCan TARGET cohorts thresholding analyses modules 
"""


__author__ = ('Eric Wafula (wafulae@chop.edu)')
__version__ = '1.0'
__date__ = '27 October 2021'


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
     p = argparse.ArgumentParser(description=("The 03-mycn-SNParray-thresholding.py script creates summary results to enable focal and segment SNP arrays copy number variation (CNV) amplification thresholding assessment for Neuroblastoma MYCN gene."), formatter_class=argparse.RawTextHelpFormatter)
     p.add_argument('FOCAL_SNPARRAY_FILE', type=str, default=None, help="TARGET MYCN gene focal SNP array file\n\n")
     p.add_argument('PAIRED_SEGMENT_SNPARRAY_FILE', type=str, default=None, help="TARGET MYCN gene paired segmnent array file\n\n")
     p.add_argument('SINGLE_SEGMENT_SNPARRAY_FILE', type=str, default=None, help="TARGET MYCN gene single segmnent array file\n\n")
     p.add_argument('HISTOLOGY_FILE', type=str, default=None, help="OPenPedCan histology file (histologies.tsv)\n\n")
     p.add_argument('-v', '--version', action='version', version="03-mycn-SNParray-thresholding.py version {} ({})".format(__version__, __date__), help="Print the current 03-mycn-SNParray-thresholding.py version and exit\n\n")
     return p.parse_args()


def get_mean_cn_ratios(mycn_focal_file, mycn_seg_paired_file, mycn_seg_single_file):
    """
    Get CN calls for MYCN genes for each Neuroblastoma sample in both the focal and segment SNP array files
    Parameters: TARGET MYCN gene ocal SNP and paired/single segmnent array files
    Returns: a dict of Pandas dataframes with Neuroblastoma samples and their associated MYCN mean Log R raitos (focal) and seg.cn ratios (segment) SNP array calls
    """

    # load the focal and segment SNP array files and compute the mean copy number ratio for each sample
    mycn_dfs = {}
    # focal array
    focal_df = pd.read_csv(mycn_focal_file, sep="\t", dtype={"Log R Ratio": "float64"})
    focal_df = focal_df[["Sample ID", "Log R Ratio"]]
    focal_df.rename(columns={"Sample ID": "sample_id", "Log R Ratio": "cn_ratio"}, inplace=True)
    focal_df = focal_df.groupby("sample_id").mean().reset_index()
    mycn_dfs["focal"] = focal_df
    # segment paired array
    seg_paired_df = pd.read_csv(mycn_seg_paired_file, sep="\t", dtype={"Seg.CN": "float64"})
    seg_paired_df = seg_paired_df[["Sample", "Seg.CN"]]
    seg_paired_df.rename(columns={"Sample": "sample_id", "Seg.CN": "cn_ratio"}, inplace=True)
    seg_paired_df = seg_paired_df.groupby("sample_id").mean().reset_index()
    mycn_dfs["paired_segment"] = seg_paired_df
    # segment paired array
    seg_single_df = pd.read_csv(mycn_seg_single_file, sep="\t", dtype={"Seg.CN": "float64"})
    seg_single_df = seg_single_df[["Sample", "Seg.CN"]]
    seg_single_df.rename(columns={"Sample": "sample_id", "Seg.CN": "cn_ratio"}, inplace=True)
    seg_single_df = seg_single_df.groupby("sample_id").mean().reset_index()
    mycn_dfs["single_segment"] = seg_single_df
    return(mycn_dfs) 

def set_mycn_status(mycn_df, cn_ratio_cutoff, histologies_file):
    """
    Set the MYCN genes amplifcation status for each Neuroblastoma sample for the copy number ration cutoff
    Parameters: a Pandas dataframe with Neuroblastoma samples and their associated MYCN mean Log R raitos (focal) and seg.cn ratios (segment) SNP array calls
    Parameter: OpendeCan clinical histologies file
    Returns: a Pandas dataframe with Neuroblastoma samples and their associated array and clinical MYCN gene status
    """

    # set mycn amplification status for the input copy number ratio threshold
    mycn_df["array_status"] = np.where(mycn_df["cn_ratio"] >= cn_ratio_cutoff, "MYCN amp", "MYCN non-amp")
    # load histlogies files and get intersect (inner join) with for TARGET tumor WXS samples
    histology_df = pd.read_csv(histologies_file, sep="\t", na_filter=False, dtype=str)
    histology_df = histology_df[(histology_df.experimental_strategy == "WXS") & (histology_df.cohort == "TARGET") & (histology_df.short_histology == "Neuroblastoma") & (histology_df.sample_type == "Tumor")]
    histology_df = histology_df[["Kids_First_Participant_ID", "molecular_subtype"]]
    histology_df = histology_df.rename(columns={"Kids_First_Participant_ID": "sample_id", "molecular_subtype": "clinical_status"})
    histology_df = histology_df[(histology_df.clinical_status == "MYCN amp") | (histology_df.clinical_status == "MYCN non-amp")]
    mycn_status_df = pd.merge(mycn_df, histology_df, on="sample_id", how="inner").reset_index(drop=True)
    return(mycn_status_df)


def compute_classification_metrics(mycn_status_df, cn_ratio_cutoff):
    """
    Compute confusion matrix to evaluate the accuracy of CNV calls for MYCN amplification status. CNV MYCN status is classified relative to the clinical MYCN status to determine a suitable copy number cutoff for calling a MYCN gene in a sample as either amplified or not amplified.
    Parameter: a Pandas dataframe with Neuroblastoma samples and their associated MYCN gene status 
    Parameter: copy number ratio cutoff for either SNP or segment array 
    Returns: a list of classification metrics for the copy number ratio cutoff, including accuracy, sensitivity (TPR), specificity (TNR), and ROC AUC score to determine a suitable cutoff for copy number MYCN amplification
    """
    # recode mycn amplification stutus for classification
    mycn_status_df = mycn_status_df.replace({"MYCN amp": 1, "MYCN non-amp": 0})
    # compute classification  metrics with clinical mycn status as ground truth
    CM = confusion_matrix(list(mycn_status_df.clinical_status), list(mycn_status_df.array_status))
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
    ROCAUC = round(roc_auc_score(list(mycn_status_df.clinical_status), list(mycn_status_df.array_status)), 6)
    mycn_metrics = OrderedDict()
    mycn_metrics["CN Ratio"] = cn_ratio_cutoff
    mycn_metrics["TNs"] = TN
    mycn_metrics["FNs"] = FN
    mycn_metrics["TPs"] = TP
    mycn_metrics["TNR"] = TNR
    mycn_metrics["TPR"] = TPR
    mycn_metrics["Accuracy"] = Accuracy
    mycn_metrics["TNR + TPR"] = TNRTPR
    mycn_metrics["ROCAUC"] = ROCAUC
    return(mycn_metrics)


def main():
# get input parameters
    args = read_parameters()

    # call functions
    mycn_dfs = get_mean_cn_ratios(args.FOCAL_SNPARRAY_FILE, args.PAIRED_SEGMENT_SNPARRAY_FILE, args.SINGLE_SEGMENT_SNPARRAY_FILE)
    for key, value in mycn_dfs.items():
        mycn_metrics_dfs = []
        cn_ratio_cutoff_list = [x/10 for x in range(1, 21)]
        for cn_ratio_cutoff in cn_ratio_cutoff_list:
            mycn_status_df = set_mycn_status(value, cn_ratio_cutoff, args.HISTOLOGY_FILE)
            status_file = "results/mycn_{}_SNParray_vs_clinical_status.tsv".format(key)
            mycn_status_df.to_csv(status_file, sep="\t", index=False, encoding="utf-8")
            mycn_metrics = compute_classification_metrics(mycn_status_df, cn_ratio_cutoff)
            df = pd.DataFrame(mycn_metrics, index=[cn_ratio_cutoff_list.index(cn_ratio_cutoff)])
            mycn_metrics_dfs.append(df)
            mycn_metrics_df = pd.concat(mycn_metrics_dfs)
            metrics_file = "results/mycn_{}_SNParray_vs_clinical_metrics.tsv".format(key)
            mycn_metrics_df.to_csv(metrics_file, sep="\t", index=False, encoding="utf-8")
    sys.exit(0)


if __name__ == '__main__':
     main()