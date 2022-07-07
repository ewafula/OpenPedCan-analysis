# load libraries
library(magrittr)
library(dplyr)
library(readr)

# base directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "independent-samples")
out_dir <- file.path(analysis_dir, "results")
dir.create(out_dir, showWarnings = F, recursive = T)

# source functions
source(file.path(analysis_dir, "util", "independent-dna-samples.R"))
source(file.path(analysis_dir, "util", "independent-rna-samples.R"))


# read histology file
histology_df <- readr::read_tsv(file.path(root_dir, 'data/histologies-base.tsv'), guess_max=100000)

# randomize rows of histology file to avoid selection bias
set.seed(100)
histology_df <- histology_df[sample(nrow(histology_df)), ]

# Filter for dna tumors samples, where composition is not "Derived Cell Line".
dna_samples <- histology_df %>%
  dplyr::filter(sample_type == "Tumor", 
                composition != "Derived Cell Line", 
                is.na(RNA_library), 
                experimental_strategy %in% c("WGS", "WXS", "Targeted Sequencing"))

# filter for dna wgs samples
wgs_primary_plus_all <- dna_samples %>% 
  dplyr::filter(experimental_strategy == "WGS") %>% 
  independent_dna_samples(tumor_types = "prefer_primary", 
                          independent_level = "all-cohorts-pre-release", seed = 2020)

# filter for dna wxs samples
wxs_primary_plus_all <- dna_samples %>% 
  dplyr::filter(experimental_strategy == "WXS") %>% 
  independent_dna_samples(tumor_types = "prefer_primary", 
                          independent_level = "all-cohorts-pre-release", seed = 2020)

# filter for dna panel samples
panel_primary_plus_all <- dna_samples %>% 
  dplyr::filter(experimental_strategy == "Targeted Sequencing") %>% 
  independent_dna_samples(tumor_types = "prefer_primary", 
                          independent_level = "all-cohorts-pre-release", seed = 2020)
# dna independent samples 
independent_dna_sample_df_all <- rbind(wgs_primary_plus_all, wxs_primary_plus_all, panel_primary_plus_all) %>%
  distinct(Kids_First_Participant_ID, .keep_all = TRUE) %>% 
  dplyr::arrange(Kids_First_Biospecimen_ID)

# Filter for tumor samples where composition is not Derived Cell Line
rnaseq_samples <- histology_df %>%
  dplyr::filter(sample_type == "Tumor", composition != "Derived Cell Line")

# write independent sample outputs for independent levels of all cohorts 
rnaseq_primary_all_file <- file.path(out_dir, "independent-specimens.rnaseqpanel.primary.pre-release.tsv")
independent_rna_primary_all <- rnaseq_samples %>%
  independent_rna_samples(independent_dna_sample_df = independent_dna_sample_df_all,
                          independent_level = "all-cohorts-pre-release",
                          histology_df = .,
                          match_type = "independent_dna_plus_only_rna",
                          tumor_description_rna_only = "primary",
                          seed = 2020) %>%
  dplyr::arrange(Kids_First_Biospecimen_ID) %>% 
  readr::write_tsv(rnaseq_primary_all_file)

rnaseq_relapse_all_file <- file.path(out_dir, "independent-specimens.rnaseqpanel.relapse.pre-release.tsv")
independent_rna_relapse_all <- rnaseq_samples %>%
  independent_rna_samples(independent_dna_sample_df = independent_dna_sample_df_all,
                          independent_level = "all-cohorts-pre-release",
                          histology_df = .,
                          match_type = "independent_dna_plus_only_rna",
                          tumor_description_rna_only = "relapse",
                          seed = 2020) %>%
  dplyr::arrange(Kids_First_Biospecimen_ID) %>% 
  readr::write_tsv(rnaseq_relapse_all_file)

rnaseq_primplus_all_file <- file.path(out_dir, "independent-specimens.rnaseqpanel.primary-plus.pre-release.tsv")
independent_rna_primary_plus_all <- rnaseq_samples %>%
  independent_rna_samples(independent_dna_sample_df = independent_dna_sample_df_all,
                          independent_level = "all-cohorts-pre-release",
                          histology_df = .,
                          match_type = "independent_dna_plus_only_rna",
                          tumor_description_rna_only = "primary_plus",
                          seed = 2020) %>%
  dplyr::arrange(Kids_First_Biospecimen_ID) %>% 
  readr::write_tsv(rnaseq_primplus_all_file)
