# Author: Komal S. Rathi
# Date: 11/26/2019
# Function: 
# 1. summarize RNA-seq to HUGO symbol x Sample matrix
# 2. tabulate corresponding gene annotations

# Example run: PolyA RNA-seq
# Rscript analyses/collapse-rnaseq/01-summarize_matrices.R \
# -i data/pbta-gene-expression-rsem-fpkm.polya.rds \
# -g data/gencode.v27.primary_assembly.annotation.gtf.gz \
# -m analyses/collapse-rnaseq/pbta-gene-expression-rsem-fpkm-collapsed.polya.rds \
# -t analyses/collapse-rnaseq/pbta-gene-expression-rsem-fpkm-collapsed_table.polya.rds

# Example run: Stranded RNA-seq
# Rscript analyses/collapse-rnaseq/01-summarize_matrices.R \
# -i data/pbta-gene-expression-rsem-fpkm.stranded.rds \
# -g data/gencode.v27.primary_assembly.annotation.gtf.gz \
# -m analyses/collapse-rnaseq/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds \
# -t analyses/collapse-rnaseq/pbta-gene-expression-rsem-fpkm-collapsed_table.stranded.rds

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(rtracklayer))

option_list <- list(
  make_option(c("-i", "--inputmat"), type = "character",
              help = "Input matrix of merged RSEM files (.RDS)"),
  make_option(c("-g", "--inputgtf"), type  = "character",
              help = "Input gtf file (.gtf)"),
  make_option(c("-m", "--outputmat"), type = "character",
              help = "Output matrix (.RDS)"),
  make_option(c("-t", "--dupstable"), type = "character",
              help = "Output gene annotation table (.RDS)"),
  make_option(c("-n", "--gtex_file"), type = "character",default=NULL,
              help = "GTEx downloaded file")
)

# parse the parameters
opt <- parse_args(OptionParser(option_list = option_list))
input.mat <- opt$inputmat
input.gtf <- opt$inputgtf
output.mat <- opt$outputmat
dups.tab <- opt$dupstable
gtex_file <- opt$gtex_file

print("Generating input matrix and drops table...!")

# read gencode (v27 for the paper) to get gene id-gene symbo-gene type annotation
# PAR_Y entries are just duplications - so remove beforehand
gtf <- rtracklayer::import(input.gtf)
gtf <- gtf %>% 
  as.data.frame() %>%
  select(gene_id, gene_name, gene_type) %>%
  mutate(gene_id = str_replace_all(gene_id, "_PAR_Y", "")) %>%
  unique()

# identify gene symbols mapped to multiple ensembl identifiers
dups <- gtf %>% 
  select(gene_id, gene_name) %>%
  unique() %>% 
  group_by(gene_name) %>%
  mutate(gene_symbol.ct = n()) %>%
  filter(gene_symbol.ct > 1) %>%
  unique()
gtf$ensembl_id <- ifelse(gtf$gene_id %in% dups$gene_id, "Multiple", "Unique")



if (!is.null(gtex_file)){
  print("Read merged GTEx data")
  # first in comment line and second just has the dimensions printed
  expr <- read_tsv(gtex_file,skip = 2) %>%
    dplyr::rename("gene_symbol"="Description",
                  "gene_id"="Name") %>%
    mutate(gene_id = str_replace_all(gene_id, "_PAR_Y", ""))
  
  # collapse to matrix of HUGO symbols x Sample identifiers
  # take mean per row and use the max value for duplicated gene symbols
  expr.collapsed <- expr %>% 
    mutate(means = rowMeans(select(.,-gene_id, -gene_symbol))) %>% # take rowMeans
    arrange(desc(means)) %>% # arrange decreasing by means
    distinct(gene_symbol, .keep_all = TRUE) %>% # keep the ones with greatest mean value. If ties occur, keep the first occurencce
    select(-means) %>%
    unique() %>%
    remove_rownames() 
  
}else{
  print("Read merged RSEM data")
  expr <- readRDS(input.mat)
  
  # split gene id and symbol
  expr <- expr %>% 
    mutate(gene_id = str_replace(gene_id, "_PAR_Y_", "_"))  %>%
    separate(gene_id, c("gene_id", "gene_symbol"), sep = "\\_", extra = "merge") %>%
    unique() 
  
  # collapse to matrix of HUGO symbols x Sample identifiers
  # take mean per row and use the max value for duplicated gene symbols
  expr.collapsed <- expr %>% 
    mutate(means = rowMeans(select(.,-gene_id, -gene_symbol))) %>% # take rowMeans
    arrange(desc(means)) %>% # arrange decreasing by means
    distinct(gene_symbol, .keep_all = TRUE) %>% # keep the ones with greatest mean value. If ties occur, keep the first occurencce
    select(-means) %>%
    unique() %>%
    remove_rownames() 
}

# matrix of HUGO symbols x Sample identifiers
expr.input <- expr.collapsed %>% 
  column_to_rownames("gene_symbol") %>%
  select(-c(gene_id)) 
print(dim(expr.input)) # 53011 for stranded and 46400 for polyA

# save matrix
print("Saving collapsed matrix...")
saveRDS(object = expr.input, file = output.mat)
print("Done!!")


