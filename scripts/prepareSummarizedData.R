###############################
### prepareSummarizedData.R ###
###############################

# This script is used to create summarized data files for the Otero et al. (2021), BrainSpan and GTEx studies
# These summarized files can be used to reproduce the analysis of the dm1 splicing correlation paper

# load libraries
suppressPackageStartupMessages( library(dplyr) )
suppressPackageStartupMessages( library(data.table) )
suppressPackageStartupMessages( library(stringr) )
suppressPackageStartupMessages( library(rtracklayer) )
suppressPackageStartupMessages( library(edgeR) )
suppressPackageStartupMessages( library(purrr) )

# directories
data_dir <- "/mnt/xomics/maxd/data"

# load sample metadata
sample_metadata <- data.table::fread("data/sample_metadata.csv")

# load annotation files
gene_annotation <- as.data.table(import.gff("lib/gencode.v26.annotation.collapsed.gtf"))
miso_annotation <- import.gff("lib/SE.hg38.annotated.gff3")

# filter annotation for skipped exon events
miso_annotation <- miso_annotation[which(miso_annotation$type == "exon" & stringr::str_detect(miso_annotation$ID, "se"))]

# utilitly function to extract column from data files
extractCol <- function(file, ids, match_col, extract_col, pattern, sample_metadata){
  data <- fread(file)
  data <- data.table(data[match(ids, data[[match_col]]), ..extract_col,])
  
  sample <- stringr::str_extract(basename(file), pattern)
  
  if (sample %in% sample_metadata$sampleID) {
    colnames(data)[ncol(data)] <- sample
  } else if (sample %in% sample_metadata$sraAccession) {
    colnames(data)[ncol(data)] <- sample_metadata$sampleID[match(sample, sample_metadata$sraAccession)]
  } else {
    stop("Could not find valid column name for sample: ", sample)
  }
  
  data
}

# utility function to normalize raw counts
normalizeCounts <- function(raw_counts, ids, group) {
  
  # create DGE object with gene count data
  dge <- edgeR::DGEList(counts = raw_counts, genes = ids, group = group)
  
  # filter genes with default parameters
  keep <- edgeR::filterByExpr(dge)
  dge <- dge[keep, ]
  
  # TMM normalization
  dge <- edgeR::calcNormFactors(dge)
  
  # calculate logCPM values
  log_CPM <- data.table(geneID = dge$genes$genes, cpm(dge, log = TRUE, prior.count = 3))
  
  log_CPM
}

####################
### prepare psi data

# get paths to all miso summary files
miso_summary_files <- c(
  list.files(file.path(data_dir, "BrainSpan/miso/summary"), full.names = T),
  list.files(file.path(data_dir, "Otero2020/miso/summary"), full.names = T),
  list.files(file.path(data_dir, "GTEx_v8/miso/summary"), full.names = T)
)

# drop GTEx frontal cortex samples since they will be added from another directory
drop <- str_detect(
  string = miso_summary_files, 
  pattern = paste(sample_metadata$sampleID[sample_metadata$dataset == "GTEx" & sample_metadata$tissue == "Frontal Cortex (BA9)"], collapse = "|")
)

miso_summary_files <- miso_summary_files[!drop]

miso_summary_files <- c(miso_summary_files,
  list.files(file.path(data_dir, "GTEx_v8/anvil/frontal_cortex/miso/summary"), full.names = T)
)

keep <- stringr::str_extract(basename(miso_summary_files), paste0("(.*)(?=\\.miso_summary)")) %in% 
  c(sample_metadata$sampleID, sample_metadata$sraAccession)

miso_summary_files <- miso_summary_files[keep]

psi <- data.table(eventID = stringr::str_replace(miso_annotation$ID, ".A.se", ""))

# note: collecting the psi data might take a while
psi_data <- cbind(psi, map_dfc(
  .x = miso_summary_files, 
  .f = extractCol, 
  ids = psi$eventID, 
  match_col = "event_name",
  extract_col = "miso_posterior_mean",
  pattern = "(.*)(?=\\.miso_summary)",
  sample_metadata = sample_metadata
))

fwrite(psi_data, "data/psi_data.csv")

#####################
### prepare gene data

gene_count_files <- c(
  list.files(file.path(data_dir, "BrainSpan/rnaseqc_output"), pattern = "gene_reads", full.names = T, recursive = T),
  list.files(file.path(data_dir, "Otero2020/rnaseqc_output"), pattern = "gene_reads", full.names = T, recursive = T)
)

keep <- stringr::str_extract(basename(gene_count_files), "(^SRR[0-9]+)(?=\\.Aligned)") %in% sample_metadata$sraAccession

gene_count_files <- gene_count_files[keep]

gene_counts <- data.table(geneID = fread(gene_count_files[1])$Name)

# note: collecting the gene count data might take a while
gene_counts <- cbind(gene_counts, map_dfc(
  .x = gene_count_files, 
  .f = extractCol, 
  ids = gene_counts$geneID, 
  match_col = "Name",
  extract_col = "Counts",
  pattern = "(^SRR[0-9]+)(?=\\.Aligned)",
  sample_metadata = sample_metadata
))

BrainSpan_counts <- subset(gene_counts, select = c("geneID", sample_metadata$sampleID[sample_metadata$dataset == "BrainSpan"]))
Otero_counts <- subset(gene_counts, select = c("geneID", sample_metadata$sampleID[sample_metadata$dataset == "Otero et al. (2021)"]))  

# load gene counts for GTEx (file was retrieved from GTEx portal)
GTEx_counts <- fread(file.path(data_dir, "GTEx_v8/public/gene_level/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct"))
GTEx_counts <- subset(GTEx_counts, select = c("Name", sample_metadata$sampleID[sample_metadata$dataset == "GTEx"]))
GTEx_counts <- dplyr::rename(GTEx_counts, geneID = Name)

# normalize and filter counts
BrainSpan_log_cpm <- normalizeCounts(
  raw_counts = as.matrix(BrainSpan_counts[, -"geneID"]), 
  ids = BrainSpan_counts$geneID, 
  group = sample_metadata$tissue[match(colnames(BrainSpan_counts)[-1], sample_metadata$sampleID)]
)  
Otero_log_cpm <- normalizeCounts(
  raw_counts = as.matrix(Otero_counts[, -"geneID"]), 
  ids = Otero_counts$geneID, 
  group = sample_metadata$tissue[match(colnames(Otero_counts)[-1], sample_metadata$sampleID)]
)
GTEx_log_cpm <- normalizeCounts(
  raw_counts = as.matrix(GTEx_counts[, -"geneID"]), 
  ids = GTEx_counts$geneID, 
  group = sample_metadata$tissue[match(colnames(GTEx_counts)[-1], sample_metadata$sampleID)]
)

log_cpm_data <- merge(
  x = BrainSpan_log_cpm,
  y = Otero_log_cpm,
  by = "geneID", all = TRUE) %>% 
  merge(
    y = GTEx_log_cpm,
    by = "geneID", all = TRUE
  )

# load genes of interest
splicing_factors <- fread("lib/splicing_factors.txt")

# prepare filtered gene-level data (log2 CPM)
log_cpm_data_filtered <- log_cpm_data %>%
  dplyr::filter(geneID %in% splicing_factors$geneID) %>%
  tidyr::gather(key = "sampleID", value = "log2CPM", -geneID) %>%
  merge(splicing_factors, by = "geneID") %>%
  merge(sample_metadata, by = "sampleID")

# save data frame
fwrite(log_cpm_data_filtered, "data/log2CPM_data_filtered_annotated.csv")
