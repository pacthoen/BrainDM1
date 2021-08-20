###############################
### prepareSampleMetadata.R ###
###############################

# this script loads, cleans and summarizes the required sample metadata for the downstream analysis

# load libraries
suppressPackageStartupMessages( library(dplyr) )
suppressPackageStartupMessages( library(data.table) )
suppressPackageStartupMessages( library(stringr) )
suppressPackageStartupMessages( library(purrr) )

# directories
repo_dir <- "/home/maxd/brainDM1"
data_dir <- "/mnt/xomics/maxd/data"
output_dir <- file.path(repo_dir, "test")

# load overall sample selection
sample_selection <- data.table::fread(file.path(repo_dir, "lib/sample_selection.txt"))


##################################
### prepare metadata for BrainSpan

# creates sample annotation table based on decrypted BrainSpan sample and subject metadata (retrieved from dbGaP)

dataset <- "BrainSpan"
metadata_dir <- file.path(data_dir, dataset, "metadata")

# required metadata files (must be retrieved from dbGaP and SRA)
sra_run_table <- "SraRunTable.txt"
sample_metadata <- "phs000755.v2.pht004992.v2.p1.BrainSpan_Atlas_Project_Sample.MULTI.txt"
subject_metadata <- "phs000755.v2.pht004993.v1.p1.c1.BrainSpan_Atlas_Project_Subject_Phenotypes.GRU.txt"

# load files
sra_run_table <- data.table::fread(file.path(metadata_dir, sra_run_table))
sample_metadata <- data.table::fread(file.path(metadata_dir, sample_metadata), skip = 10)
subject_metadata <- data.table::fread(file.path(metadata_dir, subject_metadata), skip = 10)

# merge individual tables
BrainSpan_metadata <- merge(
  x = sample_metadata,
  y = sra_run_table, 
  by.x = "SAMPID",
  by.y = "biospecimen_repository_sample_id") %>%
  merge(
    y = subject_metadata,
    by = c("dbGaP_Subject_ID", "SUBJID"))

rm(sra_run_table, sample_metadata, subject_metadata)

# get vector of BrainSpan sampleIDs
BrainSpan_samples <- sample_selection$id[sample_selection$dataset == "BrainSpan"]

# filter and clean up table
BrainSpan_metadata  <- BrainSpan_metadata %>%
  dplyr::filter(SAMPID %in% BrainSpan_samples) %>%
  dplyr::select(SUBJID, SAMPID, Run, body_site, Age, Sex) %>%
  dplyr::rename(subjectID = SUBJID, 
                sampleID = SAMPID, 
                sraAccession = Run, 
                tissue = body_site,
                age_in_days = Age,
                sex = Sex) %>%
  dplyr::mutate(dataset = "BrainSpan",
                sex = dplyr::recode(sex, "M" = "Male", "F" = "Female"),
                group = factor(ifelse(age_in_days < 280, "Prenatal", "Postnatal"), levels = c("Prenatal", "Postnatal")))


############################################
### prepare metadata for Otero et al. (2021)

# creates sample annotation table based on Otero et al. (2020) sample and subject metadata

dataset <- "Otero2020"
metadata_dir <- file.path(data_dir, dataset, "metadata")

# required metadata files (must be retrieved from GEO and SRA)
sra_run_table <- "SraRunTable.txt"
supplemental_TS1 <- "GSE157428_Supplemental_TS1.txt"

# load files
sra_run_table <- data.table::fread(file.path(metadata_dir, sra_run_table))
supplemental_TS1 <- data.table::fread(file.path(metadata_dir, supplemental_TS1))

# match GEO accession (exp) to patient id (e.g. "DM1_1")
# note: matches have been retrieved from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157428
GEOaccession_patient_match <- data.frame(GEOaccession = c(paste0("GSM47647", seq(39,59)),
                                                          paste0("GSM47647", seq(64,71))),
                                         Patient = c(paste0("DM1_", seq(1,21)),
                                                     paste0("Unaff_", seq(1,8))))

# merge individual tabes
Otero_metadata <- merge(
  x = sra_run_table, 
  y = GEOaccession_patient_match, 
  by.x = "GEO_Accession (exp)",
  by.y = "GEOaccession") %>%
  merge(
    y = supplemental_TS1,
    by = "Patient"
  )
rm(sra_run_table, supplemental_TS1, GEOaccession_patient_match)

# get vector of Otero sampleIDs
Otero_samples <- sample_selection$id[sample_selection$dataset == "Otero et al. (2021)"]

# filter and clean up table
Otero_metadata <- Otero_metadata %>%
  dplyr::filter(Run %in% Otero_samples) %>%
  dplyr::select(Run, Patient, TISSUE, Age, Sex, Status) %>%
  dplyr::rename(sampleID = Run,
                subjectID = Patient,
                tissue = TISSUE,
                age = Age,
                sex = Sex,
                group = Status) %>%
  dplyr::mutate(sraAccession = sampleID,
                dataset = "Otero et al. (2021)",
                group = dplyr::recode(group, "Unaff" = "Unaffected"))


#############################
### prepare metadata for GTEx

dataset <- "GTEx_v8"
metadata_dir <- file.path(data_dir, dataset, "metadata")

# required metadata files (must be retrieved from the GTEx AnVIL workspace after acquiring dbGaP access)
sample_metadata <- "GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt"
subject_metadata <- "GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt"

# load files
sample_metadata <- data.table::fread(file.path(metadata_dir, sample_metadata))
subject_metadata <- data.table::fread(file.path(metadata_dir, subject_metadata))

# extract subjectID from sampleID
sample_metadata$SUBJID <- substring(sample_metadata$SAMPID, 0, 10)

# merge individual tables
GTEx_metadata <- merge(
  x = sample_metadata,
  y = subject_metadata,
  by = "SUBJID"
)
rm(sample_metadata, subject_metadata)

# get vector of GTEx sampleIDs
GTEx_samples <- sample_selection$id[sample_selection$dataset == "GTEx"]

# filter and clean up table 
GTEx_metadata <- GTEx_metadata %>%
  dplyr::filter(SAMPID %in% GTEx_samples) %>%
  dplyr::select(SAMPID, SUBJID, SMTSD, AGE, SEX) %>%
  dplyr::rename(sampleID = SAMPID,
                subjectID = SUBJID,
                tissue = SMTSD,
                age = AGE,
                sex = SEX) %>%
  dplyr::mutate(tissue = str_replace(tissue, "Brain - ", ""),
                sex = dplyr::recode(sex, "1" = "Male", "2" = "Female"),
                dataset = "GTEx",
                group = "GTEx")


######################
### summarize metadata

sample_metadata <- data.table::rbindlist(list(BrainSpan_metadata, Otero_metadata, GTEx_metadata), 
                                        idcol = FALSE, use.names = TRUE, fill = TRUE)

sample_metadata <- sample_metadata[, c("sampleID", "dataset", "group", 
                                     "tissue", "sex", "age", "age_in_days", 
                                     "subjectID", "sraAccession")]

data.table::fwrite(sample_metadata, file.path(repo_dir, "data/sample_metadata.csv"))

