# Analysis code for "Alternative MAPK drivers in KRAS WT pancreatic cancer - Letter"
# Calculate cancer effect sizes in PAAD data with epistasis analysis
# Last updated: April 25, 2024

# load libraries
library(cancereffectsizeR) # using version 2.8.0
library(tidyverse)
library(ces.refset.hg19) # using version 1.2.2
library(TCGAbiolinks)
library(ggplot2)
library(data.table)
library(dplyr)

# define CESA object
paad_cesa <- cancereffectsizeR::CESAnalysis(refset = ces.refset.hg19)

# import whole genome/exome datasets
# files were manually downloaded from links provided

# generate TCGA data
# generate TCGA maf file, match records and find matching IDs
tcga_paad <- "source_data/TCGA/TCGA-PAAD.maf.gz"
tcga_clinical <- GDCquery_clinic(project = "TCGA-PAAD", type = "clinical")

if (!file.exists(tcga_paad)) {
  get_TCGA_project_MAF(project = "PAAD", filename = tcga_paad, 
                       exclude_TCGA_nonprimary = FALSE)
}

maf_tcga_PAAD <- preload_maf(maf = tcga_paad, refset = "ces.refset.hg19", chain_file = "source_data/TCGA/hg38ToHg19.over.chain", detect_hidden_mnv = FALSE)
maf_tcga_PAAD <- maf_tcga_PAAD[germline_variant_site == F][repetitive_region == F | cosmic_site_tier %in% 1:3]


# UTSW
# Witkiewicz, A. K., McMillan, E. A., Balaji, U., Baek, G., Lin, W. C., Mansour, J., et al.
# Whole-exome sequencing of pancreatic cancer defines genetic diversity and therapeutic targets. Nature communications, 6, 6744.
# https://cbioportal-datahub.s3.amazonaws.com/paad_utsw_2015.tar.gz 
maf_utsw_PAAD <- preload_maf(maf = "source_data/UTSW/data_mutations_extended.txt",
                             refset = ces.refset.hg19)
maf_utsw_PAAD = maf_utsw_PAAD[germline_variant_site == F][repetitive_region == F | cosmic_site_tier %in% 1:3]

# ICGC
# Biankin, A. V., Waddell, N., Kassahn, K. S., Gingras, M. C., Muthuswamy, L. B., Johns, A. L., et al.
# Pancreatic cancer genomes reveal aberrations in axon guidance pathway genes. Nature, 491(7424), 399–405.
# https://cbioportal-datahub.s3.amazonaws.com/paad_icgc.tar.gz 
maf_icgc_PAAD <- preload_maf(maf = "source_data/ICGC/data_mutations.txt",
                             refset = ces.refset.hg19)
maf_icgc_PAAD = maf_icgc_PAAD[germline_variant_site == F][repetitive_region == F | cosmic_site_tier %in% 1:3]

# QCMG
# Bailey, P., Chang, D. K., Nones, K., Johns, A. L., Patch, A. M., Gingras, M. C., et al.
# Genomic analyses identify molecular subtypes of pancreatic cancer. Nature, 531(7592), 47–52.
# https://cbioportal-datahub.s3.amazonaws.com/paad_qcmg_uq_2016.tar.gz 
maf_qcmg_PAAD <- preload_maf(maf = "source_data/QCMG/data_mutations_extended.txt",
                             refset = ces.refset.hg19)
maf_qcmg_PAAD = maf_qcmg_PAAD[germline_variant_site == F][repetitive_region == F | cosmic_site_tier %in% 1:3]

# CPTAC 
# Cao, L., Huang, C., Cui Zhou, D., Hu, Y., Lih, T. M., Savage, S. R., et al.
# Proteogenomic characterization of pancreatic ductal adenocarcinoma. Cell, 184(19), 5031–5052.e26.
# https://cbioportal-datahub.s3.amazonaws.com/paad_cptac_2021.tar.gz
maf_paad_cptac_2021 <- preload_maf("source_data/CPTAC/data_mutations.txt",refset = ces.refset.hg19)
maf_paad_cptac_2021 = maf_paad_cptac_2021[germline_variant_site == F][repetitive_region == F | cosmic_site_tier %in% 1:3]


# Yale-Gilead collaboration
# https://pubmed.ncbi.nlm.nih.gov/30365005/
maf_yg <- preload_maf(maf = "source_data/Yale_Gilead/mutationsTN_26_Pancreatic_Cancer.maf",chr_col = "Chrom",sample_col = "Patient_ID",refset = ces.refset.hg19)
maf_yg = maf_yg[germline_variant_site == F][repetitive_region == F | cosmic_site_tier %in% 1:3]


# import targeted sequencing dataset

# GENIE
# The AACR Project GENIE Consortium. AACR Project GENIE: Powering Precision Medicine Through An International Consortium, Cancer Discov. 2017 Aug;7(8):818-831

genie_clin <- read_tsv("source_data/GENIE/data_clinical_sample.txt", comment="#")
genie_coverage <- read_tsv("source_data/GENIE/genomic_information.txt")
genie_maf <- read_tsv("source_data/GENIE/data_mutations_extended.txt")

genie_clin_paad <- genie_clin %>% 
  filter(CANCER_TYPE == "Pancreatic Cancer" &
           (CANCER_TYPE_DETAILED == "Pancreatic Adenocarcinoma" |
              CANCER_TYPE_DETAILED == "Adenosquamous Carcinoma of the Pancreas" |
              CANCER_TYPE_DETAILED == "Acinar Cell Carcinoma of the Pancreas"))

genie_maf <- genie_maf %>% 
  filter(Tumor_Sample_Barcode %in% genie_clin_paad$SAMPLE_ID)

# exclude germline mutations in the analysis
genie_maf <- genie_maf %>%
  filter(Chromosome != "X")

genie_maf %>% 
  mutate(Chromosome = case_when(
    Chromosome == "23" ~ "X" , 
    TRUE ~ Chromosome)) -> genie_maf

# adding in the patient information
genie_maf <- left_join(x = genie_maf, y=genie_clin_paad[,c("PATIENT_ID","SAMPLE_ID")],
                       by = c("Tumor_Sample_Barcode" = "SAMPLE_ID"))

# some patients have more than one sample. It is unclear if these samples are from the same tumor. 
# We will pick the sample with the largest number of detected substitutions. 
# find the samples with the most substitutions
genie_patients_to_keep <- genie_maf %>% 
  count(PATIENT_ID, Tumor_Sample_Barcode) %>% 
  group_by(PATIENT_ID) %>% 
  slice_max(order_by = n, n=1,with_ties = F) %>% 
  pull(Tumor_Sample_Barcode)

# filter for samples with the most substitutions 
genie_maf %>% 
  filter(Tumor_Sample_Barcode %in% genie_patients_to_keep) -> 
  genie_maf

genie_maf <- preload_maf(maf = genie_maf,refset = ces.refset.hg19,keep_extra_columns = c("PATIENT_ID"))

# filter for tumor subtypes 

# QCMG

QCMG_clin <- read_tsv(file = "source_data/QCMG/paad_qcmg_uq_2016_clinical_data.tsv")

QCMG_clin %>% 
  filter(`Tumor Other Histologic Subtype` %in% "Pancreatic Ductal Adenocarcinoma") %>% 
  mutate(same_patient_sample = `Patient ID` == `Sample ID`) %>%
  pull(`Patient ID`) -> 
  qcmg_patients_to_keep

maf_qcmg_PAAD %>% 
  filter(Unique_Patient_Identifier %in% qcmg_patients_to_keep) -> 
  maf_qcmg_PAAD

## TCGA

tcga_biospecimen <- read.csv("source_data/TCGA/PAAD.clin.merged.txt",sep = "\t",header = F)

tcga_biospecimen <- as.data.frame(t(tcga_biospecimen))
colnames(tcga_biospecimen) <- tcga_biospecimen[1,]

tcga_biospecimen <- tcga_biospecimen[-1,]


tcga_biospecimen %>% 
  filter(patient.histological_type == "pancreas-adenocarcinoma ductal type") %>% 
  pull(patient.bcr_patient_barcode) %>% 
  toupper() -> 
  tcga_samples_needed 

maf_tcga_PAAD %>% 
  mutate(bcr_patient_barcode = str_sub(Unique_Patient_Identifier,1,12)) %>% 
  filter(bcr_patient_barcode %in% tcga_samples_needed) -> 
  maf_tcga_PAAD

# check for duplications

unique(maf_icgc_PAAD$Unique_Patient_Identifier)
unique(maf_qcmg_PAAD$Unique_Patient_Identifier)

maf_icgc_PAAD <- maf_icgc_PAAD %>% 
  mutate(icgc_identifier = str_replace(string = Unique_Patient_Identifier,pattern = "_TD",replacement = ""))

length(which(unique(maf_qcmg_PAAD$Unique_Patient_Identifier) %in% unique(maf_icgc_PAAD$icgc_identifier)))

# filter out identical names 
maf_icgc_PAAD <- maf_icgc_PAAD %>%
  filter(!icgc_identifier %in% maf_qcmg_PAAD$Unique_Patient_Identifier)

# assemble maf data files

all_mafs <- list(maf1 = maf_icgc_PAAD,
                 maf2 =  maf_qcmg_PAAD,
                 maf3 = maf_utsw_PAAD, 
                 maf4 = maf_tcga_PAAD, 
                 maf5  = maf_paad_cptac_2021)


combined_maf <- data.table::rbindlist(all_mafs, idcol="source",fill=T)

combined_maf <- combined_maf[germline_variant_site == F & (repetitive_region == F | cosmic_site_tier %in% 1:3)]

possible_dups <- check_sample_overlap(combined_maf)

possible_dups %>%
  filter(variants_shared > 5)

to_remove <- c("ICGC_0229",
               "PCSI0022_T",
               "ICGC_0150",
               "TCGA-US-A774-01A-21D-A32N-08",
               "PCSI0007_T",
               "ICGC_0097",
               "PCSI0024_T",
               "ICGC_0164")

# load in substitution data

# panel data 

## GENIE

genie_coverage <- read_tsv(file = "source_data/GENIE/genomic_information.txt",  col_types = list( Chromosome = col_character()))

genie_coverage <- genie_coverage %>% 
  filter(Chromosome!="Un_gl000228")

# loop to load in genie data with correct panel 
for(seq_assay_ind in 1:length(unique(genie_clin_paad$SEQ_ASSAY_ID))){
  
  # get sequencing assay data
  this_seq_assay <- unique(genie_clin_paad$SEQ_ASSAY_ID)[seq_assay_ind]
  
  # get tumor IDs that use this assay 
  these_tumor_names <- genie_clin_paad %>% 
    filter(SEQ_ASSAY_ID == this_seq_assay) %>% 
    pull(SAMPLE_ID) %>% 
    unique()
  
  # get mutation data of these tumors
  these_tumors_maf <- genie_maf %>% 
    filter(Unique_Patient_Identifier %in% these_tumor_names)
  
  # load in each set up variants with the correct panel 
  if(nrow(these_tumors_maf) > 0 & !("GENIE-YALE-TPL478-1" %in% these_tumors_maf$Unique_Patient_Identifier)){
    
    this_coverage <- genie_coverage %>% 
      filter(SEQ_ASSAY_ID == this_seq_assay)
    
    this_coverage_granges <- GRanges(seqnames = this_coverage$Chromosome, ranges = IRanges(start = this_coverage$Start_Position,end = this_coverage$End_Position))
    
    # load into cesa
    
    paad_cesa <- load_maf(cesa = paad_cesa, maf = these_tumors_maf,
                          coverage = "targeted",
                          covered_regions = this_coverage_granges,
                          covered_regions_name = this_seq_assay)
  }
}

# load whole exome data
paad_cesa <- load_maf(cesa = paad_cesa, 
                      maf = maf_icgc_PAAD %>% 
                        filter(!Unique_Patient_Identifier %in% to_remove))

paad_cesa <- load_maf(cesa = paad_cesa, maf = maf_qcmg_PAAD %>% 
                        filter(! Unique_Patient_Identifier %in% to_remove))

paad_cesa <- load_maf(cesa = paad_cesa, maf = maf_utsw_PAAD %>% 
                        filter(! Unique_Patient_Identifier %in% to_remove))

paad_cesa <- load_maf(cesa = paad_cesa, maf = maf_paad_cptac_2021 %>% 
                        filter(! Unique_Patient_Identifier %in% to_remove))

paad_cesa <- load_maf(cesa = paad_cesa, maf = maf_tcga_PAAD %>% 
                        filter(!Unique_Patient_Identifier %in% to_remove))

# mutation rate estimates

## trinucleotide rates
paad_cesa <- cancereffectsizeR::trinuc_mutation_rates(cesa = paad_cesa,signature_set = "COSMIC_v3.2",
                                                      signatures_to_remove = cancereffectsizeR::suggest_cosmic_signatures_to_remove(cancer_type = "PAAD"),
                                                      cores = 2)

## gene by gene mutation rates
paad_cesa <- cancereffectsizeR::gene_mutation_rates(cesa = paad_cesa, 
                                                    covariates = "pancreas")

# calculate scaled selection coefficients for gene variants
paad_cesa <- ces_variant(cesa = paad_cesa, run_name = "kras")

# save cesa file for visualization in "paad_plotting.R"
cancereffectsizeR::save_cesa(cesa = paad_cesa, file = "paad_cesa.rds")