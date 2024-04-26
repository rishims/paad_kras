# Analysis of alternative MAPK drivers in KRAS WT pancreatic cancer
This repository contains the code and results related to a letter to the editor prepared in response to [Singh et al.'s](https://aacrjournals.org/clincancerres/article/29/22/4627/729969/Oncogenic-Drivers-and-Therapeutic-Vulnerabilities) original article regarding alternative MAPK drivers in KRAS WT pancreatic cancer. Here, we rank the relative contribution of somatic mutations on cancer cell lineage survival and proliferation in pancreatic adenocarcinomas (PAAD). We compiled an extensive sequencing dataset consisting of over 6,000 PAAD samples to perform our analysis. We also explore the pairwise epistatic effects between PAAD driver genes and how this influences selection for KRAS. The ultimate goal of such analyses is to better guide the development of targeted therapeutics for various cancers and their subtypes.

## Guide to run the analysis

1. Open and run the `paad_analysis.R` script. Please note that you must aquire the GENIE data first and details on how to do so can be found [here](https://www.aacr.org/professionals/research/aacr-project-genie/aacr-project-genie-data/).
2. Once the analysis is run, load the saved `paad_cesa.rds` into `paad_plotting.R` and run this script to produce the visualizations. 
