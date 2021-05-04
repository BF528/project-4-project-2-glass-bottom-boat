# Project Description

Baron et al paper focuses on pancreatic cells where a droplet-based, single-cell RNA-seq method to determine the transcriptomes of over 12,000 individual pancreatic cells from four human donors and two mouse strains. The inDrop method provides a systematic approach for capturing thousands of cells without pre-sorting and uses high-throughput droplet microfluidics that barcode the RNA from thousands of individual cells, implementing a sensitive linear amplification method for single-cell transcriptome library construction. Cells were divided into 15 clusters that matched previously characterized cell types. Detailed analysis of each population separately revealed subpopulations within the ductal population, modes of activation of stellate cells, and heterogeneity in the stress among beta cells.

In this project, our attempt was to replicate the findings of the paper using current methodology and packages. The goals include processing the barcode reads of single cell RNAseq dataset, performing cell-by-gene quantification of UMI counts, performing quality control of UMI counts matrix, analyzing the UMI counts to identify clusters and marker genes for distinct cell type populations and finally providing a biological meaning to the clustered cell types followed by identification of novel marker genes associated with them.


# Contributors

Marina Natividad (Data Curator)

Vishwa Talati (Programmer)

Brad Fortunato (Analyst)

Kyrah Kotary (Biologist)

# Repository Contents
## Data Curation
Barcodecoutner.py is a python script for counting the appearance of the barcodes on each of the samples.
final_commands.qsub presents the bash script with the instructions for generating the human index as well as the UMI counts matrix.
## Programs for Programmer

### programmer.R

Script for filtering, normalization and clustering of genes.

* Input: quants_mat.gz which is output from ALevin

* Dependencies: R

* Execution: It is recommended to run this script in Rstudio by installing the following packages: Seurat, biomaRt, GenomicFeatures, tidyverse, Matrix, tximport, SeqGSEA, fishpond, EnsDb.Hsapiens.v79

   Alternatively to run on command line:

            module load R/4.0.2

            Rscript programmer.R

* Output: plots and panc_cells.rda file

