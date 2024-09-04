# Project Overview

This readme file contains information about the scripts used in this project.

## PART ONE
The analysis carried out in this part is done primarily through python scrips executed on the terminal

### create_chrom_level_fasta.py

The  project begins with the create_chrom_level_fasta.py file. It gets the needed reference genomes (lycopersicum and pennelli) and annotation files, creates the necessary directories, and uses SeqIO from the Bio library to read the reference genome files and write their records to a new folder in fasta format. Sequence records for the lycopersicum genome are written to the `Lyc` folder while sequence records for the pennelli genome are written to the `Pen` folder.

### rnaseq_mapper.py

This script performs the next step of the project. It constructs indexes for the lycopersicum and pennelli sequence records. The script uses the accession ids from the SRR_accessions.txt file to extract fastq raw sequencing data using fasterq-dump function from sra toolkit. Then it maps the raw reads to both the lycopersicum and pennelli genomes with hisat2, converts the output alignment files to sorted bam files, then extracts the reads with featureCounts. Required output files are savedd to either the `Lyc` folder or the `Pen` folder depending on whether the data is coming from the mapping to the lycopersicum genome or mapping to the pennelli genome. Extraneous files are then deleted.

### bamtobed_converter.py

This script takes the output bam files from the `Lyc` and  `Pen` folders and converts them to bed files for further processing and analysis

### bed_file_reader.py
This was an unsuccessful first attempt to read and process the bed files. This file can be skipped.

### bed_file_analyzer.py
This file successfully accomplishes what the bed_file_reader.py was supposed to. The first half of this script reads the set of bed files from the `Lyc` folder into a dictionary (`lyc_dataframes`) and reads the corresponding set of bed files from the `Pen` folder into another dictionary (`pen_dataframes`). It exctracts the keys (SRR accession) common to both dictionaries. Then for each pair of bed files, that is, a bed file coming from the `Lyc` folder and a corresponding bed file coming from the `Pen` folder (a pair of files with the same SRR accession), a union of the reads in both files is generated, and is saved to a new dictionary (`union_names_dict`) with the common key (common SRR accession id) serving as the key for each entry into the dictionary.  

The second half of the script uses a loop to iterate over all 523 unique SRR accessions (the common keys), initializes an empty dataframe with a `reads` column that is then filled the corresponding union of reads for that accession number. Two more columns, `lycopersicum` and `pennelli`  are added to the dataframe and initally filled with zeros. Then script checks if a read in the dataframe is in the `lyc_dataframes` dictionary or the `pen_dataframes` dictionary. If a read is found in the `lyc_dataframes` dictionary, the `lycopersicum` column is updated to 1, if it is found in the `pen_dataframes` dictionary, the `pennelli` column is updated to 1. Finally the dataframe is printed to a csv file

### txt_file_mover.py
This simple file moves all the counts.txt files obtained from running featureCounts to a new location to make further analysis easier. The counts.txt files from the `Lyc` folder are moved to a new folder (`Lyc/counts_folder`). The counts.txt files from the `Pen` folder are also moved to a separate folder (`Pen/counts_folder`). 

### counts_matrix_creator.py
For each folder (`Lyc/counts_folder` and `Pen/counts_folder` ), the script iterates over each counts file and reads it into a dataframe. It takes only the first and last columns, containing the genes and their counts respectively, and transposes them. Then it concatenates each successive dataframe with the preceeding one, creating a count matrix that contains the gene names as columns and the counts of each counts file as rows.

### filter_script.py
This script was an attempt to find which SRR file contains reads that mapped to lycopersicum only. It didn't work. This file can be skipped

### compare_filtered_files.py
This script accomplishes what the filter_script.py was meant to. It finds the identity of files that have no reads that mapped to pennelli at all. 

### get_lycopersicum_only.py
This script builds on the previous one. It processes the output of the compare_filtered_files.py and finds the files that have no pennelli reads at all.

### bam_mover.py
This script copies all already made bamfiles to a new folder

### index_creator.py
This script then generates index files for each of the bamfiles in their new location

### bam_filterer.py
This script makes use of the output files obtained from the `bed_file_analyzer.py` script to filterer the bamfiles so that only reads that mapped uniquely to lycopersicum are left in the file.

### filtered_bam_counter.py
This script uses featureCounts to generate the counts for the filtered bam files (the bamfiles that contain only reads that mapped uniquely to lycopersicum)

### filtered_counts_creator.py
This script generates a csv file counts matrix, containing the gene names as columns and each counts file as a row. 

### bam_mover_pen.py
This script copies all already made bamfiles to a new folder, but this time it's taking the bamfiles that were produced from mapping the raw fastq samples to the pennelli genome

### index_creator_pen.py
This script then generates index files for each of the pennelli bamfiles in their new location

### bam_filterer_pen.py
This script makes use of the output files obtained from the `bed_file_analyzer.py` script to filterer the bamfiles so that only reads that mapped uniquely to pennelli are left in the file.

### filtered_bam_counter_pen.py
This script uses featureCounts to generate the counts for the filtered pennelli bam files (the bamfiles that contain only reads that mapped uniquely to pennelli)

### filtered_counts_creator_pen.py
This script generates a csv file counts matrix, containing the gene names as columns and each counts file as a row. 

## PART TWO
The analysis carried out in this part is done with R and R studio

### R_introgression_analysis
The R_introgression_analysis folder contains three folders. `generate_cpms`, `generate_norm_cpms_no_na` and `identify_introgressed_genes`

### generate_cpms
The generate_cpms folder contains file `generate_lyc_cpm.R` which takes the csv file produced by the script `filtered_counts_creator.py` which contains the counts of reads that mapped uniquely to the lycopersicum genome and generates a corresponding lycopersicum cpm object named `lyc_cpm` and also a `lyc_gene_IDs` object.
Likewise, the folder also contains file `generate_pen_cpm.R` which takes the csv file produced by the script `filtered_counts_creator_pen.py` which contains the counts of reads that mapped uniquely to the pennelli genome and generates a corresponding pennelli cpm object named `pen_cpm` and also a `pen_geneIDs` object.

### generate_norm_cpms_no_na
This folder contains file `generate_norm_cpms_wo_na.R` which takes the lycopersicum cpm object and pennelli cpm object obtained from the corresponding `generate_lyc_cpm.R` and `generate_pen_cpm.R` and generates a normalized version of each object. The data is normalized with row means, but to avoid dividing by zero, row mean entries with values of zero are replaced with 1e-10. The file produces two new objects `norm_lyc_cpms_wo_na` and `norm_pen_cpms_wo_na`.

### identify_introgressed_genes
This is the main folder for this analysis. It contains file `introgressed_gene_matrix.R` file which uses the `norm_lyc_cpms_wo_na` object to identify which genes were introgressed from the pennelli genome.
Likewise, the folder also contains file `identify_pen_genes_that_got_introgressed.R` which uses the `norm_pen_cpms_wo_na` object to identify pennelli genes that jumped into the lycopersicum genome.
All the object required to run both R scripts are available in the R_objects folder. The outputs of both scripts are csv files with genes as rows and samples as columns, and cells with values of either 1 or 0 representing whether a gene is an introgressed gene or not.

### R_validation_analysis
This folder contains files `map_snp_to_nearest_gene.R`, which maps each snp to its closest gene, and `generate_validation_images.R` which generates validation images to compare introgressions as detected by RNA-seq vs introgressions as detected by SNP-marker data

### ML_train_lyc and ### ML_train_pen
These folders contain the  `utils.py` and `train.py` scripts that are used to build and train a convolutional neural network model on lycopersicum and pennellii genes respectively.

### ML_cross_predictions
This folder contains subfolders `apply_lyc_model_to_pen_genes` and `apply_pen_model_to_lyc_genes`, which apply the model trained on lycopersicum promoter-terminator sequences to pennellii genes and vice versa. 

### motif_discovery_modisco
This folder contains subfolders `lycopersicum` and `pennellii` which both contain the files `run_modisco.py` and `extract_motifs.py`. `run_modisco.py` is used to discover motifs the models use to make their predictions, and `extract_motifs.py` extracts those motifs.