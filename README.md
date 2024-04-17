# Danio rerio 45S-M variant Knockout raw data 
This repository contains scripts and datasets for the paper "Variant ribosomal DNA is essential for female differentiation in zebrafish, Moser TV  et al. (2024)".

## General Info
For identification of cas9 induced indel variants the custom scripts to analyse the raw data is provided.	
For analysis steps of raw PBAT data to infer 45S-M rDNA copy number shifts see methods for further details (trimming, mapping, read distribution assessment).		
Folder 'Script_Test' contains a small subset of the sequencing files from the "45S-M_KO_29dpf_01" experiment for the purpose of quickly testing whether the structure of the local environment, particularily the BBtools installation location works with the provided scripts.    	

## Prerequisits 
- The script requires the BBTools software suit (https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) to be installed under applications (mac).
- The script is currently designed for use on a Mac-OS system, running it on Windows or Linux will require adjustment of file paths.
- File structures needs to be unchanged for the cas9 indel variant pipeline to work without modification.

## Components 
### Scripts
- Bash scripts for cas9 indel variant indentification at the 45S-M locus in Danio rerio.    
- Bash scripts for listing indel variants detected in each read by sample.    
- R script for plotting cas9 indel variants.    

### Info files
- Moser_NZ_45SM-KO_Danio_SampleOverview.xlsx holds a list of all the raw sequencing files available in this repository.	
- 'HOW TO.txt' files in each amplicon data folder provide instructions for executing the scripts.
 

### Datasets
These are the raw sequencing files for various amplicon and PBAT sequencing analyses, "See Figure 1-4" for processed results.	

- 'Script_Test' contains a small subset of the sequencing files from the "45S-M_KO_29dpf_01" experiment for the purpose of quickly testing whether the structure of the local environment, 
  particularily the BBtools installation location works with the provided scripts.  
- '45S-M_KO_29dpf_01' holds sequencing files from samples taken 29dpf after fertilisation and consequent 45S-M knockout treatment.    	
  This folder is divided into PBAT and amplicon data subfolders, the later contains custom cas9 induced indel variant assessment scripts.     	
- '45S-M_KO_29dpf_03' holds sequencing files from samples taken 29dpf after fertilisation and consequent 45S-M knockout treatment.
  This folder is divided into PBAT and amplicon data subfolders, the later contains custom cas9 induced indel variant assessment scripts.     	
- '45S-M_KO_100dpf' holds sequencing files from samples taken 100dpf after fertilisation and consequent 45S-M knockout treatment.    	
  This folder contains only amplicon data, accompanied by the custom cas9 induced indel variant assessment scripts.     	
- '45S-M_KO_171dpf_SpermVSfinclip' holds sequencing files from samples taken 171dpf after fertilisation and consequent 45S-M knockout treatment. 		
  This folder contains only amplicon data, accompanied by the custom cas9 induced indel variant assessment scripts.     	
- '45S-M_KO_GuideTest_mRNAvsRNP' holds sequencing files from the initial 45S-M knockout guide verification and editing assesment.    	
    This folder contains only amplicon data (lig-amp), accompanied by the custom cas9 induced indel variant assessment scripts.
