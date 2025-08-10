# Maehara_et_al_2025

This code reproduces the figures from Maehara et al., *npj Systems Biology and Applications*, 2025,  
which identified changes in gene and protein expression associated with DNA methylation alterations  
in obese (ob/ob) mice.

## Usage
If you simply want to reproduce the figure, please run main.py.
The analysis pipeline, including table data formatting, is as follows:
1. Run RADMeth (install from https://smithlabresearch.org/software/radmeth/) to identify DMCpG with different DNA methylation ratios between WT and obese mice. Use bedtools (install from https://bedtools.readthedocs.io/en/latest/content/installation.html) to identify DMCpG near the gene transcription start site. *(Script: `methylomic_analysis/identify_DMCpG.sh`)*
2. Run edgeR to identify DEG. *(Script: `gene_expression_analysis/DEGanalysis.R`)*
3. Run limma to identify DEP. *(Script: `protein_expression_analysis/protein_normalization_and_analysis.R`)*
4. Perform intermediate data processing in MATLAB. *(Script: `main.m`)*
5. Generate the figures using Python.  *(Script: `main.py`)*
> **Note**: Matlab requires a paid license, so Matlab outputs are already included in the input folder. 
> Therefore, to reproduce the figures, you only need to run Step 5.
 
## Environment

Gene and protein expression analyses were performed using R version 4.3.1 with the following package versions:
edgeR=4.0.16, DEP=1.10.0, dplyr=1.1.2, readxl=1.4.2
Data processing were performed using **Matlab R2025a** on Apple Silicon
The Python code was executed in **Python 3.11.5** using the following version of the library:
numpy=1.24.3
