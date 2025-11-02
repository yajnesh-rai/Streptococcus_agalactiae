# Genome-Wide Association Study (GWAS)
This folder contains data and scripts to analyze GWAS results (SNPs, Indels and Kmers), map them to genes, and generate interactive Manhattan plots.

## SNPs and indels
### Extract SNPs

### Visualize interactive Manhattan plot 
#### Dependencies
- Python 3.10+
- pandas
- numpy
- plotly
- biopython

Install with:
```python
pip install pandas numpy plotly biopython
```
#### Files
- [pyseer_results_coresnp_fdr_annotation.txt](../Data/pyseer_results_coresnp_fdr_annotation.txt)
- [pyseer_results_indels_fdr_annotation.txt](../Data/pyseer_results_indels_fdr_annotation.txt)                                          
- [ref_NZ_CP012419.gb](../Data/ref_NZ_CP012419.gb)  
- [snp_indels_interactive_manhattan_plot.py](./snp_indels_interactive_manhattan_plot.py)

**Note:** Keep all files together — the script gets lonely otherwise!
#### Usage
navigate to the folder containing the script and data files, then run:
`````python
python snp_indels_interactive_manhattan_plot.py
`````
## Kmers
### Visualize interactive Manhattan plots
#### Dependencies
- Python 3.10+
- pandas
- numpy
- plotly

Install with:
```python
pip install pandas numpy plotly
```
#### Files
- [kmers_vs_ref_aquaticclade1_filtered.tsv](../../Kmers/Data/kmers_vs_ref_aquaticclade1_filtered.tsv)
- [kmers_vs_ref_aquaticclade2_filtered.tsv](../../Kmers/Data/kmers_vs_ref_aquaticclade2_filtered.tsv)
- [kmers_vs_ref_clinicalclade_filtered.tsv](../../Kmers/Data/kmers_vs_ref_clinicalclade_filtered.tsv)
- [ref_aquaticclade1_NZ_CP016501.gb](../../Kmers/Data/ref_aquaticclade1_NZ_CP016501.gb)
- [ref_aquaticclade2_NC_018646.gb](../../Kmers/Data/ref_aquaticclade2_NC_018646.gb)
- [ref_clinical_NC_007432.gb](../../Kmers/Data/ref_clinical_NC_007432.gb)
#### Usage
**Note:** Again. Dont breakapart the files — the script will start crying!
#### Usage
navigate to the folder containing the script and data files, then run it for all three references from each of the clades:
`````python
python kmers_interactive_manhattan_plot_clinicalclade.py
`````
`````python
python kmers_interactive_manhattan_plot_aquaticclade1.py
`````
`````python
python kmers_interactive_manhattan_plot_aquaticclade2.py
`````
