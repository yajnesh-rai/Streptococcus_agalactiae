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
#### Usage
`````python
python interactive_manhattan_plot.py
`````
## Kmers
### Visualize interactive Manhattan plots
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
- [kmers_vs_ref_aquaticclade1_filtered.tsv](../../Kmers/Data/kmers_vs_ref_aquaticclade1_filtered.tsv)
- [kmers_vs_ref_aquaticclade2_filtered.tsv](../../Kmers/Data/kmers_vs_ref_aquaticclade2_filtered.tsv)
- [kmers_vs_ref_clinicalclade_filtered.tsv](../../Kmers/Data/kmers_vs_ref_clinicalclade_filtered.tsv)
- [ref_aquaticclade1_NZ_CP016501.gb](../../Kmers/Data/ref_aquaticclade1_NZ_CP016501.gb)
- [ref_aquaticclade2_NC_018646.gb](../../Kmers/Data/ref_aquaticclade2_NC_018646.gb)
- [ref_clinical_NC_007432.gb](../../Kmers/Data/ref_clinical_NC_007432.gb)
#### Usage
