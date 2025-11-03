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
- [kmer_results_fdr_with_blast_gene_clinicalclade.csv](../../Kmers/Data/kmer_results_fdr_with_blast_gene_clinicalclade.csv)
- [kmer_results_fdr_with_blast_gene_aquaticclade1.csv](../../Kmers/Data/kmer_results_fdr_with_blast_gene_aquaticclade1.csv)
- [kmer_results_fdr_with_blast_gene_aquaticclade2.csv](../../Kmers/Data/kmer_results_fdr_with_blast_gene_aquaticclade2.csv)
- [kmers_interactive_manhattan_plot_clinicalclade.py](../../Kmers/Scripts/kmers_interactive_manhattan_plot_clinicalclade.py)
- [kmers_interactive_manhattan_plot_aquaticclade1.py](../../Kmers/Scripts/kmers_interactive_manhattan_plot_aquaticclade1.py)
- [kmers_interactive_manhattan_plot_aquaticclade2.py](../../Kmers/Scripts/kmers_interactive_manhattan_plot_aquaticclade2.py)

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
