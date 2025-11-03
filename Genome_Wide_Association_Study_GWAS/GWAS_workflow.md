# Genome-Wide Association Study (GWAS)
This folder contains tools and scripts to perform Genome-Wide Association Studies (GWAS) for identifying genetic markers (SNPs, Indels, and K-mers) associated with phenotypic traits, and mapping them to genes.

## SNPs
#### Extract SNPS
Extract SNPs like a boss with snippy. A bigshoutout to Torsten Seemann for making bioinformatics a little less painful.

Install Snippy; 
The easy way
```console
conda install -c conda-forge -c bioconda -c defaults snippy
```
Try other ways https://github.com/tseemann/snippy

```python
snippy-multi /path_to_/input.tab --ref /path_to_/reference.fna --cpus 4 > runme.sh
```
```python
sh ./runme.sh
```
#### GWAS with Pyseer


```python
pyseer   --vcf /path_to_/core.vcf   --phenotypes /path_to_/phenotype.txt   --no-distances   > pyseer_results_snp.txt
```


## Indels
First, lets compress all VCF files per isolate
```python
for vcf in */snps.vcf; do
    bgzip -c "$vcf" > "${vcf}.gz"
    tabix -p vcf "${vcf}.gz"
done
```
Merge the compressed individual VCFs
```python
bcftools merge */snps.vcf.gz -O z -o merged.all.vcf.gz
bcftools index merged.vcf.gz
```
Extract the indels
```python
bcftools view -v indels -Oz -o merged.indels.vcf.gz merged.vcf.gz

bcftools index merged.indels.vcf.gz
```
Decompose multiallelic sites into separate lines
```python
bcftools norm -m -any -Oz -o merged.indels.decomposed.vcf.gz merged.indels.vcf.gz

bcftools index merged.indels.decomposed.vcf.gz
```
Simplify the format - GWAS compatible
```python
bcftools annotate -x ^FORMAT/GT -Oz -o indelsdecomposed.vcf.gz merged.indels.decomposed.vcf.gz
bcftools index indelsdecomposed.vcf.gz
```
#### GWAS with Pyseer


```python
pyseer   --vcf /path_to_/indelsdecomposed.vcf   --phenotypes /path_to_/phenotype.txt   --no-distances   > pyseer_results_snp.txt
```

## Kmers

