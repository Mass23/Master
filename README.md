# Master - Genomic signs of balancing selection in a socially polymorphic ant
### Massimo Bourquin
## Summary


________________________________________________________________________________________________________________________________
## 1. Pre-processing
#### https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165
#### https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145

-----------------------------------------------------------
### 1.1 Quality control - FastQC
#### https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Code: https://github.com/Mass23/Master/blob/master/fastqc.py

-----------------------------------------------------------
###	1.2 Reads trimming - Trimmomatic
#### http://www.usadellab.org/cms/?page=trimmomatic

- Adapters trimming
- Remove leading and trailing low quality bases (below quality 3)
- Cut low quality 4-mer (below quality 15)
- Drop reads below the threshold (36 here)

Code: https://github.com/Mass23/Master/blob/master/reads_trimming.sh

-----------------------------------------------------------
### 1.3 Burrow-wheeler aligner and trimming - BWA
#### http://bio-bwa.sourceforge.net/

-----------------------------------------------------------
### 1.4 Duplicates marking - Picard
#### https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates

-----------------------------------------------------------
### 1.5 Indels realignment - GATK
#### https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_indels_IndelRealigner.php


________________________________________________________________________________________________________________________________
## 2. Phylogenomics

-----------------------------------------------------------
### 2.1 M and P individuals alignment - Mauve
#### http://darlinglab.org/mauve/mauve.html

-----------------------------------------------------------
### 2.2 SNP calling - Harvest
#### http://harvest.readthedocs.io/en/latest/content/harvest-tools.html
-----------------------------------------------------------

### 2.3 Fst, LD, Tajima D, GC content sliding-window - VCFtools
#### http://vcftools.sourceforge.net/man_latest.html
-----------------------------------------------------------

### 2.4 Phylogenetic tree - RAxML
#### https://sco.h-its.org/exelixis/software.html
-----------------------------------------------------------

### 2.5 Topology weighting - Twisst
#### https://github.com/simonhmartin/twisst

________________________________________________________________________________________________________________________________
## 3. Genomics signs of balancing selection
-----------------------------------------------------------

### 3.1 Intersect alignment and annotation - Bedtools intersect
#### http://bedtools.readthedocs.io/en/latest/content/tools/intersect.html
-----------------------------------------------------------

### 3.2 Calculate the dN, dS, pN, pS for each gene (Snipre input) - Custom script
-----------------------------------------------------------

### 3.3 Bayesian method for McDonald-Kreitman test - Snipre
#### https://bustamantelab.stanford.edu/lab-developed-software
-----------------------------------------------------------

### 3.4 Compare M and P results - Custom script
