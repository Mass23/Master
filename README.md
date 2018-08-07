# Master - Genomic signs of balancing selection in a socially polymorphic ant
### Massimo Bourquin
## Summary
This pipeline processes whole-genome re-sequencing data to find signs of balancing selection in a socially polymorphic ant. The alpine silver ant (Formica selysi) can be monogynous (M) as well as (P) polygynous and this trait is genetically based on a social chromosome. 

To find the traces of balancing selection, the following steps will be performed:

1. Pre-processing of the reads to have a good quality alignment of M and P individuals reads to their resepctive reference genome

2. Phylogenomic analysis using RAxML and the Twisst pipeline (topology weighting) to find genomic signs of phylogenetic discordance due to balancing selection.

3. Whole-genome McDonald-Kreitman test to find genes under positive selection in both M and P lineages and the ones under balancing selection.

______________________________________________________________________________________________________________________________
## 1. Pre-processing
#### https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165
#### https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145

*******************************************************
### 1.1 Quality control - FastQC
#### https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

- Raw reads quality control

Code: https://github.com/Mass23/Master/blob/master/fastqc.sh

*******************************************************
###	1.2 Reads trimming - Trimmomatic
#### http://www.usadellab.org/cms/?page=trimmomatic

- Adapters trimming
- Remove leading and trailing low quality bases
- Cut low quality 4-mer
- Drop reads below the minimal length threshold

Code: https://github.com/Mass23/Master/blob/master/reads_trimming.sh

*******************************************************
### 1.3 Burrow-wheeler aligner and trimming - BWA
#### http://bio-bwa.sourceforge.net/

- Index the reference genome
- Map the reads against it
- Output in .sam format

Code: https://github.com/Mass23/Master/blob/master/bwa_alignment.sh

*******************************************************
### 1.4 Duplicates marking - Picard
#### https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates

- Marks duplicates 

Code: https://github.com/Mass23/Master/blob/master/mark_duplicates.sh

*******************************************************
### 1.5 Indels realignment - GATK
#### https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_indels_IndelRealigner.php


______________________________________________________________________________________________________________________________
## 2. Genomics signs of balancing selection

*******************************************************

### 2.1 Intersect alignment and annotation - Bedtools intersect
#### http://bedtools.readthedocs.io/en/latest/content/tools/intersect.html
*******************************************************

### 2.2 Calculate the dN, dS, pN, pS for each gene (Snipre input) - Custom script
*******************************************************

### 2.3 Bayesian method for McDonald-Kreitman test - Snipre
#### https://bustamantelab.stanford.edu/lab-developed-software


______________________________________________________________________________________________________________________________
## 3. Phylogenomics

*******************************************************
### 3.1 M and P individuals alignment - Mauve
#### http://darlinglab.org/mauve/mauve.html

*******************************************************
### 3.2 SNP calling - Harvest
#### http://harvest.readthedocs.io/en/latest/content/harvest-tools.html
*******************************************************

### 3.3 Fst, LD, Tajima D, GC content sliding-window - VCFtools
#### http://vcftools.sourceforge.net/man_latest.html
*******************************************************

### 3.4 Phylogenetic tree - RAxML
#### https://sco.h-its.org/exelixis/software.html
*******************************************************

### 3.5 Topology weighting - Twisst
#### https://github.com/simonhmartin/twisst



### 3.4 Compare M and P results - Custom script
