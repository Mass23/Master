# Master - Balancing selection on a supergene controlling social organisation in the ant Formica selysi
#### Massimo Bourquin
## Summary
This pipeline processes whole-genome re-sequencing data to find signs of balancing selection in a socially polymorphic ant. The alpine silver ant (Formica selysi) can be monogynous (Sm allele) as well as polygynous (Sp allele) and this trait is genetically based on a supergene.

To find the signs of balancing selection, the following steps will be performed:

1. Pre-processing of the reads to have a good quality alignment of Sm/Sm and Sp/Sp individuals reads to the reference genome

2. Fst / Diversity / Tajima's D analysis

3. Whole-genome McDonald-Kreitman test to find genes under positive selection in both Sm and Sp

______________________________________________________________________________________________________________________________
## 1. Pre-processing
#### https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165
#### https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145

- Raw reads quality control
- Trim the adapters
- Map the reads to their respective (M or P) reference genome
- Mark duplicates
- Realign indels
- Get a clean sam file for M and one for P to use in the analyses

*******************************************************
### 1.1 Quality control - FastQC
#### https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

- Raw reads quality control

*******************************************************
###	1.2 Reads trimming - Trimmomatic
#### http://www.usadellab.org/cms/?page=trimmomatic

- Adapters trimming
- Remove leading and trailing low quality bases
- Cut low quality 4-mer
- Drop reads below the minimal length threshold

*******************************************************
### 1.3 Burrow-wheeler aligner and trimming - BWA
#### http://bio-bwa.sourceforge.net/

- Index the reference genome
- Map the reads against it
- Output in .sam format

*******************************************************
### 1.4 Duplicates marking - Picard
#### https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates

*******************************************************
### 1.5 Indels realignment - GATK
#### https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_indels_IndelRealigner.php


*******************************************************
### 1.6 Genotyping - GATK
#### https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_variantutils_GenotypeGVCFs.php

*******************************************************
### 1.7 Variant filtration BCFtools filter
#### https://samtools.github.io/bcftools/bcftools.html
______________________________________________________________________________________________________________________________
## 2. Fst / Diversity / Tajima's D analysis
### https://vcftools.github.io/index.html

______________________________________________________________________________________________________________________________
## 3. Positive selection inferences

- Use the annotation to extract the coding regions of the genome from the alignment
- Calculate the dN, dS, pN, Ps and other metrics needed by Snipre to find genes under positive selection
- Create the Snipre input file and launch the r code
- Compare the results to find genes under positive selection only in M or only in P

*******************************************************

### 3.1 Create individuals fasta and extract annotation
#### http://ccb.jhu.edu/software/stringtie/gff.shtml
*******************************************************

### 3.2 Calculate the dN, dS, pN, pS for each gene (Snipre input) - fasta2snipre.py
*******************************************************

### 3.3 Bayesian method for McDonald-Kreitman test - Snipre
#### https://bustamantelab.stanford.edu/lab-developed-software

*******************************************************
### 3.4 Ontology term analysis - Fisher exact test - Custom script



