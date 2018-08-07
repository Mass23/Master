#!/bin/bash
#
# initial example of a pipeline script

bwa aln -I -t 8 reference.fa s_1.txt > out.sai
bwa samse reference.fa out.sai s_1.txt > out.sam

samtools view -bSu out.sam  | samtools sort -  out.sorted
