#!/bin/sh
# LEADING: removes the leading bases under the threshold
# TRAILING: removes the trailing bases under the threshold
# MINLEN: drops the reads shorter than the threshold
# SLIDINGWINDOW: removes the k-mers under a quality threshold

# For paired reads: java -jar trimmomatic-0.35.jar PE -phred33 input_forward.fq.gz input_reverse.fq.gz output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
