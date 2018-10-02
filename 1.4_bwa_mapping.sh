#!/bin/bash

#BSUB -L /bin/bash
#BSUB -N
#BSUB -o /scratch/beegfs/monthly/mbourqui/preprocess/map_Sal4W1.out
#BSUB -e /scratch/beegfs/monthly/mbourqui/preprocess/map_Sal4W1.err
#BSUB -J Sal4W1_mapping
#BSUB -n 8
#BSUB -M 12000000
#BSUB -u massimo.bourquin@unil.ch
#BSUB -q "long"

module add UHTS/Aligner/bwa/0.7.17;

cd /scratch/beegfs/monthly/mbourqui/preprocess/M
bwa mem -M -T 8 -R '@RG\tID:Sal\tSM:Sal4W1\tPL:illumina\tLB:lib1\tPU:unit1' fsel_M.fasta Sal4W1_L7_forward_paired.fq.gz Sal4W1_L7_reverse_paired.fq.gz > Sal4W1_alignment.sam

