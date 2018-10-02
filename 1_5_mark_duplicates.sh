
#!/bin/bash

#BSUB -L /bin/bash
#BSUB -N
#BSUB -o /scratch/beegfs/monthly/mbourqui/preprocess/dedup.out
#BSUB -e /scratch/beegfs/monthly/mbourqui/preprocess/dedup.err
#BSUB -J dedup
#BSUB -u massimo.bourquin@unil.ch
#BSUB -M 8000000
#BSUB -q "long"

module add UHTS/Analysis/samtools/1.8;
module add UHTS/Analysis/picard-tools/2.9.0;
module add UHTS/Analysis/GenomeAnalysisTK/4.0.4.0;

cd /scratch/beegfs/monthly/mbourqui/preprocess/M
python3 dedup.py
