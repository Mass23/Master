#!/bin/bash

#BSUB -L /bin/bash
#BSUB -o /scratch/beegfs/monthly/mbourqui/twisst/trimo.out
#BSUB -e /scratch/beegfs/monthly/mbourqui/twisst/trimo.err
#BSUB -J twisst_trim_reads
#BSUB -u massimo.bourquin@unil.ch
#BSUB -n 8
#BSUB -M 12000000

module add UHTS/Analysis/trimmomatic/0.36;
python3 -c "import os"
python3 -c "import glob"
python3 -c "import subprocess"

cd /scratch/beegfs/monthly/mbourqui/preprocess/M
python3 trim_reads.py

cd /scratch/beegfs/monthly/mbourqui/preprocess/P
python3 trim_reads.py

