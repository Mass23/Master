import glob
import os
import subprocess

list_files = glob.glob('*_alignment.sam')
individuals = [i.split("_")[0] for i in list_files]
print(individuals)

# Mark duplicates
for ind in individuals:
    args_sam2bam = ['picard-tools', 'SortSam', 'INPUT=' + ind + '_alignment.sam', 'OUTPUT=' + ind + '_aln_sorted.bam', 'SORT_ORDER=coordinate']
    subprocess.call(' '.join(args_sam2bam), shell = True)

    args_dedup = ['picard-tools', 'MarkDuplicates', 'I=' + ind + '_aln_sorted.bam', 'O=' + ind + '_dedup.bam', 'M=' + ind + '_dedup_metrics.txt']
    subprocess.call(' '.join(args_dedup), shell = True)

    args_index = ['picard-tools', 'BuildBamIndex', 'INPUT=' + ind + '_dedup_reads.bam']
    subprocess.call(' '.join(args_index), shell = True)
