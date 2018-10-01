import glob
import os
import subprocess

curr_dir = os.getcwd()

individuals = [file.split('.')[0] for file in glob.glob(curr_dir) if file.endswith('_aln.sam')]

for ind in individuals:
    args_sam2bam = ['java', '-jar', 'picard.jar', 'SortSam', 'INPUT=' + ind + '_aln.sam', 'OUTPUT=' + ind + '_aln_sorted.bam', 'SORT_ORDER=coordinate']
    subprocess.call(' '.join(args_sam2bam))

    args_dedup = ['java', '-jar', 'picard.jar', 'MarkDuplicates', 'I=' + ind + '_aln_sorted.bam', 'O=' + ind + '_dedup.bam', 'M=' + ind + '_dedup_metrics.txt']
    subprocess.call(' '.join(args_dedup))

    args_index = ['java', '-jar', 'picard.jar', 'BuildBamIndex', 'INPUT=' + ind + 'dedup_reads.bam']
    subprocess.call(' '.join(args_index))
