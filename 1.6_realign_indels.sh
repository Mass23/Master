import glob
import os
import subprocess

curr_dir = os.getcwd()

files = glob.glob('*_dedup.bam')
print(files)

individuals = [file.split('_')[0] for file in files]

subprocess.call('picard-tools CreateSequenceDictionary R=fsel_M.fasta O=fsel_M.dict', shell = True)
subprocess.call('samtools faidx fsel_M.fasta', shell = True)

# Realign indels
for ind in individuals:
    # Index bam
    args_index = ['samtools', 'index', ind + '_dedup.bam']
    subprocess.call(' '.join(args_index), shell = True)

    # Intervals
    intervals = ['GenomeAnalysisTK', '-T', 'RealignerTargetCreator', '-R', 'fsel_M.fasta', '-I', ind + '_dedup.bam', '-o', ind + 'forIndelR$
    subprocess.call(' '.join(intervals), shell = True)
 
    # Realign indels
    realign = ['GenomeAnalysisTK', '-T', 'IndelRealigner', '-R', 'fsel_M.fasta', '-I', ind + '_dedup.bam', '-targetIntervals', ind + 'forIn$
    subprocess.call(' '.join(realign), shell = True)
 
    # Create fasta files
    args_2fasta = ['samtools', 'bam2fq', ind + '_realigned.bam', '|', 'seqtk', 'seq', '-A', '-', '>', ind + '.fasta']
    subprocess.call(' '.join(args_2fasta), shell = True)

