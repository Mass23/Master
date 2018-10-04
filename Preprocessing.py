#!/bin/bash
import subprocess
import glob
import multiprocessing

def trimmomatic(fastqfile, individual):

    if "_R1" in fastqfile:
        r1 = fastqfile
        r2 = r1.replace("_R1", "_R2")

    elif "_R2" in fastqfile:
        r2 = fastqfile
        r1 = r2.replace("_R2", "_R1")

    args = ["trimmomatic", "PE", "-threads 8", "-phred33",
            # Input R1, R2
            r1 , r2,
            # Output forward/reverse, paired/unpaired
            individual + "_forward_paired.fq.gz",
            individual + "_forward_unpaired.fq.gz",
            individual + "_reverse_paired.fq.gz",
            individual + "_reverse_unpaired.fq.gz",
            "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10",
            "LEADING: 3",
            "TRAILING: 3",
            "SLIDINGWINDOW:4:15",
            "MINLEN:36"]

    subprocess.call(' '.join(args), shell = True)

def bwa_map(individual, ref):

    # Merge forward and reverse files
    args_forward = ["cat", individual + "*forward_paired.fq.gz", ">", individual + "_forward.fastq.gz"]
    args_reverse = ["cat", individual + "*reverse_paired.fq.gz", ">", individual + "_reverse.fastq.gz"]

    subprocess.call(" ".join(args_forward), shell = True)
    subprocess.call(" ".join(args_reverse), shell = True)

    # Map reads
    map_args = ["bwa", "mem", "-T", "8", "-M", "-R", "'" + r"@RG\tID:" + individual + r"\t'" ,ref, individual + "_forward.fastq.gz", individual + "_reverse.fastq.gz", ">", individual + "_alignment.sam"]
    subprocess.call(" ".join(map_args), shell = True)

def picard_dedup(individual):

    # Sam to bam
    args_sam2bam = ['picard-tools', 'SortSam', 'INPUT=' + individual + '_alignment.sam', 'OUTPUT=' + individual + '_aln_sorted.bam', 'SORT_ORDER=coordinate']
    subprocess.call(' '.join(args_sam2bam), shell = True)

    # Mark duplicates
    args_dedup = ['picard-tools', 'MarkDuplicates', 'I=' + individual + '_aln_sorted.bam', 'O=' + individual + '_dedup.bam', 'M=' + individual + '_dedup_metrics.txt']
    subprocess.call(' '.join(args_dedup), shell = True)

    # Build bam index
    args_index = ['picard-tools', 'BuildBamIndex', 'INPUT=' + individual + '_dedup.bam']
    subprocess.call(' '.join(args_index), shell = True)

def realign_indels(individual, ref):


    # Index bam
    args_index = ['samtools', 'index', individual + '_dedup.bam']
    subprocess.call(' '.join(args_index), shell = True)

    # Intervals
    intervals = ['GenomeAnalysisTK', '-T', 'RealignerTargetCreator', '-R', ref, '-I', individual + '_dedup.bam', '-o', individual + 'forIndelRealigner.intervals']
    subprocess.call(' '.join(intervals), shell = True)
 
    # Realign indels
    realign = ['GenomeAnalysisTK', '-T', 'IndelRealigner', '-R', ref, '-I', individual + '_dedup.bam', '-targetIntervals', individual + 'forIndelRealigner.intervals', '-o', individual + '_realigned.bam']
    subprocess.call(' '.join(realign), shell = True) 

def call_snp(individual, ref):

    # Call snp
    call_args = ['GenomeAnalysisTK -T HaplotypeCaller -R', ref, '-I', individual + '_realigned.bam' , '-ERC GVCF -o', individual + '.g.vcf']
    subprocess.call(' '.join(call_args), shell = True)

################################################################
# 1. Index genome, list files
ref = 'fsel_M.fasta'
subprocess.call('bwa index ' + ref, shell = True)
subprocess.call('picard-tools CreateSequenceDictionary R=' + ref + ' O=' + ref.split('.')[0] + '.dict', shell = True)
subprocess.call('samtools faidx ' + ref, shell = True)

fastq_list = glob.glob('*.fastq.gz')

################################################################
# 2. Preprocess
def preprocess(fastqfile, ref):

    individual = fastqfile.split('_')[0]

    # 2.1 trimming
    trimmomatic(fastqfile, individual)

    # 2.2 mapping
    bwa_map(individual, ref)

    # 2.3 mark duplicates
    picard_dedup(individual)

    # 2.4 Realign indels
    realign_indels(individual, ref)

    # 2.4 call snp
    call_snp(individual, ref)

p = multiprocessing.Pool(64)
p.map(preprocess, fastq_list)
p.close()

