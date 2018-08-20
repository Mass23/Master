# A script that processes alle the fastq files in a directory with trimmomatic.
# LEADING = removes the leading bases under the threshold
# TRAILING = removes the trailing bases under the threshold
# MINLEN = drops the reads shorter than the threshold
# SLIDINGWINDOW = removes the k-mers under a quality threshold
import os
#import subprocess
import glob

# List all the fastq files:
print(os.getcwd())
M = os.getcwd() + os.sep + "M" + os.sep 
P = os.getcwd() + os.sep + "P" + os.sep 

fastq_list = []

# Process all the files in the list and group them by pairs:
def process_individuals(dirMP):
    for i in os.walk(dirMP):
        fastq_list.append(i)
    
    for fastq_file in fastq_list:

        if fastq_file.endswith(".fastq.gz"):

            file_dir = str(dirMP) + "/" + fastq_file

            if "R1" in file_dir:
                r1_dir = file_dir
                r2_dir = r1_dir.replace("R1", "R2")
            
            elif "R2" in file_dir:
                r2_dir = file_dir
                r1_dir = r2_dir.replace("R2", "R1")
            
            fastq_list.remove(r1_dir)
            fastq_list.remove(r2_dir)

            print(r1_dir, r2_dir)

process_individuals(M)
process_individuals(P)
