# LEADING: removes the leading bases under the threshold
# TRAILING: removes the trailing bases under the threshold
# MINLEN: drops the reads shorter than the threshold
# SLIDINGWINDOW: removes the k-mers under a quality threshold
import os
import subprocess

def treat_individual(individual, M_or_P):
    dir_1 = os.chdir(M_or_P + os.path.sep + individual + os.path.sep)
    fastq_files = os.listdir(dir_1)
    individual = set()

    for i in fastq_files:
        i = i.split("." or "_")
        i[4].add(individual)
    
    for seq in individual:
        R1 = "R1_" + seq + ".fastq.gz"
        R2 = "R2_" + seq + ".fastq.gz"

        output = M_or_P + "_" + individual + "_" + seq + "_trim_"

        subprocess.call("java -jar trimmomatic-0.35.jar PE -phred33 " 
                        + "input_forward.fastq.gz " #TO CHANGE
                        + "input_reverse.fastq.gz " #TO CHANGE
                        + output + "forward_paired.fastq.gz " 
                        + output + "forward_unpaired.fastq.gz " 
                        + output + "reverse_paired.fastq.gz " 
                        + output + "reverse_unpaired.fastq.gz " 
                        + "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36")

individual_M_list = os.listdir("/M")

individual_P_list = os.listdir("/P")

for i in individual_M_list:
    treat_individual(i, "M")

for i in individual_P_list:
    treat_individual(i, "P")

print("Trimming done!")
