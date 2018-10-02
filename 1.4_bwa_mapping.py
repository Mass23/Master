#!/bin/bash

def bwa_map(ref):
    import subprocess
    import os

    files = [f for f in os.listdir(os.getcwd()) if f.endswith("_paired.fq.gz")]

    individuals_set = set()

    for i in files:
        name = i.split("_")[0]
        individuals_set.add(name)

    for individual in individuals_set:
        args_forward = ["cat", individual + "*forward_paired.fq.gz", ">", individual + "_forward.fastq.gz"]
        args_reverse = ["cat", individual + "*reverse_paired.fq.gz", ">", individual + "_reverse.fastq.gz"]
        subprocess.call(" ".join(args_forward), shell = True)
        subprocess.call(" ".join(args_reverse), shell = True)

        map_args = ["bwa", "mem", "-T", "8", "-M", "-R", "'" + r"@RG\tID:" + individual + r"\t'" ,ref, individual + "_forward.fastq.gz", individual + "_reverse.fastq.gz", ">", individual + "_alignment.sam"]
        subprocess.call(" ".join(map_args), shell = True)

bwa_map("fsel_M.fasta")
