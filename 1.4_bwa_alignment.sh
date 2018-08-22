#!/bin/bash
import subprocess

file_list = [i for i in glob.glob(os.getcwd()) if i.endswith(".fq.gz")]

subprocess.call("bwa mem " + "fsel*.fasta " + " ".join(mlist) + " > alignment.sam")
