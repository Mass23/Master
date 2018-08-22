#!/bin/bash
import subprocess
import os

files = [f for f in os.listdir(os.getcwd()) if f.endswith("paired.fq.gz")]

subprocess.call("bwa mem " + "fsel*.fasta " + " ".join(files) + " > alignment.sam", shell = True)
