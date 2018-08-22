#!/bin/bash
import subprocess

mdir = "M/trimmed_reads"
pdir = "P/trimmed_reads"

mref = "fsel_M.fasta"
pref = "fsel_P.fasta"

os.chdir(mdir)
subprocess.call("bwa index " + mref)

os.chdir(pdir)
subprocess.call("bwa index " + pref)
