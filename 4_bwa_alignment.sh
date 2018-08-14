#!/bin/bash

bwa mem -M fsel_M.fasta s_1.txt > out.sam

bwa mem -M fsel_P.fasta s_1.txt > out.sam
