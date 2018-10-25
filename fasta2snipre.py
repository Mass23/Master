# Structure of the script arguments:
# python3 snipre.py annotation reference outgroup genome1 genome2 etc.
import subprocess
import sys
import glob
from Bio import SeqIO
import numpy
import operator

def get_amino_acid(codon):
    """ A function that returns the amino acid corresponding to a codon."""

    codon_table = {
    # A**
    'AAA' : 'LYS',
    'AAC' : 'ASN',
    'AAG' : 'LYS',
    'AAT' : 'ASN',

    'ACA' : 'THR',
    'ACC' : 'THR',
    'ACG' : 'THR',
    'ACT' : 'THR',

    'AGA' : 'ARG',
    'AGC' : 'SER',
    'AGG' : 'ARG',
    'AGT' : 'SER',

    'ATA' : 'ILE',
    'ATC' : 'ILE',
    'ATG' : 'MET',
    'ATT' : 'ILE',

    # C**
    'CAA' : 'GLN',
    'CAC' : 'HIS',
    'CAG' : 'GLN',
    'CAT' : 'HIS',

    'CCA' : 'PRO',
    'CCC' : 'PRO',
    'CCG' : 'PRO',
    'CCT' : 'PRO',

    'CGA' : 'ARG',
    'CGC' : 'ARG',
    'CGG' : 'ARG',
    'CGT' : 'ARG',

    'CTA' : 'LEU',
    'CTC' : 'LEU',
    'CTG' : 'LEU',
    'CTT' : 'LEU',

    # G**
    'GAA' : 'GLU',
    'GAC' : 'ASP',
    'GAG' : 'GLU',
    'GAT' : 'ASP',

    'GCA' : 'ALA',
    'GCC' : 'ALA',
    'GCG' : 'ALA',
    'GCT' : 'ALA',

    'GGA' : 'GLY',
    'GGC' : 'GLY',
    'GGG' : 'GLY',
    'GGT' : 'GLY',

    'GTA' : 'VAL',
    'GTC' : 'VAL',
    'GTG' : 'VAL',
    'GTT' : 'VAL',

    # T**
    'TAA' : 'STOP',
    'TAC' : 'TYR',
    'TAG' : 'STOP',
    'TAT' : 'TYR',

    'TCA' : 'SER',
    'TCC' : 'SER',
    'TCG' : 'SER',
    'TCT' : 'SER',

    'TGA' : 'STOP',
    'TGC' : 'CYS',
    'TGG' : 'TRP',
    'TGT' : 'CYS',

    'TTA' : 'LEU',
    'TTC' : 'PHE',
    'TTG' : 'LEU',
    'TTT' : 'PHE'}

    return codon_table[codon]

def get_nout_npop(reference_genome, individuals_genome):
    with open(reference_genome) as f:
        nout = len(SeqIO.parse(f, 'fasta'))
    with open(individuals_genome) as f:
        npop = len(SeqIO.parse(f, 'fasta'))
    return([nout,npop])

def get_geneset_size(file):
    gene_count = 0
    for i in SeqIO.parse(open(file), "fasta"):
        if len(i.seq)/ 3 == len(i.seq) // 3:
            gene_count +=1
    return(gene_count)

class gene:
    def __init__(self, individual, gene_id, gene_sequence):
        self.individual = individual
        self.gene_id = gene_id
        self.gene_sequence = gene_sequence
class result:
    def __init__(self, gene_id, pr, ps, fr, fs, tsil, trepl):
        self.gene_id = gene_id
        self.pr = pr
        self.ps = ps
        self.fr = fr
        self.fs = fs
        self.tsil = tsil
        self.trepl = trepl

##### 1. Parse arguments ####################################################################################################
print("Parsing ...")
arg_list = list(sys.argv)
annotation_file = arg_list[1]
reference_genome = arg_list[2]

# With the outgroup as first genome of the individuals list (to be on the top of the matrix)
outgroup_genome = arg_list[3]
individuals_genome = arg_list[3:]
tested_individuals = arg_list[4:]

##### 2. Get nout/npop ######################################################################################################
nout_npop = get_nout_npop(reference_genome, individuals_genome[0])
nout = nout_npop[0]
npop = nout_npop[1]

##### 3. Extract transcripts ################################################################################################
# Subprocess samtools faidx to index the reference genome and run subgffread to extract coding sequences
print("Indexing genomes, extracting transcripts ...")
args_faidx = ["samtools", "faidx", reference_genome]
subprocess.call(" ".join(args_faidx), shell = True)

for individual in individuals_genome:

    args_faidx = ["samtools", "faidx", individual]
    subprocess.call(" ".join(args_faidx), shell = True)

    args_subgffread = ("gffread", "-C", "-J", "-w", individual + "_transcripts.fa", "-g", individual, annotation_file)
    subprocess.call(" ".join(args_subgffread), shell = True)

transcript_files = glob.glob("*_transcripts.fa")

##### 4. Create genes matrix ################################################################################################
# Create empty matrix with width of the first fasta file length
data = numpy.empty(shape=[0, get_geneset_size(transcript_files[0])])

# Fill in the matrix with the gene class instances. Each iteration adds an individual on the next row with all the gene class instances on the columns.
print("Creating gene matrix ...")
for individual in transcript_files:

    name = individual.split(".")[0]
    gene_list = numpy.array(())
    
    with open(individual) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if len(record.seq)/ 3 == len(record.seq) // 3:
                gene_name = record.id
                gene_sequence = record.seq
                individual = name
                gene_list = numpy.hstack((gene_list, gene(individual, gene_name, gene_sequence)))
        
    data = numpy.vstack((data, gene_list))

##### 5. Prepare results output #############################################################################################
# Create a dictionary for the results, iterate over the columns, compare the sequences and calculate the metrics
results = {}

rows_number = numpy.shape(data)[0]
columns_number = numpy.shape(data)[1]
print("Data matrix created!")
print(rows_number, " individual(s)")
print(columns_number, " gene(s)")

##### 6. Analysis, print results #############################################################################################
def get_dndspnps(matrix, sequence_length):
    outgroup_sequence = matrix[0]
    individuals_sequences = matrix[1:]

    pr_count = 0
    fr_count = 0
    ps_count = 0
    fs_count = 0
    trepl = 0
    tsil = 0

    for i in range(0, sequence_length, 3):
        outgroup_codon = outgroup_sequence[i] + outgroup_sequence[i + 1] + outgroup_sequence[i + 2]

        # tsil/trep
        nucleotides = ['A', 'T', 'G', 'C']

        pos1, pos2, pos3 = outgroup_codon

        position1 = list()
        position2 = list()
        position3 = list()

        for nuc in nucleotides:
            position1.append(str(get_amino_acid(nuc + pos2 + pos3)))
            position2.append(str(get_amino_acid(pos1 + nuc + pos3)))
            position3.append(str(get_amino_acid(pos1 + pos2 + nuc)))

        aa1 = len(position1) - position1.count(str(get_amino_acid(outgroup_codon)))
        aa2 = len(position2) - position2.count(str(get_amino_acid(outgroup_codon)))
        aa3 = len(position3) - position3.count(str(get_amino_acid(outgroup_codon)))

        trepl += (aa1/3) + (aa2/3) + (aa3/3)
        tsil += 3 - ((aa1/3) + (aa2/3) + (aa3/3))

        # PR/FR/PS/FS
        individuals_codons = set()
        for individual in individuals_sequences:
            codon =  individual[i] + individual[i + 1] + individual[i + 2]
            individuals_codons.add(codon)
        
        # If No change or substitution
        if len(individuals_codons) == 1:
            
            individual_codon = str(''.join(individuals_codons))

            # No change
            if individual_codon == outgroup_codon:
                tsil += 3
            
            # Fixed
            if individual_codon != outgroup_codon:

                # Fixed replacement
                if get_amino_acid(individual_codon) == get_amino_acid(outgroup_codon):
                    fr_count += 1

                # Fixed silent
                elif get_amino_acid(individual_codon) == get_amino_acid(outgroup_codon):
                    fs_count += 1
        
        # If polymorphic
        elif len(individuals_codons) != 1:

            amino_acids_set = set()
            for i in individuals_codons:
                amino_acids_set.add(get_amino_acid(i))
            
            # If polymorphic silent
            if len(amino_acids_set) == 1:
                ps_count += 1
            
            # If polymorphic replacement
            elif len(amino_acids_set) != 1:
                pr_count += 1
        
    return([pr_count, fr_count, ps_count, fs_count, tsil, trepl])

results = {}

# Iterate over the genes
with open("snipre_input.csv", "w") as f:
    # Writes the header
    f.write("gene_id\tdn\tds\tpn\tps\ttsil\trepl\tnout\tnpop")

    for gene in range(0, columns_number):

        # Create matrix
        print("Gene: ", data[0][gene].gene_id, "(", gene + 1, "/",columns_number + 1 ,")")
        current_gene_size = len(data[0][gene].gene_sequence)
        current_gene_matrix = numpy.empty(shape=[0, current_gene_size])

        for individual in range(0, rows_number):
            individual_seq_list = [nucleotide for nucleotide in data[individual][gene].gene_sequence]
            current_gene_matrix = numpy.vstack((current_gene_matrix, individual_seq_list))
            
        # Analysis
        print("Analysis ...")
        mkr_stats = get_dndspnps(current_gene_matrix, current_gene_size)

        gene_results = result(data[0][gene].gene_id, mkr_stats[0], mkr_stats[1],mkr_stats[2],mkr_stats[3],mkr_stats[4],mkr_stats[5])

        results_line = "\t".join([str(gene_results.gene_id), str(gene_results.pr), str(gene_results.fr), str(gene_results.ps), str(gene_results.fs), str(gene_results.trepl), str(gene_results.tsil), str(nout), str(npop) + "\n"])
        f.writelines(results_line)

print("Counting done!")       
