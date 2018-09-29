# Structure of the script arguments:
# python3 fasta_to_dnds.py annotation reference outgroup genome1 genome2 etc.
import subprocess
import sys
import glob
from Bio import SeqIO
import numpy
import operator

def custom_division(number1, number2):
    try:
        out = float(number1 / number2)
    except:
        out = "Na"
    return(out)


class gene:
    def __init__(self, individual, gene_id, gene_sequence):
        self.individual = individual
        self.gene_id = gene_id
        self.gene_sequence = gene_sequence

class result:
    def __init__(self, gene_id, dn, ds, pn, ps, no_change, undetermined, polymorphism_count, substitution_count, a_content, c_content, g_content, t_content, n_content, length, identity):
        self.gene_id = gene_id
        self.dn = dn
        self.ds = ds
        self.pn = pn
        self.ps = ps
        self.dnds = custom_division(dn, ds)
        self.pnps = custom_division(pn, ps)
        self.mkr = custom_division(self.dnds, self.pnps)

        self.no_change = no_change
        self.undetermined = undetermined
        self.polymorphism_count = polymorphism_count
        self.polymorphism = custom_division(polymorphism_count, a_content + c_content + g_content + t_content)
        self.substitution_count = substitution_count
        self.substitution = custom_division(substitution_count, a_content + c_content + g_content + t_content)

        self.a_content = a_content
        self.c_content = c_content
        self.g_content = g_content
        self.t_content = t_content
        self.n_content = n_content
        self.gc_content = custom_division(g_content + c_content, a_content + c_content + g_content + t_content)
        self.uncertainty = custom_division(n_content, length)

        self.length = length
        self.identity = identity

# Parse arguments
print("Parsing ...")
arg_list = list(sys.argv)

annotation_file = arg_list[1]
reference_genome = arg_list[2]

# With the outgroup as first genome of the individuals list (to be on the top of the matrix)
outgroup_genome = arg_list[3]

individuals_genome = arg_list[3:]
tested_individuals = arg_list[4:]

# Subprocess samtools faidx to index the reference genome and run subgffread to extract coding sequences
print("Indexing genomes, extracting transcripts ...")
for individual in individuals_genome:

    args_faidx = ["samtools", "faidx", reference_genome]
    subprocess.call(" ".join(args_faidx), shell = True)

    args_subgffread = ("gffread", "-C", "-J", "-w", individual + "_transcripts.fa", "-g", individual, annotation_file)
    subprocess.call(" ".join(args_subgffread), shell = True)

transcript_files = glob.glob("*_transcripts.fa")

# Create empty matrix with width of the first fasta file length
gene_count = 0
for i in SeqIO.parse(open(transcript_files[0]), "fasta"):
    if len(i.seq)/ 3 == len(i.seq) // 3:
        gene_count +=1

data = numpy.empty(shape=[0, gene_count])

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

# Create a dictionary for the results, iterate over the columns, compare the sequences and calculate the metrics
results = {}

rows_number = numpy.shape(data)[0]
columns_number = numpy.shape(data)[1]
print("Data matrix created!")
print(rows_number, " individual(s)")
print(columns_number, " gene(s)")

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

def get_consensus_identity(matrix):
    data = matrix
    rows = numpy.shape(data)[0]
    columns = numpy.shape(data)[1]

    consensus = []
    identity_count = 0

    for nucleotide in range(0, columns):
        alleles_list = []
        identity_set = set()

        for individual in range(0, rows):
            alleles_list.append(data[individual][nucleotide])
            identity_set.add(data[individual][nucleotide])

        count_dict = {
            "A": alleles_list.count("A"),
            "C": alleles_list.count("C"),
            "G": alleles_list.count("G"),
            "T": alleles_list.count("T"),
            "N": alleles_list.count("N")}


        if len(identity_set) == 1:
            identity_count += 1
        
        consensus_allele = max(count_dict.items(), key = operator.itemgetter(1))[0]
        consensus.append(consensus_allele)
    
    return([consensus, float(identity_count / columns)])

def get_dndspnps(matrix, sequence_length):
    outgroup_sequence = matrix[0]
    individuals_sequence = matrix[1:]

    dn_count = 0
    ds_count = 0
    substitution_count = 0

    pn_count = 0
    ps_count = 0
    polymorphism_count = 0

    no_change = 0
    undetermined = 0

    for codon in range(0 ,sequence_length ,3):

        individuals_codon_set = set()
        individuals_amino_acid_set = set()

        # Create sets for the individuals codon and amino acids, extract the outgroup codon
        outgroup_codon = str(outgroup_sequence[codon] + outgroup_sequence[codon + 1] + outgroup_sequence[codon + 2])
        try:
            outgroup_amino_acid = get_amino_acid(outgroup_codon)
        # If the outgroup codon contains N's, the amino acid is unknown: NNN
        except:
            outgroup_amino_acid = "NNN"

        for individual in individuals_sequence:
            individual_codon = str(individual[codon] + individual[codon + 1] + individual[codon + 2])
            individuals_amino_acid = get_amino_acid(individual_codon)
            individuals_codon_set.add(individual_codon)
            individuals_amino_acid_set.add(individuals_amino_acid)

        # Substitution change 
        if len(individuals_codon_set) == 1:
            # No change
            if individuals_codon_set == outgroup_codon:
                no_change += 1
            # ds change
            elif individuals_amino_acid[0] == outgroup_amino_acid and individual_codon[0] != outgroup_codon:
                ds_count += 1
                substitution_count += 1
            # dN change
            elif individuals_amino_acid[0] != outgroup_amino_acid:
                dn_count += 1
                substitution_count += 1
            
        # Polymorphic change
        elif len(individuals_codon_set) > 1:
            polymorphism_count += 1

            if len(individuals_amino_acid_set) > 1:
                pn_count += 1
            else:
                ps_count += 1
            
    return([dn_count, ds_count, pn_count, ps_count, no_change, undetermined, polymorphism_count, substitution_count])

def get_nucleotide_count(matrix):
    data = matrix
    rows = numpy.shape(data)[0]
    columns = numpy.shape(data)[1]

    acount = 0
    ccount = 0
    gcount = 0
    tcount = 0
    ncount = 0

    for nucleotide in range(0, columns):
        for individual in range(0, rows):
            if data[individual][nucleotide] == "A":
                acount += 1
            elif data[individual][nucleotide] == "C":
                ccount += 1
            elif data[individual][nucleotide] == "G":
                gcount += 1
            elif data[individual][nucleotide] == "T":
                tcount += 1
            elif data[individual][nucleotide] == "N":
                ncount += 1
            else:
                print("Warning: character not recognised (A/C/G/T/N) !")

    return([acount, ccount, gcount, tcount, ncount])

results = {}
output_file = input("Output name: ")

# Iterate over the genes
with open(output_file + ".csv", "w") as f:
    f.write("gene_id\tdn\tds\tpn\tps\tdnds\tpnps\tmkr\tno_change\tundetermined\tpolymorphism_count\tpolymorphism\tsubstitution_count\tsubstitution\ta_content\tc_content\tg_content\tt_content\tn_content\tgc_content\tuncertainty\tlength\tidentity\n")

    for gene in range(0, columns_number):

        # Define the size of the current gene and create an empty matrix of the gene set width
        print("Gene: ", data[0][gene].gene_id, "(", gene + 1, "/",columns_number + 1 ,")")

        current_gene_size = len(data[0][gene].gene_sequence)
        current_gene_matrix = numpy.empty(shape=[0, current_gene_size])

        for individual in range(0, rows_number):
            individual_seq_list = [nucleotide for nucleotide in data[individual][gene].gene_sequence]
            current_gene_matrix = numpy.vstack((current_gene_matrix, individual_seq_list))
            
        # Get consensus sequence
        print("Getting consensus sequence ...")
        consensus_sequence = get_consensus_identity(current_gene_matrix)

        # Get nucleotide count
        print("Getting nucleotide counts ...")
        nucleotide_count = get_nucleotide_count(current_gene_matrix)

        # Get MKR stats
        print("Getting MKR stats ...")
        mkr_stats = get_dndspnps(current_gene_matrix, current_gene_size)

        gene_results = result(data[0][gene].gene_id, mkr_stats[0], mkr_stats[1], mkr_stats[2], mkr_stats[3], mkr_stats[4], mkr_stats[5], mkr_stats[6], mkr_stats[7], nucleotide_count[0], nucleotide_count[1], nucleotide_count[2], nucleotide_count[3], nucleotide_count[4], current_gene_size, consensus_sequence[1])
        results[data[0][gene].gene_id] = gene_results

        results_line = "\t".join([str(gene_results.gene_id), str(gene_results.dn), str(gene_results.ds), str(gene_results.pn), str(gene_results.ps), str(gene_results.dnds), str(gene_results.pnps), str(gene_results.mkr), str(gene_results.no_change), str(gene_results.undetermined), str(gene_results.polymorphism_count), str(gene_results.polymorphism), str(gene_results.substitution_count), str(gene_results.substitution), str(gene_results.a_content), str(gene_results.c_content), str(gene_results.g_content), str(gene_results.t_content), str(gene_results.n_content), str(gene_results.gc_content), str(gene_results.uncertainty), str(gene_results.length), str(gene_results.identity) + "\n"])
        f.writelines(results_line)

print("Analysis done!")

       
