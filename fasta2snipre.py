# Tested individuals fasta files finish with *1.fa and *2.fa corresponding to the haplotypes
# Outgroup sequence finishes with *A.fa
import glob
from Bio import SeqIO

class results:
    def __init__(self, pr, fr, ps, fs, tsil, trep, sample_size, identity,length,ncount):
        self.pr = pr
        self.fr = fr
        self.ps = ps
        self.fs = fs
        self.tsil = tsil
        self.trep = trep
        self.sample_size = sample_size
        self.identity = identity
        self.length = length
        self.ncount = ncount

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

    return(codon_table[codon])

def process_genes(individuals_list, gene_dict, outgroup_dict):

    gene_list = outgroup_dict.keys()
    results_dict = {}

    for gene in gene_list:

        gene_length = len(outgroup_dict[gene])

        for codon in range(0,gene_length,3):

            pos1 = codon
            pos2 = codon + 1
            pos3 = codon + 2

            outgroup_codon = outgroup_dict[gene][pos1] + outgroup_dict[gene][pos2] + outgroup_dict[gene][pos3]
            outgroup_amino_acid = get_amino_acid(outgroup_codon)

            population_codons = set()
            population_amino_acids = set()

            for individual in individuals_list:
                individual_codon = gene_dict[individual][gene][pos1] + gene_dict[individual][gene][pos2] + gene_dict[individual][gene][pos3]

                if 'N' or '*' or '-' in individual_codon:
                    continue
                else:
                    population_codons.add(individual_codon)
                    population_amino_acids.add(get_amino_acid(individual_codon))



        else:
            continue

# Create outgroup dictionary
outgroup = glob.glob('*_transcriptsA.fa')[0]
outgroup_dict = {}

for index, gene in enumerate(SeqIO.parse(open(outgroup), 'fasta')):
    if len(str(gene.seq)) // 3 == len(str(gene.seq)) / 3:

        gene_name = str(gene.id.replace('-RA','').replace('-RB','').replace('-RC','').replace('-RD','').replace('-RE','').replace('-RF','').replace('-RG','').replace('-RH','').replace('-RI','').replace('-RJ','').replace('-RK','').replace('-RL',''))

        other_transcripts = [i for i in outgroup_dict.keys() if i.startswith(gene_name)]

        try:
            last_length = len(outgroup_dict[other_transcripts[0]])
        except:
            last_length = 0

        if last_length > len(gene.seq):
            continue
        elif last_length == len(gene.seq):
            continue
        else:

            for k in list(outgroup_dict.keys()):
                if k.startswith(gene_name):
                    outgroup_dict.pop(k)

            outgroup_dict[gene.id] = str(gene.seq)
    else:
        continue

# Create individuals dictionaries
individual_files = list(glob.glob('*_transcripts1.fa')) + list(glob.glob('*_transcripts2.fa'))

gene_dict = {}

for individual in individual_files:
    individual_name = individual.replace('_transcripts1.fa','_1').replace('_transcripts2.fa','_2')

    for index, gene in enumerate(SeqIO.parse(open(individual), 'fasta')):
        if len(str(gene.seq)) // 3 == len(str(gene.seq)) / 3:
            if individual_name in gene_dict.keys():
                gene_dict[individual_name][gene.id] = str(gene.seq)
            else:
                gene_dict[individual_name] = dict()
                gene_dict[individual_name][gene.id] = str(gene.seq)
        else:
            continue

print('Outgroup: ', outgroup.replace('_transcripts1.fa',''), ' - ', len(outgroup_dict.keys()), ' genes')
for individual in gene_dict.keys():
    print('  - ', individual, ' - ', len(gene_dict[individual].keys()), ' genes')

# Process genes
genes_list = list(outgroup_dict.keys())

with open('snipre_input.csv', 'w') as output:
    nout = 27
    npop = 27

    output.write('\t'.join(['gene_id','PR','FR','PS','FS','Tsil','Trepl','nout','npop']) + '\n')

    count = 0
    total = len(genes_list) + 1

    for gene in genes_list:
        count += 1

        print('Gene: ', gene, ' - ', count, '/', total)
        outgroup_sequence = outgroup_dict[gene]

        individuals_sequence = []

        for individual in list(gene_dict.keys()):
            try:
                individuals_sequence.append(gene_dict[individual][gene])
            except:
                continue

        sample_size = len(individuals_sequence)
        gene_length = len(outgroup_sequence)

        ############################# MKR stats #############################

        dn = 0
        ds = 0
        pn = 0
        ps = 0
        tsil = 0
        trepl = 0
        identity = 0
        ncount = 0

        # 1. Get codons
        for i in range(0,gene_length,3):

            pos1 = i
            pos2 = i+1
            pos3 = i+2

            outgroup_codon = outgroup_sequence[pos1] + outgroup_sequence[pos2] + outgroup_sequence[pos3]

            if 'N' in outgroup_codon or '-' in outgroup_codon or '*' in outgroup_codon:
                ncount += 1
                continue
            else:
                outgroup_amino_acid = get_amino_acid(outgroup_codon)


            individuals_codons = set()
            individuals_amino_acids = set()

            for individual in individuals_sequence:
                try:
                    current_individual_codon = individual[pos1] + individual[pos2] + individual[pos3]
                except:
                    continue

                if 'N' in current_individual_codon or '-' in current_individual_codon or '*' in current_individual_codon:
                    continue
                else:
                    current_individual_amino_acid = get_amino_acid(current_individual_codon)
                    individuals_codons.add(current_individual_codon)
                    individuals_amino_acids.add(current_individual_amino_acid)

            # 2. Calculate stats

            # TSIL/TREP
            nucleotides = ['A', 'T', 'G', 'C']

            pos1, pos2, pos3 = outgroup_codon

            position1 = []
            position2 = []
            position3 = []

            for nuc in nucleotides:
                position1.append(str(get_amino_acid(nuc + pos2 + pos3)))
                position2.append(str(get_amino_acid(pos1 + nuc + pos3)))
                position3.append(str(get_amino_acid(pos1 + pos2 + nuc)))

            aa1 = len(position1) - position1.count(str(get_amino_acid(outgroup_codon)))
            aa2 = len(position2) - position2.count(str(get_amino_acid(outgroup_codon)))
            aa3 = len(position3) - position3.count(str(get_amino_acid(outgroup_codon)))

            trepl += (aa1/3) + (aa2/3) + (aa3/3)
            tsil += 3 - ((aa1/3) + (aa2/3) + (aa3/3))

            # Na or No change
            if len(individuals_codons) == 0:
                ncount +=1
            elif len(individuals_codons) == 1 and ''.join(individuals_codons) == outgroup_codon:
                identity += 1

            # Substitutions
            elif len(individuals_codons) == 1 and outgroup_codon not in individuals_codons:
                # Non-synonymous
                if ''.join(individuals_amino_acids) == outgroup_amino_acid:
                    ds += 1
                # Synonymous
                else:
                    dn += 1
            # Polymorphisms
            elif len(individuals_codons) > 1 :
                # Non-synonymous
                if len(individuals_amino_acids) > 1:
                    pn += 1
                # Synonymous
                elif len(individuals_amino_acids) == 1:
                    ps += 1
            else:
                print('PROBLEM:')
                print('gene: ', gene)
                print('outgroup:')
                print(outgroup_codon, outgroup_amino_acid)
                print('individuals:')
                print(individuals_codons, individuals_amino_acids)

        try:
            raw_mkr = (dn / ds) / (pn / ps)
        except:
            raw_mkr = 'NA'
        # Get summary stats
        output.write('\t'.join([str(gene),str(pn),str(dn),str(ps),str(ds),str(tsil),str(trepl),str(nout),str(npop)]) + '\n')
    
