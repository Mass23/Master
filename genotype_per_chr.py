# Example of use: python3 genotype_per_chr.py -r reference_genome -b bam_file -n cores_number
from multiprocessing import Pool
import subprocess

parser = argparse.ArgumentParser()

parser.add_argument('-r', '--ReferenceFasta', help='Reference genome fasta file used for mapping', type=str, action = 'store', required = True)
parser.add_argument('-b', '--BamFile', help='Bam file with aligned reads', type=str, action = 'store', required = True)
parser.add_argument('-n', '--NumberCores', help='Number of cores to use', type=str, action = 'store', required = True)

args = parser.parse_args()

ref_fasta = args.ReferenceFasta
bam_file = args.BamFile
n_cores = args.NumberCores

################################################################################

subprocess.call('bamtools split -in ' + bam_file + ' -reference', shell = True)

data_inputs = glob.glob('merged.REF_*.bam')

def process_scaffold(scaffold):
    args_freebayes = ['freebayes', '-f', ref_fasta, scaffold, '>', scaffold + '.vcf']

pool = Pool(n_cores)
pool.map(process_scaffold, data_inputs)
pool.close()
pool.join()
