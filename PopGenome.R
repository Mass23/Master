library(PopGenome)
library(ape)

# Create PopGenome dataset
GENOME.class <- readVCF('data/SC.vcf.gz', tid = 'Scaffold03', numcols = 20000, frompos = 1800000, topos = 12900000,  gffpath="data/f_selysi_M_v02.gff", include.unknown=TRUE)
GENOME.class <- set.synnonsyn(GENOME.class, ref.chr="fsel_M.fasta")
GENOME.class <- set.populations(GENOME.class,list(c("701W1","677W1","674W1","700W3","De115W1","De125W1","De172W1","De397W1","De434W1","De44W1","De67W1","De84W1"),c("De446W2","De299W1","De267W1","De259W2","710W2","716W1","722W1","733W1","748W1","715W2","508W1","706W4","703W4","750W2","713W3")), diploid = TRUE)

# McDonald-Kreitman test
genes <- splitting.data(GENOME.class, subsites="gene")
mkt <- MKT(genes)
mkt_results <- get.MKT(mkt)
mkt_results = as.data.frame(mkt_results)

# Select results to keep
mkt_results <- cbind(mkt_results,gene_pos=genes@region.names)
mkt_results$start <- sapply(strsplit(as.character(mkt_results$gene_pos),' - '), "[", 1)
keeps <- c("P2_nonsyn", "P2_syn", "D_nonsyn", "D_syn", "start")
mkt_results <- subset(mkt_results, select = keeps)

# Changes names
colnames(mkt_results)[which(names(mkt_results) == "P2_nonsyn")] <- "PR"
colnames(mkt_results)[which(names(mkt_results) == "P2_syn")] <- "PS"
colnames(mkt_results)[which(names(mkt_results) == "D_nonsyn")] <- "FR"
colnames(mkt_results)[which(names(mkt_results) == "D_syn")] <- "FS"

# Parse GFF
gff_fsel = read.gff("data/f_selysi_M_v02.gff")
gff_fsel = gff_fsel[which(gff_fsel$type=='gene'),]
gff_fsel$attributes = substr(gff_fsel$attributes, start = 4, stop = 14)

# Replace start by gene ID
mkt_results$gene_id = NA
for (row in 1:nrow(mkt_results)) {
  print(mkt_results[row, 'start'])
  gene_name= gff_fsel[which(gff_fsel$start == mkt_results[row, 'start']), 'attributes']
  mkt_results[row, 'gene_id'] <- gene_name
  }
mkt_results$start = NULL

# Add Tsil, Trep, outgroup and species chromosomes number
data_others = read.csv('tsil:trep.csv')
mkt_results$Tsil = NA
mkt_results$Trepl = NA
mkt_results$nout = 27
mkt_results$npop = 27

for (row in 1:nrow(mkt_results)) {
  print(mkt_results[row, 'gene_id'])
  mkt_results[row, 'Tsil'] = data_others[which(data_others$gene_id == mkt_results[row, 'gene_id']), 'Tsil']
  mkt_results[row, 'Trepl'] = data_others[which(data_others$gene_id == mkt_results[row, 'gene_id']), 'Trepl']
}
# FINAL DATASET
mkt_results = mkt_results[0:522,]

#################################################################
# SnIPRE
#################################################################
## Part (2)  Bayesian Implementation (R2WinBUGS package,
##          B_SnIPRE_source.R, and WinBUGS or OpenBUGS) necessary
#################################################################
source("B_SnIPRE_source.R")
source("my.jags2.R")
library(lme4)
library(R2jags)
library(arm)

BSnIPRE.run(mkt_results, burnin = 10000, thin = 4, iter = 15000)
load("samples")
res.mcmc <- samples
b.res <- BSnIPRE(res.mcmc, mkt_results)
bres = b.res$new.dataset
write.table(bres, file = "bayesian_results_SCP.csv", sep  = ",", row.names = FALSE)
