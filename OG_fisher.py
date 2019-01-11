import pandas as pd
import scipy.stats as stats

# Parse annotation file
with open('f_selysi_M_v02.gff') as gff:
    data = pd.read_csv(gff, delimiter = '\t', header = None)

    ontologies_dict = {}

    for index, row in enumerate(data.itertuples()):
        if row._3 == 'gene':
            split = row._9.split(';')
            gene_id = split[0].replace('ID=','')
            ontologies = []
            for i in range(0,len(split)):
                if split[i].startswith('Ontology_term='):
                    go_terms = split[i].replace('Ontology_term=','').split(',')
                    ontologies.extend(go_terms)

            ontologies_dict[gene_id] = ontologies

#########################################################
#             |  Positively_selected | Not_significant  #
#------------------------------------|------------------#
# Go-term     |        Pos-GO        |       Neg-GO     #
# No GO-term  |       Pos-No-Go      |     Neg-No-Go    #
#########################################################

gene_list = []
ontologies_to_test = []
with open('gene_list_P.txt') as in_list:
    for line in in_list.readlines():
        gene_id = line.replace('\n','')
        gene_list.append(gene_id)
        ontologies_to_test.extend(ontologies_dict[gene_id])

ontologies_to_test = set(ontologies_to_test)

list_ontologies = list(ontologies_to_test)

cols = ['OntologyTerm','PosGO','NegGO','PosNoGo','NegNoGo','p-value','odds_ratio']
lst = []

for term in list_ontologies:

    ontology_term = str(term)

    PosGo = 0
    PosNoGo = 0

    NegGo = 0
    NegNoGo = 0

    for gene in ontologies_dict.keys():
        if term in ontologies_dict[gene]:
            if gene in gene_list:
                PosGo += 1
            else:
                NegGo += 1
        else:
            if gene in gene_list:
                PosNoGo += 1
            else:
                NegNoGo += 1

    oddsratio, pvalue = stats.fisher_exact([[PosGo, NegGo], [PosNoGo, NegNoGo]])

    lst.append([term, PosGo, NegGo, PosNoGo, NegNoGo, pvalue, oddsratio])

import statsmodels.stats.multitest as mltcorr

df1 = pd.DataFrame(lst, columns=cols)
reject, corr_p_value, alphas, alphab = mltcorr.multipletests(df1['p-value'], 0.05, 'fdr_bh')
df1['FDR_p-value'] = corr_p_value
# Output file
df1.to_csv('out_M.csv', sep = ',')
