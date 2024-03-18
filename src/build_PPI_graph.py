import matplotlib.pyplot as plt
import pickle 
import pandas as pd
import numpy as np 
import networkx as nx
from typing import List, Tuple, Dict
import tqdm
import sys 


def find_unambiguous_nodes(
	aliases:pd.DataFrame,
	measured_genes:List[str]
	) -> Tuple[List[str],List[str]]:

	common_genes =  list(set(measured_genes).intersection(set(aliases['alias'])))	
	ambiguous_genes = []
	good_genes = []
	good_string_ids = []
	empty_genes = []
	protein_gene_map = {}
	print("\n reconciling measured genes and STRING aliases\n")

	for gene in tqdm.tqdm(common_genes):
		temp = aliases[aliases['alias'] == gene]

		if len(pd.unique(temp['#string_protein_id']))>1:
			ambiguous_genes.append(gene)
		elif len(pd.unique(temp['#string_protein_id']))==1:
			good_genes.append(gene)
			# print(pd.unique(temp['#string_protein_id']))
			good_string_ids.append( list(pd.unique(temp['#string_protein_id']))[0])
			protein_gene_map[list(pd.unique(temp['#string_protein_id']))[0]] = gene
		elif temp.shape[0]==0:
			empty_genes.append(gene)
			print('empty')

	stat_string = "\n\t Started with {cg} genes\n\t Removed {ag} ambiguous genes\n\t {mg} missing genes\n\t keeping {gg}".format(
		cg = len(common_genes),
		ag = len(ambiguous_genes),
		mg = len(empty_genes),
		gg = len(good_genes))
	
	print(stat_string)
	
	

	return good_string_ids, good_genes, protein_gene_map



df = pd.read_csv('../data/tmb_mskcc_2018/data_mutations.txt',sep = "\t")


VALID_STRING_SOURCES = ['BLAST_KEGG_NAME','BioMart_HUGO']

STRING_aliases = pd.read_csv("../data/networks/STRING/9606.protein.aliases.v12.0.txt",sep="\t")
STRING_aliases = STRING_aliases[STRING_aliases['source'].isin(VALID_STRING_SOURCES)]
STRING_links = pd.read_csv("../data/networks/STRING/9606.protein.links.v12.0.txt",sep = ' ')


genes = list(pd.unique(df['Hugo_Symbol']))
good_pids, keep_genes, pid_gene_map = find_unambiguous_nodes(STRING_aliases,genes)

STRING_links = STRING_links[(STRING_links['protein1'].isin(good_pids)) & (STRING_links['protein2'].isin(good_pids))]
G = nx.Graph()

nodes1 = STRING_links['protein1'].values
nodes2 = STRING_links['protein2'].values
scores = STRING_links['combined_score'].values

for n1, n2, score in tqdm.tqdm(zip(nodes1, nodes2, scores),total = STRING_links.shape[0]):
	# print(row)
	# n1, n2, score = row['protein1'], row['protein2'], row['combined_score']
	node1 = pid_gene_map[n1]
	node2 = pid_gene_map[n2]
	if score >= 900:
		G.add_edge(node1, node2)


rem = [x for x in G.nodes() if x not in max(nx.connected_components(G), key=len)]

G.remove_nodes_from(rem)
with open("../data/networks/PPI_Graph.pk","wb") as ostream:
	pickle.dump(G,ostream)
