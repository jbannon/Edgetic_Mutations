import pickle 
import pandas as pd
import numpy as np 
import networkx as nx
from typing import List, Tuple, Dict
import tqdm



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
