import time 
from collections import defaultdict
import argparse
from typing import Dict, Union, List
import sys
import pandas as pd
import numpy as np 
import networkx as nx
import ot
import os
import yaml 
import NetworkCurvature as nc
import pickle
import tqdm



def unpack_parameters(
    D:Dict
    ):
    if len(D.values())>1:
        return tuple(D.values())
    else:
        return tuple(D.values())[0]


CANCER_TYPE_MAP = {
	'PANCAN':['PANCAN'],
	'BRCA':['Breast Cancer'],
	'LUNG':['Non-Small Cell Lung Cancer'],
	'SKCM':['Melanoma'],
	'BLCA':['Bladder Cancer'],
	'KIRC':['Renal Cell Carcinoma'],
	'ESO':['Esophagogastric Cancer'],
	'COAD':['Colorectal Cancer'],
	'GBM':['Glioma']}

def main(config:Dict):

	

	drug:str 
	tissue:str
	mutation_types:List[str]
	month_cutoff:int 
	exp_suffix:str
	vaf_cutoff:float
	primary_only:bool 
	mutation_file:str 
	sample_file:str
	patient_file:str
	network_file:str



	drug, tissue, mutation_types, month_cutoff, exp_suffix, vaf_cutoff,primary_only, num_bins =\
	 	unpack_parameters(config['EXPERIMENT_PARAMS'])

	mutation_file, sample_file, patient_file, network_file = unpack_parameters(config['FILES'])

	mutation_data = pd.read_csv(mutation_file,sep = "\t")
	
	mutation_data['VAF'] = mutation_data['t_alt_count']/(mutation_data['t_ref_count']+mutation_data['t_alt_count'])
	
	mutation_data = mutation_data[['Hugo_Symbol','Consequence','Variant_Classification','VAF','Tumor_Sample_Barcode']]

	
	with open(network_file, "rb") as istream:
		PPI_Graph = pickle.load(istream)
	primary_string = "PRIMARY" if primary_only else "METAST"

	cancer_types = CANCER_TYPE_MAP[tissue]
	path = f"../results/featurized/{primary_string}/{tissue}/{drug}/{exp_suffix}/"
	os.makedirs(path, exist_ok = True)
	mutation_data = mutation_data[mutation_data['Hugo_Symbol'].isin(PPI_Graph.nodes())]

	counter = 0
	mut2idx = {}
	mut_sub = mutation_data[mutation_data['Variant_Classification'].isin(mutation_types)]

	bins = np.linspace(-2,1,num_bins+1)
	for gene in pd.unique(mut_sub['Hugo_Symbol']):
		temp = mut_sub[mut_sub['Hugo_Symbol']==gene]
		for vclass in pd.unique(temp['Variant_Classification']):
			mut2idx[gene+"-"+vclass]= counter
			counter+=1
		
	sample = pd.read_csv(sample_file, sep = "\t",skiprows = 4)

	patient = pd.read_csv(patient_file, sep = "\t",skiprows = 4)
	
	patient = patient[patient['OS_MONTHS']>month_cutoff]
	
	if drug!="ALL":
		patient = patient[patient['DRUG_TYPE']==drug]

	if primary_only:
		sample = sample[(sample['SAMPLE_TYPE']=='Primary') & (sample['METASTATIC_SITE']=="Not Applicable")]

	if cancer_types[0]!= "PANCAN":
		sample = sample[sample['CANCER_TYPE'].isin(cancer_types)]

	samples = sample[['PATIENT_ID','SAMPLE_ID','TMB_NONSYNONYMOUS']]

	info = samples.merge(patient, on = 'PATIENT_ID')
	# info.reset_index(inplace=True, drop = True)
	

	vaf_vectors = []
	tmb_precomputed = []
	survival_status = []
	survival_duration = []
	curvature_features = []
	scaled_curvature_features = []
	

	for idx, row in tqdm.tqdm(info.iterrows(),total = info.shape[0]):
		
		patient_graph = PPI_Graph.copy()
		sample = row['SAMPLE_ID']
		survival_status.append(int(row['OS_STATUS'][0]))
		survival_duration.append(row['OS_MONTHS'])
		tmb_precomputed.append(row['TMB_NONSYNONYMOUS'])

		sample_mutations = mutation_data[mutation_data['Tumor_Sample_Barcode']==sample]
		potential_removals = sample_mutations[sample_mutations['Variant_Classification'].isin(mutation_types)]
		
		vaf_vect = np.zeros(len(mut2idx.keys()))
		genes_to_delete = []
		
		for rem_idx, rem_row in potential_removals.iterrows():
			vaf = rem_row['VAF']
			mut_key = rem_row['Hugo_Symbol']+"-"+rem_row['Variant_Classification']
			
			vaf_vect[mut2idx[mut_key]] = vaf
			if vaf>vaf_cutoff:
				genes_to_delete.append(rem_row['Hugo_Symbol'])
		

		vaf_vectors.append(vaf_vect)

		if len(genes_to_delete)>0:
			patient_graph.remove_nodes_from(genes_to_delete)
			
			# if not nx.is_connected(patient_graph):
			# 	largest_cc = sorted(list(max(nx.connected_components(this_G), key=len)))
			# 	drop_nodes = [x for x in this_G.nodes() if x not in largest_cc]
	
		orc = nc.OllivierRicciCurvature(patient_graph)	
		orc.compute_edge_curvatures()
		patient_graph = orc.G.copy()
		_,curv_vec =  zip(*nx.get_edge_attributes(patient_graph,'ricci_curvature').items())
		
		curv_feat, bins = np.histogram(curv_vec,bins)

		curvature_features.append(curv_feat)
		scaled_curvature_features.append(curv_feat/len(patient_graph.edges()))
		
		
		
	

	vaf_vectors = np.array(vaf_vectors)
	survival = np.hstack((np.array(survival_duration),np.array(survival_status)))
	curvature_features = np.array(curvature_features)
	scaled_curvature_features = np.array(scaled_curvature_features)
	tmb_precomputed = np.array(tmb_precomputed)


	np.save(path +"curvature.npy", curvature_features)
	np.save(path + "scaled_curvature.npy", scaled_curvature_features)
	np.save(path + "tmb_features.np",tmb_precomputed)
	np.save(path+"survival.npy",survival)
	np.save(path+"vaf_features.npy",vaf_vectors)
	

		



if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("-config",help="The config file for these experiments")
	args = parser.parse_args()
	
	with open(args.config) as file:
		config = yaml.safe_load(file)

	main(config)