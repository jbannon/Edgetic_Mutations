import sys
import pickle 
import networkx as nx 
import NetworkCurvature as nc
import numpy as np 
from lifelines import KaplanMeierFitter
import tqdm
from sklearn.cluster import AffinityPropagation, KMeans
import pandas as pd

with open("data.pk","rb") as istream:
	dataset = pickle.load(istream)


tmbs = []
surv_times = []
surv_status = []

for i in tqdm.tqdm(range(len(dataset))):
	tmbs.append(dataset[i][4])
	surv_times.append(dataset[i][1])
	surv_status.append(int(dataset[i][2][0]))

X = np.array(tmbs).reshape(-1,1)
clustering = KMeans(n_clusters = 2,random_state=5).fit(X)

res_df = pd.DataFrame({"OS_MONTHS":surv_times,"OS_STATUS":surv_status, "CLUSTER":clustering.labels_})
res_df.to_csv("label_KM_tmb.csv")
sys.exit(1)


num_bins = 20

bins = np.linspace(-2,1,num_bins+1)

X = np.zeros((len(dataset),num_bins))
surv_times = []
surv_status = []

for i in tqdm.tqdm(range(len(dataset))):
	surv_times.append(dataset[i][1])
	surv_status.append(int(dataset[i][2][0]))
	G = dataset[i][0]
	orc = nc.OllivierRicciCurvature(G)
	orc.compute_edge_curvatures()
	G_ = orc.G.copy()
	x = np.array([j for j in nx.get_edge_attributes(G_,'ricci_curvature').values()])
	hist, bins = np.histogram(x,bins)
	X[i,:] = hist

clustering = AffinityPropagation(random_state=5).fit(X)

print(clustering.labels_)
res_df = pd.DataFrame({"OS_MONTHS":surv_times,"OS_STATUS":surv_status, "CLUSTER":clustering.labels_})
res_df.to_csv("label_AP.csv")
clustering = KMeans(n_clusters = 2,random_state=5).fit(X)
res_df = pd.DataFrame({"OS_MONTHS":surv_times,"OS_STATUS":surv_status, "CLUSTER":clustering.labels_})
res_df.to_csv("label_KM.csv")


