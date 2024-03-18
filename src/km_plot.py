import lifelines
import pandas as pd
from lifelines import KaplanMeierFitter
import matplotlib.pyplot as plt
import seaborn as sns
from lifelines.statistics import multivariate_logrank_test


kmf = KaplanMeierFitter()


df = pd.read_csv("label_KM.csv", index_col = 0)

ax = plt.subplot(111)
print(df['CLUSTER'].value_counts())
for clust in pd.unique(df['CLUSTER']):
	sub = df[df["CLUSTER"] == clust]

	kmf.fit(sub['OS_MONTHS'], event_observed=sub["OS_STATUS"], label="Cluster {d}".format(d=clust))
	kmf.survival_function_.plot(ax=ax)

# kmf.survival_function_.plot(ax=ax)
plt.title("Cluster KMs - Curve")
plt.savefig("curvature.png")
plt.close()

res = multivariate_logrank_test(df['OS_MONTHS'],df['CLUSTER'],df['OS_STATUS'])
res.print_summary()
df = pd.read_csv("label_KM_tmb.csv", index_col = 0)

ax = plt.subplot(111)
print(df['CLUSTER'].value_counts())
for clust in pd.unique(df['CLUSTER']):
	sub = df[df["CLUSTER"] == clust]

	kmf.fit(sub['OS_MONTHS'], event_observed=sub["OS_STATUS"], label="Cluster {d}".format(d=clust))
	kmf.survival_function_.plot(ax=ax)

# kmf.survival_function_.plot(ax=ax)
plt.title("Cluster KMs - TMB")
plt.savefig("tmb.png")
plt.close()
res = multivariate_logrank_test(df['OS_MONTHS'],df['CLUSTER'],df['OS_STATUS'])
res.print_summary()