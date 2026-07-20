import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

d1 = pd.read_csv("data/blattner2020_csf_concentration.csv")
#d1 = pd.read_csv("data/lucey2018_csf_concentration.csv")
#d1 = pd.read_csv("data/liu2022_csf_concentration.csv")
#d1 = pd.read_csv("data/liu2022_plasma_concentration.csv")
d1.head()

x = np.arange(19)
labels = np.arange(0, 37, 2)

plt.figure(figsize=(18, 10))
plt.fill_between(x, d1["LSD"], d1["USD"], alpha=0.1, color="r")
plt.plot(x, d1["LSD"], color="k", linewidth=3, alpha=0.8)
plt.plot(x, d1["USD"], color="k", linewidth=3, alpha=0.8)
plt.plot(x, d1["Concentration"], linestyle="--", linewidth=5, color="r")
plt.scatter(x, d1["Concentration"], s=200, color="r", edgecolor="k", linewidth=3)
plt.axvline(x=8, linestyle="--", alpha=0.5, linewidth=3)
plt.axvline(x=12, linestyle="--", alpha=0.5, linewidth=3)
plt.xticks(ticks=x, labels=labels, fontsize=25)
plt.yticks(fontsize=25)
plt.xlim(0, 18)
plt.savefig("blattner_plot.png", dpi=300, bbox_inches="tight")
plt.show()
