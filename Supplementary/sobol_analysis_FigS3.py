import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

csv_file = "sobol_results_brain.csv"
#csv_file = "sobol_results_csf.csv"
#csv_file = "sobol_results_plasma.csv"

df = pd.read_csv(csv_file)
params = df["Parameter"]
Si = df["FirstOrder"]
ST = df["TotalOrder"]

interaction = ST - Si

x = np.arange(len(params))
width = 0.35

fig, axes = plt.subplots(2, 1, figsize=(12, 10))

# First Order vs Total Order
axes[0].bar(x - width/2, Si, width, label='First Order', edgecolor="k", linewidth=2)
axes[0].bar(x + width/2, ST, width, label='Total Order', edgecolor="k", linewidth=2)
axes[0].set_xticks(x)
axes[0].set_xticklabels(params, rotation=45)
axes[0].set_ylabel("Sensitivity Index")
axes[0].set_title("First Order vs Total Order Sobol Indices")
axes[0].legend()
axes[0].tick_params(axis="both", labelsize=25)
#axes[0].set_xticks([])
axes[0].grid(True, alpha=0.3, linestyle="--")

# First Order vs Interaction
axes[1].bar(x - width/2, Si, width, label='First Order', edgecolor="k", linewidth=2)
axes[1].bar(x + width/2, interaction, width, label='Interaction (ST - Si)', edgecolor="k", linewidth=2)
axes[1].set_xticks(x)
axes[1].set_xticklabels(params, rotation=45)
axes[1].set_ylabel("Sensitivity Index")
axes[1].set_title("First Order vs Interaction Effects")
axes[1].legend()
axes[1].tick_params(axis="both", labelsize=25)
#axes[1].set_xticks([])
#axes[1].set_ylim(0, 1.05)
axes[1].grid(True, alpha=0.3, linestyle="--")

plt.tight_layout()
#plt.savefig("sobol_barplots_brain.png", dpi=300, bbox_inches='tight')
plt.show()
