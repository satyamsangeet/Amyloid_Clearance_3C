import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Load data
d1 = pd.read_csv("data/blattner2020_csf_concentration.csv")
d2 = pd.read_csv("data/lucey2018_csf_concentration.csv")
d3 = pd.read_csv("data/liu2022_csf_concentration.csv")
d4 = pd.read_csv("data/liu2022_plasma_concentration.csv")

d5 = pd.read_csv("data/blattner_ab40_concentration.csv")
d6 = pd.read_csv("data/lucey_ab40_concentration.csv")
d7 = pd.read_csv("data/liu_csf_ab40_concentration.csv")
d8 = pd.read_csv("data/liu_plasma_ab40_concentration.csv")

x = np.arange(0, 38, 2)
plt.figure(figsize=(10, 6))
plt.plot(
    x,
    #d5['Concentration'] / d1['Concentration'], 
    #d6['Concentration'] / d2['Concentration'], 
    d3['Concentration'] / d7['Concentration'], 
    #d8['Concentration'] / d4['Concentration'], 
    marker="o", 
    markersize=8, 
    color="teal", 
    linewidth=3
)

# Custom x-ticks: 0 to 36 in steps of 2
plt.xticks(ticks=np.arange(0, 37, 2), labels=np.arange(0, 37, 2), fontsize=15)
plt.yticks(fontsize=15)
plt.axvline(x=16, color='gray', linestyle='--', linewidth=2)
plt.axvline(x=24, color='gray', linestyle='--', linewidth=2)
#plt.savefig("liu_csf_ab40_ab42.png", dpi=300, bbox_inches="tight")
plt.show()
