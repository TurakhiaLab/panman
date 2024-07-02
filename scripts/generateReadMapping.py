import matplotlib.pyplot as plt
import numpy as np


categories = ["SARS-CoV-2 (2,000)", "RSV (4,000)", "HIV (2,000)"]
panman_gfa = [11,16,17]
pggb_gfa = [27,313,51]


# Plotting with updated categories
x = np.arange(len(categories))

fig, ax = plt.subplots(figsize=(10, 5))

bar_width = 0.3
opacity = 0.8

shades = ['#000000', '#404040', '#808080', '#BFBFBF', '#E0E0E0', '#FFFFFF']

rects2 = ax.bar(x - 0.5*bar_width, pggb_gfa, bar_width, alpha=opacity, color=shades[1],edgecolor='black', label='PGGB-GFA')
rects4 = ax.bar(x + 0.5*bar_width, panman_gfa, bar_width, alpha=opacity,color=shades[4] ,edgecolor='black',hatch="/", label='PanMAN-GFA')



ax.set_xlabel('Species (Number of genomes)', fontsize=20, fontweight='bold')
ax.set_ylabel('File Size (MB)', fontsize=20, fontweight='bold')
# ax.set_title('Comparison of Different Metrics Across Datasets')
ax.set_xticks(x)
ax.set_xticklabels(categories, fontsize=16, fontweight='bold')
ax.tick_params(axis='y', labelsize=16, width=2, labelcolor='black', length=6)
ax.set_yscale('log')
ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.3), ncol=6, fontsize=16)
ax.margins(y=0.3)

fig.tight_layout()
plt.savefig("../figures/gfasize.svg", format="svg")
plt.savefig("../figures/gfasize.png")

categories = ["Sars-CoV-2 (2,000)", "RSV (4,000)", "HIV (2,000)"]
# panman_gfa = [89.76, 37.55, 53.32]
# pggb_gfa = [2.54, 17.92, 54.7]

panman_gfa = [100, 100, 99.7]
pggb_gfa = [100,100, 100]

x = np.arange(len(categories))

fig, ax = plt.subplots(figsize=(10, 5))

bar_width = 0.3
opacity = 0.8

shades = ['#000000', '#404040', '#808080', '#BFBFBF', '#E0E0E0', '#FFFFFF']

rects2 = ax.bar(x - 0.5*bar_width, pggb_gfa, bar_width, alpha=opacity, color=shades[1],edgecolor='black', label='PGGB-GFA')
rects4 = ax.bar(x + 0.5*bar_width, panman_gfa, bar_width, alpha=opacity,color=shades[4] ,edgecolor='black',hatch="/", label='PanMAN-GFA')


ax.set_xlabel('Species (Number of genomes)', fontsize=20, fontweight='bold')
ax.set_ylabel('Percentage of reads\nmapped', fontsize=20, fontweight='bold')
# ax.set_title('Comparison of Different Metrics Across Datasets')
ax.set_xticks(x)
yy=[0,25,50,75,100]
ax.set_yticks(yy)
ax.set_xticklabels(categories, fontsize=16, fontweight='bold')
ax.tick_params(axis='y', labelsize=16, width=2, labelcolor='black', length=6)
# ax.set_yscale('log')
ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.3), ncol=6, fontsize=16)
ax.margins(y=0.3)

fig.tight_layout()
plt.savefig("../figures/Mapping.svg", format="svg")

