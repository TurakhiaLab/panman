import matplotlib.pyplot as plt
import numpy as np


# Updated categories
categories = ["SARS-CoV-2\n(20,000)", "RSV\n(4,000)", "Escherichia\nColi\n(1,000)", "Mycobacterium\nTuberculosis\n(400)"]
panman_fitch = [2, 114, 0.46, 5]
panman_mppa = [4.46, 250.8, 1.84, 35]
gfa = [1064.51, 2457.6, 313, 744]
vg = [307, 747, 102, 398]
gbz = [34, 900, 2.4, 210]
pangraph = [1228.8, 3481.6, 108, 402]

# Plotting with updated categories
x = np.arange(len(categories))

fig, ax = plt.subplots(figsize=(16, 5))

bar_width = 0.15
opacity = 0.8

shades = ['#000000', '#404040', '#808080', '#BFBFBF', '#E0E0E0', '#FFFFFF']

rects1 = ax.bar(x - 2*bar_width, gfa, bar_width, alpha=opacity, color=shades[0],edgecolor='black', label='GFA')
rects2 = ax.bar(x - bar_width, vg, bar_width, alpha=opacity, color=shades[1],edgecolor='black', label='VG')
rects3 = ax.bar(x, gbz, bar_width, alpha=opacity, color=shades[2] ,edgecolor='black',label='GBZ')
rects4 = ax.bar(x + bar_width, pangraph, bar_width, alpha=opacity,color=shades[3] ,edgecolor='black', label='PanGraph')
rects5 = ax.bar(x + 2*bar_width, panman_mppa, bar_width, alpha=opacity, color=shades[4] ,edgecolor='black', hatch="/", label='PanMAN-MPPA')
rects6 = ax.bar(x + 3*bar_width, panman_fitch, bar_width, alpha=opacity,color=shades[5], edgecolor='black', label='PanMAN-Fitch')



ax.set_xlabel('Species\n(Number of genomes)', fontsize=20, fontweight='bold')
ax.set_ylabel('File Size (MB)', fontsize=20, fontweight='bold')
# ax.set_title('Comparison of Different Metrics Across Datasets')
ax.set_xticks(x)
ax.set_xticklabels(categories, fontsize=16, fontweight='bold')
ax.tick_params(axis='y', labelsize=16, width=2, labelcolor='black', length=6)
ax.set_yscale('log')
ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.3), ncol=6, fontsize=16)
ax.margins(y=0.3)

fig.tight_layout()
plt.savefig("../figures/LKPlot.svg", format="svg")


