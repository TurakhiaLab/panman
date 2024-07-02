import matplotlib.pyplot as plt
import numpy as np



# Plotting
fig, axs = plt.subplots(2, 3, figsize=(14, 12))


# Flatten the array of axes
axs_flat = axs.flatten()


# SARS DATA
num_sequences = [20, 200, 2000, 20000]
GFA = [0.3, 2.1, 27, 1064.51]
VG = [0.075, 0.54, 7.4, 307]
GBZ = [0.12, 0.71, 5.5, 34]
Pangraph = [0.109, 1.7, 52, 1228.8]
PanMAN = [0.016, 0.042, 0.225, 2.0]


axs_flat[0].plot(num_sequences, GFA, '-o', label='GFA')
axs_flat[0].plot(num_sequences, VG, '-o', label='VG')
axs_flat[0].plot(num_sequences, GBZ, '-o', label='GBZ')
axs_flat[0].plot(num_sequences, Pangraph, '-o', label='Pangraph')
axs_flat[0].plot(num_sequences, PanMAN, '-o', label='PanMAN')

axs_flat[0].set_xscale('log')
axs_flat[0].set_yscale('log')

axs_flat[0].tick_params(axis='x', labelsize=16)
axs_flat[0].tick_params(axis='y', labelsize=16)

axs_flat[0].set_xticks(num_sequences)
axs_flat[0].set_xticklabels([str(n) for n in num_sequences])
axs_flat[0].set_title("SARS-CoV-2", fontsize=16,fontweight='bold')

#  HIV DATA
num_sequences = [20, 200, 2000, 20000]
GFA = [0.54, 1.8, 51, 1228.8]
VG = [0.21, 0.50, 19, 384]
GBZ = [0.195, 0.643, 13, 62.0]
Pangraph = [0.22, 2.6, 41, 451]
PanMAN = [0.03, 0.1, 1.5, 9.1]


axs_flat[1].plot(num_sequences, GFA, '-o', label='GFA')
axs_flat[1].plot(num_sequences, VG, '-o', label='VG')
axs_flat[1].plot(num_sequences, GBZ, '-o', label='GBZ')
axs_flat[1].plot(num_sequences, Pangraph, '-o', label='Pangraph')
axs_flat[1].plot(num_sequences, PanMAN, '-o', label='PanMAN')

axs_flat[1].set_xscale('log')
axs_flat[1].set_yscale('log')

axs_flat[1].tick_params(axis='x', labelsize=16)
axs_flat[1].tick_params(axis='y', labelsize=16)

axs_flat[1].set_xticks(num_sequences)
axs_flat[1].set_xticklabels([str(n) for n in num_sequences])
axs_flat[1].set_title("HIV", fontsize=16,fontweight='bold')

#  RSV DATA
num_sequences = [40, 400, 4000]
GFA = [1.1,18,313]
VG = [0.39,6,102]
GBZ = [0.14,0.5,2.4]
Pangraph = [0.28,5.4,108]
PanMAN = [0.03,0.1,0.46]


axs_flat[2].plot(num_sequences, GFA, '-o', label='GFA')
axs_flat[2].plot(num_sequences, VG, '-o', label='VG')
axs_flat[2].plot(num_sequences, GBZ, '-o', label='GBZ')
axs_flat[2].plot(num_sequences, Pangraph, '-o', label='Pangraph')
axs_flat[2].plot(num_sequences, PanMAN, '-o', label='PanMAN')

axs_flat[2].set_xscale('log')
axs_flat[2].set_yscale('log')

axs_flat[2].tick_params(axis='x', labelsize=16)
axs_flat[2].tick_params(axis='y', labelsize=16)

axs_flat[2].set_xticks(num_sequences)
axs_flat[2].set_xticklabels([str(n) for n in num_sequences])
axs_flat[2].set_title("RSV", fontsize=16, fontweight='bold')

# E.Coli DATA
num_sequences = [10, 100, 1000]
GFA = [43, 245, 2457.6]
VG = [14, 75, 747]
GBZ = [13, 94, 900]
Pangraph = [19, 173, 3481.6]
PanMAN = [3.1, 10, 114]


axs_flat[3].plot(num_sequences, GFA, '-o', label='GFA')
axs_flat[3].plot(num_sequences, VG, '-o', label='VG')
axs_flat[3].plot(num_sequences, GBZ, '-o', label='GBZ')
axs_flat[3].plot(num_sequences, Pangraph, '-o', label='Pangraph')
axs_flat[3].plot(num_sequences, PanMAN, '-o', label='PanMAN')
axs_flat[3].set_xscale('log')
axs_flat[3].set_yscale('log')
axs_flat[3].tick_params(axis='x', labelsize=16)
axs_flat[3].tick_params(axis='y', labelsize=16)
axs_flat[3].set_xticks(num_sequences)
axs_flat[3].set_xticklabels([str(n) for n in num_sequences])
axs_flat[3].set_title("E.Coli", fontsize=16, fontweight='bold')






#  TB DATA
num_sequences = [40, 400]
GFA = [72, 744]
VG = [22, 398]
GBZ = [28, 210]
Pangraph = [17, 402]
PanMAN = [1.6, 5]


axs_flat[4].plot(num_sequences, GFA, '-o', label='GFA')
axs_flat[4].plot(num_sequences, VG, '-o', label='VG')
axs_flat[4].plot(num_sequences, GBZ, '-o', label='GBZ')
axs_flat[4].plot(num_sequences, Pangraph, '-o', label='Pangraph')
axs_flat[4].plot(num_sequences, PanMAN, '-o', label='PanMAN')

axs_flat[4].set_xscale('log')
axs_flat[4].set_yscale('log')

axs_flat[4].tick_params(axis='x', labelsize=16)
axs_flat[4].tick_params(axis='y', labelsize=16)

axs_flat[4].set_xticks(num_sequences)
axs_flat[4].set_xticklabels([str(n) for n in num_sequences])
axs_flat[4].set_title("M. Tuberculosis", fontsize=16, fontweight='bold')


#  Klebs DATA
num_sequences = [10, 100, 1000]
GFA = [50, 470, 3891]
VG = [15, 139, 1228]
GBZ = [20, 183, 1536]
Pangraph = [21, 252, 5222]
PanMAN = [4.8, 22, 200]


axs_flat[5].plot(num_sequences, GFA, '-o', label='GFA')
axs_flat[5].plot(num_sequences, VG, '-o', label='VG')
axs_flat[5].plot(num_sequences, GBZ, '-o', label='GBZ')
axs_flat[5].plot(num_sequences, Pangraph, '-o', label='Pangraph')
axs_flat[5].plot(num_sequences, PanMAN, '-o', label='PanMAN')

axs_flat[5].set_xscale('log')
axs_flat[5].set_yscale('log')

axs_flat[5].tick_params(axis='x', labelsize=16)
axs_flat[5].tick_params(axis='y', labelsize=16)

axs_flat[5].set_xticks(num_sequences)
axs_flat[5].set_xticklabels([str(n) for n in num_sequences])
axs_flat[5].set_title("Klebsiella pneumoniae", fontsize=16, fontweight='bold')


def add_annotations(ax, num_sequences, data):
    ax.annotate(f"{data[0]}", (num_sequences[0], data[0]), textcoords="offset points", xytext=(-15,-10), ha='center')
    ax.annotate(f"{data[-1]}", (num_sequences[-1], data[-1]), textcoords="offset points", xytext=(15,10), ha='center')


fig.supxlabel('Number of Genome Sequences', fontsize=20, fontweight='bold', y=0.02)
fig.supylabel('File Size (MB)', fontsize=20, fontweight='bold', x=0.03)

legend_properties = {'weight':'bold', 'size':15}
fig.legend(['GFA', 'VG', 'GBZ', 'PanGraph', 'PanMAN'], loc='upper center', bbox_to_anchor=(0.5, 0.97), ncol=6, prop=legend_properties)


plt.savefig("../figures/scalingPlot.svg", format="svg")
