import glob, os
import numpy as np
from numpy import genfromtxt
from matplotlib.colors import ListedColormap
from skbio.stats.ordination import pcoa
from skbio.diversity import beta_diversity
import matplotlib.pyplot as plt
import pandas as pd
#
def otu_table(files):
    '''
        At last time point, species abundance
        '''
    table = [] # each sample
    for file in files:
        f2 = open(file, 'r')
        abd = [float(line.rstrip()) for line in f2]
        table.append(abd)
    return table

# Part1: Calculate PCoA matrix
adr = "/temp"
style = "input_pcoa"
layers = [1, 5, 10, 20, 50, 100]
repeat = 3 # repeat number
for layer in layers:
    for repeat in range(1, 4):
        files1 = glob.glob("%s/fmt_abd*"%adr)
        files2 = glob.glob("%s/disease_abd_*"%adr)
        files3 = glob.glob("%s/donor_abd_*"%adr)
        files = list(files1) + list(files2) + list(files3)
        # dX disease_id
        ids = ["aft_%s"%(x.split("_")[-1].split(".")[0]) for x in files1]
        ids += ["bft_%s"%(x.split("_")[-1].split(".")[0]) for x in files2]
        ids += ["donor_%s"%(x.split("_")[-1].split(".")[0]) for x in files3]
        table = otu_table(files)
        bc_dm = beta_diversity("braycurtis", table, ids)
        bc_pc = pcoa(bc_dm).samples
        bc_pc = bc_pc[bc_pc.columns[:4]] #only save top 3 PCoA
        bc_pc.to_csv("%s/%s_%s_%s.csv"%(adr, style, layer, repeat))
# Part2: Visualization
dataset, colormap = [], [] # each layer
for layer in layers:
    data, colors = [], []
    data = genfromtxt("%s/%s_%i_%i.csv"%(adr, style, layer, repeat), delimiter=',')
    data = data[1:, 1:3] #skip header, [PC1,PC2]
    custom_len = int(len(list(data))/3) # aft,bft,donor
    colors = [0]*custom_len + [1]*custom_len + [2]*custom_len
    dataset.append(data)
    colormap.append(colors)
#Figure2: PCoA
fig, axes = plt.subplots(2, 3, figsize=(20, 18), tight_layout=True)
axes = axes.flatten()
pallete = ["blue", "red", "green"]
for index in range(len(dataset)):
    df = dataset[index]
    ax1 = axes[index]
    ax1.set(xlabel=None, ylabel=None)
    ax1.set_title('%i-layer of Network'%layers[index], fontsize=18)
    scatter = ax1.scatter(df[:, 0], df[:, 1], c=np.asarray(colormap[index]), s=15,\
                          edgecolors="black", cmap=ListedColormap(pallete))
    handles1, labels1 = scatter.legend_elements(num=len(pallete))
    labels1 = ["ADT", "Disease", "Donor"]
    legend1 = ax1.legend(handles1, labels1, loc="best", \
                         title="Sample Type", borderpad=0.5, prop={"size":15}, title_fontsize=18, \
                         ncol=1, borderaxespad=0, frameon=False)
fig.text(0.5, 0.0, 'PCoA1', ha='center', fontsize=18)
fig.text(0.0, 0.5, 'PCoA2', va='center', rotation='vertical', fontsize=18)
plt.savefig('%s/%s_repeat%i.png'%(adr, style, repeat), dpi=300, bbox_inches='tight')
plt.show()
