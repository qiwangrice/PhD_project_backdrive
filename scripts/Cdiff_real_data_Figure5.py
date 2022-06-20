#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt
import pandas as pd
import seaborn.apionly as sns
from matplotlib.colors import LinearSegmentedColormap
from collections import defaultdict
import glob
import numpy as np
from skbio.stats.ordination import pcoa
import skbio
from skbio.stats.ordination import pcoa
from skbio.diversity import beta_diversity
import numpy as np
import pandas as pd
from skbio.stats.distance import anosim
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import ListedColormap

def parse(file):
    #parse metaddata.txt
    dict1 = {} 
    file1 = open(file)
    next(file1)
    for line in file1:
        temp = line.rstrip().split("\t")
        dict1[temp[0]] = temp[2]
    return dict1
#
def recovery_calculation(file,dict1):
    file1 = open(file)
    next(file1)
    output = {}
    for line in file1:
        temp = line.rstrip().split("\t")
        if temp[0] in pathogens:
            normalized_disease = float(temp[2])
            driver_fmt = float(temp[4])
            recovery= (normalized_disease-driver_fmt)/normalized_disease
            output[temp[0]]=recovery
            if recovery < 0: # set recovery degree <0 as 0
                print(temp[0],recovery)
                recovery=0
            dict1[temp[0]].append(recovery)
    return dict1,output

#Part1 Generate Figure5b

adr = "temp"
meta = parse("%s/metadata.txt"%adr)
# only drivers
files = sorted(glob.glob("%s/SRR*_comparison.*"%adr))
dict1=defaultdict(list)
pathogens=["Escherichia_coli","Klebsiella_oxytoca","Klebsiella_pneumoniae"]
handle=open("%s/output_recovery_degree.txt"%adr,"w")
handle.write("Sample\t%s\n"%("\t".join(pathogens)))
for file in files:
    srr = file.split("/")[-1].split("_")[0]
    sample = meta[srr]
    # drivers 
    dict1,output = recovery_calculation(file,dict1)
    handle.write("%s\t"%sample)
    for path in pathogens:
        if path in output:
            handle.write("%f\t"%output[path])
        else:
            handle.write("NA\t")
    handle.write("\n")
data,xaxis=[],[]
handle.write("Average\t")
for k,v in dict1.items():
    data.append(v)
    xaxis.append(k)
    handle.write("%f\t"%(sum(v)/len(v)))


#Figure 5b Bar plot
fig, ax1 = plt.subplots(figsize=(5, 5))
ax1.boxplot(data)
for i, d in enumerate(data):
    x = np.random.normal(i + 1, 0.04, len(d))
    ax1.scatter(x, d)
ax1.set_xticklabels(xaxis, rotation=-45)
ax1.set_xlabel("Pathogens", fontsize=13)
ax1.set_ylabel("Recovery Degree", fontsize=13)
plt.savefig('%s/real_data_recovery_degree'%adr, dpi=300, bbox_inches='tight')
plt.show()


# Part2 Generate PCoA


def parse_file(file):
    dict1={}
    file1=open(file)
    next(file1)
    total=0
    for line in file1:
        temp=line.rstrip().split("\t")
        dict1[temp[0]]=[float(x) for x in temp[1:]]
        total+=float(temp[-2]) #fmt
    for k,v in dict1.items():
        v[-2]=v[-2]/total*100
        dict1[k]=v
    return dict1


# only drivers
files = sorted(glob.glob("%s/SRR*_comparison.*"%adr))
dict2, species, samples = {}, [], []
for file in files:
    srr = file.split("/")[-1].split("_")[0]
    sample = meta[srr] 
    samples.append(sample)
    dict1=parse_file(file)
    species+=dict1.keys()
    dict2[sample]=dict1
species=sorted(list(set(species)))
data,labelss,ample_md,colors=[],[],{},[]
pallete = ["blue","red","green"]
for name in samples:
    sample=dict2[name]
    dis,fmt,drive=[],[],[] # length equals number of samples (12)
    for sp in species:        
        if sp in sample.keys():
            dis.append(sample[sp][1])
            fmt.append(sample[sp][2])
            drive.append(sample[sp][3])
        else:
            dis.append(0)
            fmt.append(0)
            drive.append(0)
    data.append(dis)
    data.append(fmt)
    data.append(drive)
    colors+=[0,1,2]
    labels+=["%s_disease"%name,"%s_fmt"%name,"%s_driver"%name]
#
sample_md = pd.DataFrame.from_dict(sample_md,orient="index")
bc_dm = beta_diversity("braycurtis", data, labels)
np.savetxt("%s/braycurtis_real_cdiff.csv"%adr,np.asarray(bc_dm.data),fmt='%.3f')
bc_pc = pcoa(bc_dm)
#PCoA
fig, ax1 = plt.subplots(figsize=(5, 5))
scatter=ax1.scatter(list(bc_pc.samples["PC1"]), list(bc_pc.samples["PC2"]), c = np.asarray(colors),                      s=25, edgecolors="black", cmap=ListedColormap(pallete))
handles1, labels1 = scatter.legend_elements(num=len(pallete))
labels1 = ["Disease","FMT", "Driver"]
legend1 = ax1.legend(handles1, labels1, loc="best", title="Condition", prop={"size":10},title_fontsize=10)
ax1.set_xlabel("PCoA1", fontsize=13)
ax1.set_ylabel("PCoA2", fontsize=13)
plt.savefig('%s/pcoa_read_cdiff.png'%adr,bbox_inches='tight', dpi=300)
plt.show()


