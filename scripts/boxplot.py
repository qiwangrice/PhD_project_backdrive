import glob
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.stats import ranksums
# random
def num_meta(file):
    '''
        parse metafile: {layer number_repeat: number of drivers}
        output: count_drivers
        '''
    count_drivers = {}
    file2 = open(file)
    for line in file2:
        temp = line.rstrip().split(": ")
        count = int(temp[1]) # number of drivers
        count_drivers[temp[0]] = count
    return count_drivers

def parse_random(file):
    # rnd: number of random species is introduced
    rnd = int(file.split("/")[-1].split("_")[1]) #0-10
    pers = []
    file1 = open(file)
    next(file1)
    for line in file1:
        temp = line.rstrip().split("\t")
        disease = float(temp[1]) # disease sample C.diff abundance
        post_fmt = float(temp[2]) # aft-fmt C.diff abundance
        rec = max(-1, (disease - post_fmt)/(disease - 0)) # recovery degree
        pers.append(rec)
    return rnd, pers
#
def extract_data(file):
    '''
        parse driver species after-FMT results
        '''
    dict1 = defaultdict(list)
    file1 = open(file)
    next(file1)
    for line in file1:
        temp = line.rstrip().split("\t")
        num_fmt = int(temp[1]) # no. species introduced
        disease = float(temp[2]) # C. diff abundance in diseased
        post_fmt = float(temp[3]) # C. diff abundance in after-fmt
        recovery = max(-1, (disease - post_fmt)/(disease - 0))
        dict1[num_fmt].append(recovery)
    return dict1
#
def build_dataframe(dict1, rnd_dict):
    '''
        build pandas dataframe
        input: driver fmt dictionary, random fmt dictionary
        abundance: list of recovery degree of C. diff
        group: Random & MDSM
        group2: FMT size, no. of species
        '''
    abundance, group, group2 = [], [], []
    pvalues = []
    for num_fmt in sorted(dict1.keys()):
        # random data w/ same number of FMT
        rnd_rec = rnd_dict[num_fmt]
        driver_rec = dict1[num_fmt]
        # ranksum: pvalues
        _, pval = ranksums(rnd_rec, driver_rec)
        if pval < 0.05:
            pvalues.append(1)
        else:
            pvalues.append(0)
        #
        abundance += rnd_rec
        group += ["Random"]*len(rnd_rec)
        abundance += driver_rec
        group += ["MDSM"]*len(driver_rec)
        group2 += [num_fmt]*(len(rnd_rec)+len(driver_rec))
    data = {"FMT size, number of species":group2, "recovery degree":abundance, "Group":group}
    df = pd.DataFrame(data)
    return df, pvalues

# Figure1: Boxplot
adr = "/temp"
custom = "input_driver_fmt_result"
layers = [1, 5, 10, 20, 50, 100]
repeat = 1 # repeat no.
count_drivers = num_meta("%s/driver_metadata.txt"%adr)
files = glob.glob("%s/random_*"%adr)
rnd_dict = {} # {no. species: [recovery disease 1,...,recovery disease n]}
for file in files:
    rnd, pers = parse_random(file)
    rnd_dict[rnd] = pers
# Random Species
dataset, num_drivers = [], []
for layer in layers:
    file = "%s/%s_%i_%i.txt"%(adr, custom, layer, repeat)
    dict1 = extract_data(file) # {no. fmt species: [after-fmt C. diff abundance]}
    #dict1: driver fmt 1000 disease
    #rnd_dict: random fmt 1000 disease
    df, pvalues = build_dataframe(dict1, rnd_dict)
    dataset.append([df, pvalues])
    num_drivers.append(count_drivers["Layer_%i_repeat_%i"%(layer, repeat)])
# Visualization
fig, axes = plt.subplots(2, 3, figsize=(12, 15), tight_layout=True)
axes = axes.flatten()
for index, item in enumerate(dataset):
    df = item[0]
    ax1 = axes[index]
    sns.boxplot(ax=ax1, x=df["FMT size, number of species"],\
                y=df["recovery degree"],\
                hue=df["Group"], palette='Set2')
    ax1.legend(loc='lower left')
    ax1.set(xlabel=None, ylabel=None)
    ax1.set_title('%i-layer of Network (n=%i)'%(layers[index], num_drivers[index]), fontsize=14)
    pvalues = dataset[index][1]
    for index, pvalue in enumerate(pvalues):
        if pvalue == 1:
            # statistical annotation
            x1, x2 = index-0.2, index+0.2
            y, h, col = df['recovery degree'].max()*1.1, 0.001, 'k'
            ax1.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
            ax1.text((x1+x2)*.5, y+h, "*", ha='center', va='bottom', color=col)
fig.text(0.5, 0.0, "FMT size, number of species", ha='center', fontsize=13)
fig.text(0.0, 0.5, "Recovery Degree (μ)", va='center', rotation='vertical', fontsize=13)
plt.savefig('%s/fmt_boxplot_repeat%i'%(adr, repeat), dpi=300, bbox_inches='tight')
plt.show()
#
#


