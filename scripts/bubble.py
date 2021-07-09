import matplotlib.pyplot as plt
import pandas as pd

def parse(file):
    '''
        parse metadata
        '''
    dict1 = {} # {bft:aft}
    file1 = open(file)
    next(file1)
    for line in file1:
        temp = line.rstrip().split("\t")
        dict1[temp[0]] = temp[2]
    return dict1
#
def parse_file(file):
    '''
        parse SRRXXXX_comparison file
        '''
    dict1, agreement = {}, 0 # agreement
    file1 = open(file)
    next(file1)
    for line in file1:
        temp = line.rstrip().split("\t")
        disease = float(temp[1]) # disease (real)
        org_fmt = float(temp[3]) # aft-fmt (real)
        normalized_disease = float(temp[2]) # normalized disease
        driver_fmt = float(temp[4]) # adt (simulation)
        if (disease != 0) or (org_fmt != 0):
            real_trend = np.sign(org_fmt - disease) # real data abundance shift
            simulation_trend = np.sign(driver_fmt - normalized_disease) # simulated data abundance shift
            if real_trend == simulation_trend: # a) consistent b/t real and simulated data (agree)
                if real_trend > 0: # 1. increase
                    dict1[temp[0]] = ("red", normalized_disease*100)
                else: # 2. decrease
                    dict1[temp[0]] = ("blue", normalized_disease*100)
                agreement += normalized_disease
            else: # b) inconsistent
                dict1[temp[0]] = ("grey", normalized_disease*100)
    return dict1, agreement
#
adr = "/temp"
meta = parse("%s/metadata.txt"%adr)
tag = "tag"
# only driver species
files = sorted(glob.glob("%s/SRR*_comparison.%s.driver.txt"%(adr, tag)))
# disease + driver species
files2 = sorted(glob.glob("%s/SRR*_comparison.%s.disease.txt"%(adr, tag)))
pathogen = ["Streptococcus_agalactiae", "Campylobacter_jejuni"]
#
species1, data1, agree_dict = [], {}, {}
species2, data2 = [], {}
for file, file2 in zip(files, files2):
    srr = file.split("/")[-1].split("_")[0]
    sample = meta[srr]
    # drivers
    dict1, _ = parse_file(file)
    species1 += list(dict1.keys())
    data1[sample] = dict1
    # disease
    dict2, agreement = parse_file(file2)
    species2 += list(dict2.keys())
    data2[sample] = dict2
    # agreement
    agree_dict[sample] = agreement
#
species1 = sorted(["* %s"%x for x in set(species1) if x not in pathogen], reverse=True)
# remove driver species: only species in disease sample
species2 = sorted([x for x in set(species2) if "* %s"%x not in species1], reverse=True)
samples = sorted(list(data2.keys()))
species = species2 + species1
# Build dataframe
marker_dict = {"disease":"o", "driver":"X"}
df_dict = {"species":[], "samples":[], "Count":[], "color":[], "category":[]}
pallete = ["blue", "red", "lightgrey", "white"]
for sp in species:
    for sample in samples:
        df_dict["species"].append(sp)
        df_dict["samples"].append(sample)
        if sp in species1:
            if sp[2:] in data1[sample].keys():
                df_dict["Count"].append(30)
                df_dict["color"].append(data1[sample][sp[2:]][0])
            else:
                df_dict["Count"].append(0)
                df_dict["color"].append("white")
            df_dict["category"].append("driver")
        if sp in species2:
            if sp in data2[sample].keys():
                df_dict["Count"].append(data2[sample][sp][1])
                df_dict["color"].append(data2[sample][sp][0])
            else:
                df_dict["Count"].append(0)
                df_dict["color"].append("white")
            df_dict["category"].append("disease")

df = pd.DataFrame(df_dict)
df["padd"] = 200* (df.Count - df.Count.min()) / (df.Count.max() - df.Count.min())+5
###############Visualization###################
#Figure1 Bubble Plot
fig, ax = plt.subplots(figsize=(6, 11))
plt.rcParams['axes.facecolor'] = 'white'
for kind in marker_dict.keys():
    d = df[df.category==kind]
    if kind == "disease":
        scatter = ax.scatter(d.samples, d.species, \
                             s=d.padd, c=d.color,\
                             marker=marker_dict[kind])
    else:
        scatter = ax.scatter(d.samples, d.species, \
                             s=40, c=d.color,\
                             marker=marker_dict[kind], linewidths=0.05)
ax.grid(False)
plt.xticks(rotation=90)
plt.ylim([-0.5, len(species)-0.5])
plt.tight_layout()
plt.savefig('%s/bubble.%s.png'%(adr, tag), bbox_inches='tight', dpi=300)
plt.show()

#Figure2 Bar plot
# agreement y axis
agree_list = [agree_dict[sample] for sample in samples]
fig = plt.figure(figsize=(6, 1.5))
ax = fig.add_subplot(111)
ax.bar(samples, agree_list, color="white", edgecolor="black")
plt.xlim([-0.5, len(samples)-0.5])
plt.setp(ax.get_xticklabels(), visible=False)
ax.set_ylabel('Agreement (%)')
ax.set_facecolor("white")
plt.tight_layout()
plt.savefig('%s/agreement.%s.png'%(adr, tag), bbox_inches='tight', dpi=300)
plt.show()
