import glob, csv, os
import pandas as pd
from collections import defaultdict
from micom.workflows import build
from micom import load_pickle
from micom.data import test_medium
from micom.qiime_formats import load_qiime_medium
from micom.workflows import grow
from micom.workflows import tradeoff
from micom import Community, data
from pathlib import Path
import logging
import argparse
import subprocess
import numpy as np
from pulp import *
import os.path

# MODE innteraction inference
def interaction_inference(input_file,model_address,medium_file,output,threshold,flag):
    #metabolic models
    models = glob.glob("%s/*"%model_address)
    # parameters
    handle = open("%s/micom.interaction.log"%output,"w")
    handle.write("Sample\tFraction of match\tUnfound Species\n")
    #
    filelist = [line.rstrip() for line in open(input_file)]
    for file in filelist:
        sample = Path(file).stem
        # total_frac: percentage of species > threshold
        # frac: percentage of species with metabolic models
        total_frac, frac = 0,0
        # ids: model id; model_files: model address; 
        # unmatched: list of species without metabolic models
        ids, model_files, species, abundance, unmatched = [], [], [], [], [] 
        file1 = open(file)
        for line in file1:
            temp = line.rstrip().split("\t")
            per = float(temp[1])
            if len(temp) >= 3:
                reads = float(temp[2])
            else:
                reads = 200*per*0.01
            if per > threshold:
                name = convert(temp[0])
                match = list(filter(lambda x: \
                        name == "_".join(x.split("/")[-1].split("_")[:2]), models))
                if len(match) > 0:
                    ids.append(match[0].split("/")[-1].split(".")[0])
                    species.append(name)
                    model_files.append(match[0])
                    abundance.append(reads)
                    frac += per
                else:
                    unmatched.append(name)
                total_frac += per
        handle.write("%s\t%f\t%s\n"%(sample,(frac/total_frac*100),",".join(unmatched)))
        data = {'id': ids, "species": species, \
                'file': model_files, "abundance": abundance}
        df = pd.DataFrame(data=data)
        build_community(df,output,sample,medium_file,flag)
    return

# covert name format
def convert(item):
    dict1 = {"Clostridioides_difficile":"Clostridium_difficile"}
    name = item.replace(" ","_")\
        .replace("]","").replace("[","")
    name = name.replace("Lachnoclostridium","Clostridium")
    if name in dict1.keys():
        name = dict1[name]
    return name

def build_community(df,output,prefix,medium_file,flag):
    # extract medium 
    medium = {}
    file2 = open(medium_file)
    next(file2)
    for line in file2:
        temp = line.rstrip().split(",")
        medium[temp[-1]] = float(temp[0])
    # save MICOM input information
    df.to_csv("%s/%s_micom_input.tsv"%(output,prefix), index=False)
    # build community 
    com = Community(df)
    com.medium = medium
    com.to_pickle("%s/%s_community.pickle"%(output,prefix))
    com = load_pickle("%s/%s_community.pickle"%(output,prefix))
    if flag:
        # community growth rate
        sol = com.optimize()
        gr = sol.members
        gr.to_csv("%s/%s_growth_rate.tsv"%(output,prefix))
    # bacterial interaction inference
    ko = com.knockout_taxa(fraction=0.5, method="relative change")
    ko.to_csv("%s/%s_ko.tsv"%(output,prefix), index=False)
    return

#mode MDSM
#STEP1 convert interaction matrix into undirected network
def parse_micom_undirected(file, strength, output):
    sample = str(Path(file).stem)
    network = defaultdict(set)
    file1 = open(file)
    first = file1.readline().rstrip()
    # species in the network
    species = ["_".join(x.split("_")[:2]) for x in first.split(",")]
    # write down undirected netowrk
    handle2 = open("%s/%s.str%s.undirected.txt"%(output,sample,str(strength).replace(".","")),"w")
    for count, line in enumerate(file1):
        sp1 = species[count]
        numbers = line.split(",") # list of interaction strength
        for index in range(len(numbers)):
            num = float(numbers[index])
            if (abs(num) > strength) and (index != count):
                sp2 = species[index]
                network[sp1].add(sp2)
                network[sp2].add(sp1)
    for k,v in network.items():
        handle2.write("%s\t%s\n"%(k,",".join(v)))
    return network
#STEP2 Find drivers
def MDSM(input_folder, strength, output, prefix):
    prob = LpProblem("DMS",LpMinimize) #objective 
    filelist = glob.glob("%s/*ko.tsv"%input_folder)
    V, networks = set(), []
    for file in filelist:
        g1 = parse_micom_undirected(file, strength, output)
        V = V.union(set(g1.keys()))
        networks.append(g1)
    x = LpVariable.dicts("nodes", V, lowBound=0, upBound=1, cat='Integer')
    handle = open("%s/%s.%ilayer.str%s.txt"%(output,prefix,len(filelist),str(strength).replace(".","")),"w")
    handle.write("@Total Number of Nodes: %i\n"%len(V))
    for counter, Vk in enumerate(networks):
        for vi in Vk.keys():
            temp = [x[i] for i in Vk[vi]] 
            temp.append(x[vi])
            prob += (lpSum(temp) >=1, "%s_constraint%i"%(vi,counter))
    obj = lpSum(x.values())
    prob += obj
    # Solve the optimization problem
    status = prob.solve()
    driver = [] 
    for var in prob.variables():
        if var.value() == 1:
            driver.append(str(var.name)[6:])
    handle.write("@Total number of drivers: %i\n"%len(driver))
    handle.write("\n".join(driver))
#
# FMT simulation: time-dependent GLV model
def glv_Euler_type(x,A,r):
    time = np.arange(0,200,0.01)
    # x initial values: abundance
    N = len(A[:,0])
    # nt: total time points
    nt = len(time)
    dt = time[1] - time[0];
    dxx = np.asarray(x)
    dx = np.asarray(x)
    # each time point
    for i in range(1,nt):
        tA = np.transpose(A)
        # A'*dxx: np.matmul(tA,np.transpose(dxx))
        # np.matmul(tA,np.transpose(dxx)) + r : row-wise
        # np.multiply(*,dxx) element-wise multiplication
        # dt: scaler
        dxx = dxx + np.multiply((np.matmul(tA,np.transpose(dxx)) + r),dxx)*dt
        # matrix -> single row array
        dx = np.vstack((dx,dxx))
    return dx, dxx
# FMT simulation: directed interaction networks (A)
def parse_micom(file, strength):
    A = []
    file1 = open(file)
    next(file1)
    for count, line in enumerate(file1):
        entity = []
        numbers = line.rstrip().split(",") # list of interaction strength
        intra_strength = float(numbers[count])
        for count2, num in enumerate(numbers):
            # intra_strength: -1
            if count2 == count:
                entity.append(num)
            #competition coefficient
            #inter_strength/intra_strength
            else:
                num = float(num)
                if abs(num) > strength:
                    entity.append(float(num/intra_strength))
                else:
                    entity.append(0)
        A.append(entity)
    cA = [list(map(float, lst)) for lst in A]
    return np.asarray(cA)
# FMT simullation: growth rate (r)
def parse_growth(file):
    file1 = open(file)
    next(file1)
    abs_r = []
    for line in file1:
        temp = line.rstrip().split(",")
        if "medium" == temp[0]:
            break
        abs_r.append(float(temp[2])) # growth rate
    #MICOM suggests divide predicted growth rate with
    #the overall weight of microbiota (200g)
    r = np.asarray([x/200 for x in abs_r]) # relative growth rate
    return r

#parse fmt_all input
#STEP1 parse input
def parse_fmt_input(input_file):
    # first column: disease sample
    # second column: donor sample
    file1 = open(input_file)
    diseases,donors = [],[]
    for line in file1:
        temp = line.rstrip().split(",")
        if len(temp) >= 2:
            diseases.append(temp[0])
            donors.append(temp[1])
    return diseases, donors

# pre_fmt bowel cleansing
def parse_taxa(file, antibiotic, threshold, abs_dict,model_address):
    models = glob.glob("%s/*"%model_address)
    f1 = open(file, 'r')
    match_frac, unmatch = 0, []
    for line in f1:
        temp = line.rstrip().split("\t")
        per = float(temp[1])
        if per > threshold:
            name = convert(temp[0])
            # check sp in models
            match = list(filter(lambda x: name == "_".join(str(Path(x).stem).split("_")[:2]) in x, models))
            if (len(match) > 0) and (name not in abs_dict.keys()):
                # strain_id, model_file, abundance
                strain = str(Path(match[0]).stem)
                abs_dict[name] = [strain, match[0], per/100*200*antibiotic]
                if (len(match) > 0) and (name in abs_dict.keys()):
                    abs_dict[name][2] += (per/100*200*antibiotic)
                    if len(match) > 0:
                        match_frac += per
                    else:
                        unmatch.append(name)
    unmatch_list = ",".join(unmatch)
    return abs_dict,match_frac,unmatch_list

#STEP2a build community (Donor)
def fmt_community(dis_sample,donor_sample,output,prefix,threshold, models):
    abs_dict = {}
    ids, model_files, abundance, species = [], [], [], []
    handle = open("%s/%s.micom.donor.fmt.log"%(output,prefix), "w")
    handle.write("Disease\tMatch%\tUnfound\tDonor\tMatch%\tUnfound\n")
    #disease metagenome: after bowel cleansing
    abs_dict, match_frac, unmatch = parse_taxa(dis_sample, 0.01, threshold, abs_dict, models)
    handle.write("%s\t%f\t%s\t"%(dis_sample,match_frac,unmatch))
    # fecal_transplant
    abs_dict, match_frac, unmatch = parse_taxa(donor_sample, 0.1, threshold, abs_dict, models)
    handle.write("%s\t%f\t%s\n"%(donor_sample,match_frac,unmatch))
    #convert to dataframe
    for sp, items in abs_dict.items():
        ids.append(items[0]) # strain
        model_files.append(items[1]) # model file address
        species.append(sp) #species
        abundance.append(items[2]) # abundance
    # Infer bacteria interaction network (MICOM)
    data = {'id': ids, "species": species, \
            'file': model_files, "abundance": abundance}
    df = pd.DataFrame(data=data)
    return df, abs_dict
#STEP2 build community (Driver)
def fmt_driver_community(dis_sample,driver_file,amount,output,prefix,threshold,model_address):
    abs_dict = {}
    ids, model_files, abundance, species = [], [], [], []
    handle = open("%s/%s.micom.driver.fmt.log"%(output,prefix), "w")
    handle.write("Disease\tMatch%\tUnfound\n")
    # disease metagenome: after bowel cleansing
    abs_dict, match_frac, unmatch = parse_taxa(dis_sample, 0.01, threshold, abs_dict, model_address)
    models = glob.glob("%s/*"%model_address)
    # all species in the disease sample
    patient_species = list(abs_dict.keys()) 
    handle.write("%s\t%f\t%s\n"%(dis_sample,match_frac,unmatch))
    drivers = [x.strip() for x in open(driver_file) if x[0] != "@"]
    #driver species transplantation
    for sp in drivers:
        name = convert(sp)
        match = list(filter(lambda x: name == "_".join(x.split("/")[-1].split("_")[:2]), models))
        # species not found in disease
        if (len(match) > 0) and (name not in patient_species):
            strain = match[0].split("/")[-1].split(".")[0]
            #add amount driver species
            # *0.1: pre-FMT bowel cleansing
            abs_dict[name] = (strain, match[0], amount*0.1)
    for sp, items in abs_dict.items():
        ids.append(items[0]) # strain
        model_files.append(items[1]) # model file address
        species.append(sp) #species
        abundance.append(items[2]) # abundance
    # Infer bacteria interaction network (MICOM)
    data = {'id': ids, "species": species, \
            'file': model_files, "abundance": abundance}
    df = pd.DataFrame(data=data)
    return df, abs_dict

#STEP3 Infer bacteria interaction
#STEP4 FMT process
def fmt_process(recovery,output,prefix,strength):
    A = parse_micom("%s/%s_ko.tsv"%(output,prefix), strength)
    # relative growth rate
    r = parse_growth("%s/%s_growth_rate.tsv"%(output,prefix))
    xx_fmt, x_fmt = glv_Euler_type(recovery,A,r)
    neg_count = len(list(filter(lambda x: (round(x,3) <= 0),x_fmt)))
    extreme_large = len(list(filter(lambda x: (x > 200),x_fmt)))
    if ("nan" in list(map(str, x_fmt))) or ("inf" in list(map(str, x_fmt))) or (neg_count == 0) or (extreme_large == 0):
        for item in xx_fmt[::-1]:
            neg_count = len(list(filter(lambda x: (round(x,3) <= 0), item)))
            extreme_large = len(list(filter(lambda x: (x > 200), item)))
            if (neg_count == 0) and (extreme_large == 0) \
                    and ("nan" not in list(map(str,item))) \
                    and ("inf" not in list(map(str,item))):
                x_fmt = item
                break
    np.savetxt('%s/fmt_abd_%s.txt'%(output,prefix),x_fmt,fmt='%.3f')
    np.savetxt('%s/fmt_timepoints_%s.txt'%(output,prefix),xx_fmt,fmt='%.3f')

## fmt_only
def fmt_only_process(input_file,output_folder,prefix, strength):
    abd_dict = {}
    file1 = open("%s_micom_input.tsv"%input_file)
    next(file1)
    for line in file1:
        temp = line.rstrip().split(",")
        abd_dict[temp[1]] = float(temp[-1])
    recovery = []
    for key in sorted(abd_dict.keys()):
        recovery.append(abd_dict[key])
    A = parse_micom("%s_ko.tsv"%input_file, strength)
    # relative growth rate
    r = parse_growth("%s_growth_rate.tsv"%input_file)
    xx_fmt, x_fmt = glv_Euler_type(recovery,A,r)
    neg_count = len(list(filter(lambda x: (round(x,3) <= 0),x_fmt)))
    extreme_large = len(list(filter(lambda x: (x > 200),x_fmt)))
    if ("nan" in list(map(str, x_fmt))) or ("inf" in list(map(str, x_fmt))) or (neg_count == 0) or (extreme_large == 0):
        for item in xx_fmt[::-1]:
            neg_count = len(list(filter(lambda x: (round(x,3) <= 0), item)))
            extreme_large = len(list(filter(lambda x: (x > 200), item)))
            if (neg_count == 0) and (extreme_large == 0) \
                    and ("nan" not in list(map(str,item))) \
                    and ("inf" not in list(map(str,item))):
                x_fmt = item
                break
    np.savetxt('%s/fmt_abd_%s.txt'%(output_folder,prefix),x_fmt,fmt='%.3f')
    np.savetxt('%s/fmt_timepoints_%s.txt'%(output_folder,prefix),xx_fmt,fmt='%.3f')
    return


##RUN PROGRAM##
def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="subparser_name", help='sub-command help')
    #
    interaction_command = subparsers.add_parser("interaction", help="Bacterial interaction inference using MICOM")
    interaction_command.add_argument("input", metavar="input_file", type=str, help="Input file of a list of taxonomic classification file addresses")
    interaction_command.add_argument("-m", "--medium", default="medium.csv", help="Medium CSV file, default medium.csv")
    interaction_command.add_argument("-d", "--model", default="dbs", help="Metabolic model database, default dbs")
    interaction_command.add_argument("-p", "--percentage", default=0.1, type=float, help="Percentage of species removed, default 0.1")
    interaction_command.add_argument("-f", "--flag", default=True, help="Calculate growth rate, default True")
    interaction_command.add_argument("-o", "--output", default="output_interaction", help="output folder, default output_interaction")
    #
    driver_command = subparsers.add_parser("driver", help="Driver nodes detection using MDSM")
    driver_command.add_argument("input", metavar="input_folder", type=str, help="Input folder of bacteria interaction networks")
    driver_command.add_argument("-s", "--strength", default=0.2, type=float, help="Threshold of interaction strength, default 0.2")
    driver_command.add_argument("-p", "--prefix", default="driver_nodes", help="Output file prefix, default driver_nodes")
    driver_command.add_argument("-o", "--output", default="output_driver", help="Output file folder, default output_driver")
    #
    fmt_command = subparsers.add_parser("fmt_donor", help="After-FMT community construction and simulation following the GLV model")
    fmt_command.add_argument("input", metavar="input_file", help="Input disease and donor sample file addresses")
    fmt_command.add_argument("-m", "--medium", default="medium.csv", help="Medium CSV file, default medium.csv")
    fmt_command.add_argument("-d", "--model", default="dbs", help="Metabolic model database, default dbs")
    fmt_command.add_argument("-p", "--percentage", default=0.1, type=float, help="Percentage of species removed, default 0.1")
    fmt_command.add_argument("-s", "--strength", default=0.2, type=float, help="Threshold of interaction strength, default 0.2")
    fmt_command.add_argument("-o", "--output", default="fmt_output", help="Output file folder, default fmt_output")
    #
    #
    fmt_driver_command = subparsers.add_parser("fmt_driver", help="Afte-driver species transplantation (ADT) community consturction and simulation following the GLV model")
    fmt_driver_command.add_argument("input", metavar="input_file", help="Input a list of disease sample file addresses")
    fmt_driver_command.add_argument("-i", "--driver", help="Input Driver Species")
    fmt_driver_command.add_argument("-a", "--amount", default=40, type=float, help="Input amount of driver species, default 40g")
    fmt_driver_command.add_argument("-m", "--medium", default="medium.csv", help="Medium CSV file, default medium.csv")
    fmt_driver_command.add_argument("-d", "--model", default="dbs", help="Metabolic model database, default dbs")
    fmt_driver_command.add_argument("-p", "--percentage", default=0.1, type=float, help="Percentage of species removed, default 0.1")
    fmt_driver_command.add_argument("-s", "--strength", default=0.2, type=float, help="Threshold of Interaction Strength, default 0.2")
    fmt_driver_command.add_argument("-o", "--output", default="fmt_driver_output", help="Output file folder, default fmt_driver_output")

    #
    fmt_only_command = subparsers.add_parser("fmt_only", help="After-FMT or ADT simulation following the GLV model")
    fmt_only_command.add_argument("input_file", metavar="input-file", help="Input file prefix")
    fmt_only_command.add_argument("-s", "--strength", default=0.2, type=float, help="Threshold of Interaction Strength, default 0.2")
    fmt_only_command.add_argument("-p", "--prefix", default="fmt_only_output", help="Output file prefix")
    fmt_only_command.add_argument("-o", "--output", default="fmt_output", help="Output file prefix")
    #
    args = parser.parse_args()
    #
    if args.subparser_name == "interaction":
        if os.path.exists(args.output) == False:
            os.mkdir(args.output)
        medium_file = args.medium
        flag = args.flag
        interaction_inference(args.input,args.model,medium_file,args.output,float(args.percentage),flag)
    #
    elif args.subparser_name == "driver":
        if os.path.exists(args.output) == False:
            os.mkdir(args.output)
        MDSM(args.input,float(args.strength),args.output,args.prefix)
    #
    elif args.subparser_name == "fmt_donor":
        if os.path.exists(args.output) == False:
            os.mkdir(args.output)
        output = args.output
        percentage = args.percentage
        medium_file = args.medium
        flag = True
        diseases,donors = parse_fmt_input(args.input)
        for file1, file2 in zip(diseases,donors):
            dis_sample = str(Path(file1).stem)
            don_sample = str(Path(file2).stem)
            prefix = "%s_%s"%(dis_sample,don_sample)
            df, abs_dict = fmt_community(file1,file2,output,prefix,args.percentage,args.model)
            build_community(df,output,prefix,medium_file,flag)
            # absolute abundance after immediate FMT
            recovery = [abs_dict[x][-1] for x in sorted(abs_dict.keys())]
            fmt_process(recovery,output,prefix,args.strength)
    #
    elif args.subparser_name == "fmt_driver":
        # build output folder
        if os.path.exists(args.output) == False:
            os.mkdir(args.output)
        output = args.output
        percentage = args.percentage
        medium_file = args.medium
        flag = True
        diseases = [x.strip() for x in open(args.input)]
        for file1 in diseases:
            dis_sample = str(Path(file1).stem)
            prefix = "%s_drivers"%dis_sample
            df, abs_dict = fmt_driver_community(file1,args.driver,args.amount,output,prefix,args.percentage,args.model)
            build_community(df,output,prefix,medium_file,flag)
            # absolute abundance after immediate FMT
            recovery = [abs_dict[x][-1] for x in sorted(abs_dict.keys())]
            fmt_process(recovery,output,prefix,args.strength)
#
    elif args.subparser_name == "fmt_only":
        if os.path.exists(args.output) == False:
            os.mkdir(args.output)
        fmt_only_process(args.input_file,args.output,args.prefix,args.strength)

main()









