from numpy import arange
from array import *
import math
import random
import numpy as np
def generate_A(N,C,delta,diag,VarianceType):
    #N: no. of species
    #C: probability of effect from species-i to species-j
    # delta: vaariaance of a_ij
    # diag = -1.0 a_ii
    r = np.random.uniform(low=0, high=1, size=(1,N))
    A = np.zeros((N,N))
    if VarianceType == 1:
        zigma = 1/((N*(2+delta))**(1/2))
    elif VarianceType == 2:
        zigma = delta
    for i in range(N):
        for j in range(N):
            # C probability species-i & species-j has interaction
            if np.random.random_sample() < C:
                # j species impact on i 
                # np.random.randn(): normally distributed random  
                A[i,j] =  round(zigma * np.random.randn(),4)
        A[i,i] = diag
    return A, r


def glv_Euler_type(x,A,r,time,Type,h1,h2):
    # x initial values: abundance
    N = len(A[:,0])
    # nt: total time points
    nt = len(time)
    dt = time[1] - time[0];
    # 1*ncol
    dxx = np.asarray(x)
    dx = np.asarray(x)
    if Type == 1:
        # each time point
        for i in range(1,nt):
            tA = np.transpose(A)
            # A'*dxx: np.matmul(tA,np.transpose(dxx))
            # np.matmul(tA,np.transpose(dxx)) + r : row-wise
            # np.multiply(*,dxx) element-wise multiplication
            # dt: scaler
            dxx = dxx + np.multiply((np.matmul(tA,np.transpose(dxx)) + r),dxx)*dt
            # matrix -> single row array
            dxx = dxx[0,:]
            dx = np.vstack((dx,dxx))
    if Type == 2:
        for i in range(1,nt):
            tA = np.transpose(A) # the (complex) conjugate transpose
            tr =  np.transpose(r)
            # np.divide(dxx,(1+h1*dxx))
            # np.matmul(tA,*) + r
            # np.multiply(*,dxx) element-wise multiplication
            # dxx + * dt
            dxx = dxx + np.multiply((np.matmul(tA,np.transpose(np.divide(dxx,(1+h1*dxx)))) + r), dxx)*dt
            dxx = dxx[0,:]
            dx = np.vstack((dx,dxx))
    if Type == 3:
        for i in range(1,nt):
            for j in range(0,N):
                #tA = A(:,j)'
                tA = np.transpose(A[:,j])
                # 1+h1*dxx+h2*dxx[j]
                a1 = 1+h1*dxx+h2*dxx[j]
                # np.divide(dxx,*)
                a2 = np.transpose(np.divide(dxx,a1))
                #np.matmul(tA,*)
                a3 = np.matmul(tA,a2) + r[0,j]
                dxx[j] = dxx[j] + a3*dxx[j]*dt
            dx = np.vstack((dx,dxx))
    if Type == 4:
        for i in range(1,nt):
            for j in range(0,N):
                tA = np.transpose(A[:,j])
                #(1+h1*dxx)*(1+h2*dxx(j,1))
                #(1+h2*dxx[j]): scalar
                a1 = (1+h1*dxx)*(1+h2*dxx[j])
                #
                a2 = np.divide(dxx, a1)
                #
                a3 = np.matmul(tA, a2) + r[0,j]
                dxx[j] = dxx[j] + a3*dxx[j]*dt
            dx = np.vstack((dx,dxx))
    dx = np.transpose(dx)
    return dx, dxx

def Generate_Network_A_of_BandW(N,C,delta,diag,VarianceType,\
                                time,FunctionType,h1,h2,Cdiff,\
                                Cdiff_disease_abundance,\
                                Cdiff_health_abundance,\
                                select_white_black_mixed):
    while 1:
        A, r = generate_A(N,C,delta,diag,VarianceType)   
        if select_white_black_mixed == 'mixed':
            #species (+,-) impact on Cdiff
            promoters  = list(np.where(A[:,Cdiff]>0)[0])
            inhibitors = list(np.where(A[:,Cdiff]<0)[0])
            neutrals = list(np.where(A[:,Cdiff]== 0)[0])
            #healthy:  bacteria + on Cdiff < bacteria - on Cdiff
            # uniformaly distributed random integer from (0, max)
            ri = random.randint(0,len(inhibitors))           
            rp = random.randint(0,math.ceil(len(promoters)*0.1)) #math.ceil: round up to the nearest integer
            rn = random.randint(0,len(neutrals))
            # k valaues sample uniformly at random w/o replacement 
            set1 = set(random.choices(inhibitors,k=ri))
            set2 = set(random.choices(promoters,k=rp))
            set3 = set(random.choices(neutrals,k=rn)) 
            # unique list of indices
            local_W_index = list(set([Cdiff])|set1|set2|set3)
            #Diseased: inhibited bacteria < promoted bacteria
            ri2 = random.randint(0,math.ceil(len(inhibitors)*0.1))
            rp2 = random.randint(0,len(promoters))
            rn2 = random.randint(0,len(neutrals))
            set1 = set(random.choices(inhibitors,k=ri2))
            set2 = set(random.choices(promoters,k=rp2))
            set3 = set(random.choices(neutrals,k=rn2)) 
            local_B_index = list(set([Cdiff])|set1|set2|set3) 
        initial = [0]*N
        for i in local_W_index:
            initial[i] += 0.2
        [dx_W,dxx_W]=glv_Euler_type(initial,A,r,time,FunctionType,h1,h2);
    
        initial = [0]*N
        for i in local_B_index:
            initial[i] += 0.2
        #dx_B: abundance at diff. timepoints
        #dxx_B: abundance at last timepoints
        [dx_B,dxx_B]=glv_Euler_type(initial,A,r,time,FunctionType,h1,h2);
        #last time poin Cdiff abundance
        if list(dxx_W)[Cdiff]<Cdiff_health_abundance \
        and list(dxx_B)[Cdiff]>Cdiff_disease_abundance:
            print("universal diseased Cdiff",list(dxx_B)[Cdiff])
            print("universal health Cdiff",list(dxx_W)[Cdiff])
            break
    return A,r 

def Generate_donor_samples(N,A,r,disordered_species,X_health_target,time,FunctionType,h1,h2,min_donor_species,max_donor_species):
    Cdiff = 0
    # A: row: species receive impact; column: species give impact
    # species impact on Cdiff
    promoters  = list(np.where(A[:,Cdiff]>0)[0])
    inhibitors = list(np.where(A[:,Cdiff]<0)[0])
    neutrals = list(np.where(A[:,Cdiff]== 0)[0])
    flag = 1
    while flag == 1:
        ri = random.randint(0,len(inhibitors))
        rp = random.randint(0,len(promoters))
        rn = random.randint(0,len(neutrals))
        set1 = set(random.choices(inhibitors,k=ri))
        set2 = set(random.choices(promoters,k=rp))
        set3 = set(random.choices(neutrals,k=rn))
        index = list(set([Cdiff])|set1|set2|set3)
        #abundance
        initial = [0]*N
        for i in index:
            initial[i] += 0.2;
        # X_donor: abundance at last time poin
        XX_donor,X_donor = glv_Euler_type(initial,A,r,time,FunctionType,h1,h2);
        donor_species_richness = sum(1 for x in X_donor if x>0)
        if X_donor[disordered_species]<X_health_target and\
        donor_species_richness>min_donor_species and \
        donor_species_richness<max_donor_species:
            flag = 0
    return XX_donor,X_donor


def Generate_disease_sample(A,r,time,FunctionType,h1,h2,min_threshold,max_threshold,Disease_threshold):
    # Cdiff index
    Cdiff = 0
    # N: number of species
    N = len(A[:,0])
    # species (+,-) impact on Cdiff
    promoters  = list(np.where(A[:,Cdiff]>0)[0])
    inhibitors = list(np.where(A[:,Cdiff]<0)[0])
    neutrals = list(np.where(A[:,Cdiff]== 0)[0])
    disordered_species = Cdiff
    flag = 1
    while flag == 1:
        count = 0
        XX_health,X_health = Generate_donor_samples(N,A,r,disordered_species,1e-4,time,FunctionType,h1,h2,50,100);
        while 1:
            # after antibiotic administration,
            # most of commensal species in the recipient's intial healthy microbiota are removed
            # including. species can inhibit the growth of Cdiff,
            ri = random.randint(0,len(inhibitors))
            rn = random.randint(0,len(neutrals))
            set1 = set(random.choices(inhibitors,k=ri))
            set3 = set(random.choices(neutrals,k=rn))
            # set of species removed (species inhibit Cdiff)
            Delete = set1.union(set3) - set([Cdiff])
            X_anti = list(X_health)
            for i in Delete:
                X_anti[i] = 0
            # XX_disease: abundance at differen time points
            # X_disease: abundance at last time points
            XX_disease,X_disease =glv_Euler_type(X_anti,A,r,time,FunctionType,h1,h2)
            Patient_species_richness = sum(1 for x in X_disease if x != 0)
            count += 1;
            if ((X_disease[disordered_species] - X_health[disordered_species]) > Disease_threshold) and\
            (Patient_species_richness>min_threshold) and Patient_species_richness<max_threshold:
                flag = 0
                break
            if count>200:
                break
    return XX_disease,X_disease,XX_health,X_health

N = 100
C = 0.4
delta = 0.2 #variance of a_ij 
diag = -1 
VarianceType = 2
time = np.arange(0,30,0.1)
h1 = 0.1
h2 = 0.1
FunctionType = 1
Cdiff_disease_abundance = 0.5;
Cdiff_health_abundance = 1e-4;
select_white_black_mixed = 'mixed';
Cdiff = 0;
Disease_threshold = Cdiff_disease_abundance;
min_donor_species = 60;
max_donor_species = 100;
min_threshold_rCDI = 10;
max_threshold_rCDI = 25;
output = "simulated_network"
for iteration in range(20,30):
    print("universal network")
    A,r = Generate_Network_A_of_BandW(N,C,delta,diag,VarianceType,\
                                      time,FunctionType,h1,h2,Cdiff,\
                                      Cdiff_disease_abundance,\
                                      Cdiff_health_abundance,\
                                      select_white_black_mixed)
    print("donor sample")
    XX_donor,X_donor = Generate_donor_samples(N,A,r,Cdiff,Cdiff_health_abundance,time,FunctionType,h1,h2,min_donor_species,max_donor_species);
    print("diseased sample")
    XX_disease,X_disease,XX_health,X_health = Generate_disease_sample(A,r,time,FunctionType,h1,h2,min_threshold_rCDI,max_threshold_rCDI,Disease_threshold);
    np.savetxt('%s/interaction_network_%i.txt'%(output,iteration),A,fmt='%.3f')
    np.savetxt('%s/growth_%i.txt'%(output,iteration),r,fmt='%.3f')
    np.savetxt('%s/donor_abd_%i.txt'%(output,iteration),X_donor,fmt='%.3f')
    np.savetxt('%s/donor_timepoints_%i.txt'%(output,iteration),XX_donor,fmt='%.3f')
    np.savetxt('%s/disease_abd_%i.txt'%(output,iteration),X_disease,fmt='%.3f')
    np.savetxt('%s/disease_timepoints_abd_%i.txt'%(output,iteration),XX_disease,fmt='%.3f')    
    np.savetxt('%s/health_abd_%i.txt'%(output,iteration),X_health,fmt='%.3f')
    np.savetxt('%s/health_timepoints_abd_%i.txt'%(output,iteration),XX_health,fmt='%.3f')
