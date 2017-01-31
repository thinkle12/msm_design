# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 23:26:11 2017

@author: Trent
"""


L=76


#k=9.83e-25
k=1
t=298.15
chosen_state_number = 97
chosen_state = 'State97'
mutation_threshold = L-6


import os
import re
import numpy
import random
import math

filelist = os.listdir("bad_mad")

#print filelist

namelist = []
for names in filelist:
    namelist.append(re.sub("[^0-9]",'',names))

MADStates = {}
for q in range(len(filelist)):
    MADStates['MADState'+str(namelist[q])] = []
    
    
sumlist = []
for j in range(len(filelist)):
    h = open('C:\Users\Trent\Desktop\PythonFiles\Project\\bad_mad\\'+str(filelist[j]),'r+')
    text = h.read()
    text = text.split('\n')
    #text = pickle.load(open('seq_energies_state0_ubiquitin.p'))
    del text[-1]
    floatedlist = []
    for items in text:
        floatedlist.append(float(items))
    if floatedlist != []:
        MADStates['MADState'+str(namelist[j])] = floatedlist
        thousandlist = []
        for i in range(1,1001):
            thousandlist.append(floatedlist[-i])
        sumlist.append(sum(thousandlist)/1000)

newfilelist = os.listdir("bad_wfiles")

newnamelist = []
for newnames in newfilelist:
    newnamelist.append(re.sub("[^0-9]",'',newnames))

wStates = {}
for q in range(len(filelist)):
    wStates['wState'+str(newnamelist[q])] = []


for j in range(len(filelist)):
    hh = open('C:\Users\Trent\Desktop\PythonFiles\Project\\bad_wfiles\\'+str(newfilelist[j]),'r+')
    newtext = hh.read()
    newtext = newtext.split('\n')
    del newtext[-1]
    newfloatedlist = []
    for numbers in newtext:
        newfloatedlist.append(float(numbers))
    wStates['wState'+str(newnamelist[j])] = numpy.array(newfloatedlist)
    

AAs = ['A', 'C', 'D', 'E', 'F', 'G', 'H','I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W','Y'] #all 20 amino acids
L = 76 # the length of the protein

# this is the amino acids allowed at each postion
h = open('smallalignment.txt','r+')
text = h.read()
text = text.split('\n')
del text[-1]

overall_solution_space = []
for topaminoacids in text:
    tempsolspace = []
    for i in range(len(topaminoacids)):
        tempsolspace.append(topaminoacids[i])
    overall_solution_space.append(tempsolspace)
    
for top5 in overall_solution_space:
    del top5[-1]    
print(overall_solution_space)

wildtype = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
wtlist = [wildtype[i] for i in range(len(wildtype))]
reference_seq = wtlist
print(reference_seq)

    
def sigma2Phi(sigma,i):
    """This takes in an amino acid and a position and returns the 1x19 binary indicator vector Phi"""
    AAchanges = [aa for aa in overall_solution_space[i] if aa!=reference_seq[i]] # these are all 19 possible AA changes (i.e. AAs that are not the reference sequence AA)
    Phi = []
    for aa in AAchanges:
        if sigma==aa:
            Phi.append(1)
        else:
            Phi.append(0)
    return Phi
    



def binary_conversion(sequence):
    seq_Phi = [1] 
    for i in range(L):
        Phi_i = sigma2Phi(sequence[i],i)
        seq_Phi.extend(Phi_i)
    #    for pair in position_pairs:
    #        firstpos = pair[0]
    #        secondpos = pair[1]
    #        pairs_j = pairs(lowESeqtest[n][firstpos],lowESeqtest[n][secondpos],firstpos,secondpos)
    #        seq_Phi.extend(pairs_j)
    return seq_Phi







def single_energy_estimation(binary, State_to_optimize):
        E_estimate = math.exp(-(.57*numpy.dot(binary,wStates['w'+str(State_to_optimize)]))/(k*t))
        return E_estimate
        
def multi_energy_estimation(binary):
    State_E_list = []
    for a in range(len(namelist)):
        E_estimate_mult = math.exp(-(.57*numpy.dot(binary,wStates['wState'+str(namelist[a])]))/(k*t))
        State_E_list.append(E_estimate_mult)
    return sum(State_E_list)




def percentile_function(x):
    return math.exp(-float(x)/2500)



#Might need to take out WT mutations from overall_solution_space to avoid mutating to wt
rand_seq_list = []
boltz_list = []
#Test Optimization 
for n in range(100000):
    if n==0:
        aminoacidchain = [wildtype[i] for i in range(len(wildtype))]
        random_seq = ''.join(aminoacidchain)
        rand_seq_list.append(aminoacidchain)
        Phi = binary_conversion(random_seq)
        single_E = single_energy_estimation(Phi,chosen_state)
        multi_E = multi_energy_estimation(Phi)
        boltz_dist = single_E/multi_E
        boltz_list.append(boltz_dist)
    if n>0:
        randAAlist = []
        for i in range(1):
            randAAlist.append(random.randint(0,75))
        if len([a for a, b in zip(rand_seq_list[-1], reference_seq) if a == b])<=mutation_threshold:
            aachaintobemutated = list(aminoacidchain)
        if len([a for a, b in zip(rand_seq_list[-1], reference_seq) if a == b])>mutation_threshold:
            aachaintobemutated = list(rand_seq_list[-1])
        for mutantpositions in randAAlist:
            aachaintobemutated[mutantpositions] = random.choice(overall_solution_space[mutantpositions])
        random_seq = ''.join(aachaintobemutated)
        listseq = list(random_seq)
        threshold = percentile_function(n)
        if len([a for a, b in zip(rand_seq_list[-1], reference_seq) if a == b])<=mutation_threshold:
            rand_seq_list.append(aminoacidchain)
            Phi = binary_conversion(rand_seq_list[-1])
            single_E = single_energy_estimation(Phi,chosen_state)
            multi_E = multi_energy_estimation(Phi)
            boltz_dist = single_E/multi_E
            boltz_list.append(boltz_dist)
        if len([a for a, b in zip(rand_seq_list[-1], reference_seq) if a == b])>mutation_threshold:
            rand_seq_list.append(listseq)
            Phi = binary_conversion(rand_seq_list[-1])
            single_E = single_energy_estimation(Phi,chosen_state)
            multi_E = multi_energy_estimation(Phi)
            boltz_dist = single_E/multi_E
            boltz_list.append(boltz_dist)
            if boltz_list[-1]>boltz_list[-2]:
                pass
            if boltz_list[-1]<=boltz_list[-2]:
                del rand_seq_list[-1]
                del boltz_list[-1]
#            probability = random.random()
#            if probability<threshold:
#                if boltz_list[-1]>=boltz_list[-2]:
#                    pass
#                if boltz_list[-1]<boltz_list[-2]:
#                    pass
#            if probability>=threshold:
#                if boltz_list[-1]>boltz_list[-2]:
#                    pass
#                if boltz_list[-1]<=boltz_list[-2]:
#                    del rand_seq_list[-1]
#                    del boltz_list[-1]
                
        
#print threshold
    
maxindex = boltz_list.index(max(boltz_list))
candidate_sequence = ''.join(rand_seq_list[maxindex])
max_boltz = boltz_list[maxindex]

print candidate_sequence
print max_boltz




    