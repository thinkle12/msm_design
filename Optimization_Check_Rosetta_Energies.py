# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 02:39:38 2017

@author: Trent
"""

state_number = '97'
chosen_state_number = 97
chosen_state = 'State97'


import re
import math
import numpy
import os


hh = open('one_sequence_100structures_energies_4mutations_new.txt','r+')
textt = hh.read()
textt = textt.split('\n')
del textt[-1]
newtextt = [tt.split(',') for tt in textt]
numlist = []
for states in newtextt:
    numlist.append(re.sub("[^0-9]",'',states[1]))
for j in range(len(newtextt)):
    newtextt[j][1]=numlist[j]
for j in range(len(newtextt)):
    newtextt[j][1]=int(newtextt[j][1])
sortedtextt = sorted(newtextt, key=lambda x: x[1])
Energies = [float(sortedtextt[t][0]) for t in range(len(sortedtextt))]
States = [str(sortedtextt[t][1]) for t in range(len(sortedtextt))]
 

k=1
t=1


    


#This is Energies[-3] because its the one we are optimizing...
def rosetta_energy_single(select_state):
        E_estimate = math.exp(-(Energies[States.index(select_state)])/(k*t))
        return E_estimate
        
def rosetta_energy_multi():
    State_E_list = []
    for a in range(len(Energies)):
        E_estimate_mult = math.exp((-Energies[a])/(k*t))
        State_E_list.append(E_estimate_mult)
    return sum(State_E_list)
    
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




    
def single_energy_estimation(binary, State):
        E_estimate = math.exp(-(numpy.dot(binary,wStates['w'+str(chosen_state)]))/(k*t))
        return E_estimate
        
def multi_energy_estimation(binary):
    State_E_list = []
    for a in range(len(namelist)):
        E_estimate_mult = math.exp(-(numpy.dot(binary,wStates['wState'+str(namelist[a])]))/(k*t))
        State_E_list.append(E_estimate_mult)
    return sum(State_E_list)
  

for s in range(len(namelist)):
    namelist[s]=int(namelist[s])
#
state_e_rosetta_list = []
updated_States = []
updated_Energies = []
for kk in range(len(States)):
    if int(States[kk]) in namelist:
        ind = States.index(States[kk])
        updated_States.append(States[ind])
        updated_Energies.append(Energies[ind])
        state_e_rosetta_list.append([Energies[ind],int(States[ind])])
    



#opt_seq = 'LQIFVKTPTGKTITIEVEPSDTIENVRAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
#opt_seq = list('LQIFVKTPTGKTITIEVEPSDTIENVRAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG') 
#opt_seq = 'MQIFVKTPTGKTITIEVEPSDTIENVKAKIQDKEGIPLDQQHLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
#opt_seq = list('MQIFVKTPTGKTITIEVEPSDTIENVKAKIQDKEGIPLDQQHLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG')
#opt_seq = 'LQIFVKTPTGKTITIEVEPSDTIENVKAKIQDKEGIPPDQQRLIIAGKQLEDGRTLSDYNIQKESTLHLVLRLRGS'
#opt_seq = list('LQIFVKTPTGKTITIEVEPSDTIENVKAKIQDKEGIPPDQQRLIIAGKQLEDGRTLSDYNIQKESTLHLVLRLRGS')
#opt_seq = 'LQIFLETPTGKTITLEVEPSDTIENVKAKIQDKVGIPLDEQRLIFAGKQLADGRTVSDYNIQKESTLYLVLRLRGG'
#opt_seq = list('LQIFLETPTGKTITLEVEPSDTIENVKAKIQDKVGIPLDEQRLIFAGKQLADGRTVSDYNIQKESTLYLVLRLRGG')
#opt_seq = 'LQIFAEMPTEKTIEVEVEPSDTLANMKAKIHDKVGLPPVQQRLIFSGKQLEDGRTLSDYNIQKESTLHLVFRLCGG'
#opt_seq = list('LQIFAEMPTEKTIEVEVEPSDTLANMKAKIHDKVGLPPVQQRLIFSGKQLEDGRTLSDYNIQKESTLHLVFRLCGG')
opt_seq = 'MQIFVKTLTEKTITVEVEPSDTIENVKAKIQDKVGIPPDQQRLIIAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
opt_seq = list('MQIFVKTLTEKTITVEVEPSDTIENVKAKIQDKVGIPPDQQRLIIAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG')

bin_seq = binary_conversion(opt_seq)
keylist = wStates.keys()
state_e_est_list = []
for stringname in keylist:
    state_e_est_list.append([numpy.dot(bin_seq,wStates[stringname]),stringname])

for l in range(len(state_e_est_list)):
    state_e_est_list[l][1]=int(re.sub("[^0-9]",'',state_e_est_list[l][1]))

sorted_state_e_est_list = sorted(state_e_est_list, key=lambda x: x[1])
sorted_state_e_rosetta_list = sorted(state_e_rosetta_list, key=lambda x: x[1])

single_E = rosetta_energy_single(state_number)
multi_E = rosetta_energy_multi()
boltz_dist = single_E/multi_E
print 'Probability = '+str(boltz_dist)


joined_opt_seq = ''.join(opt_seq)

IDseq = []
for jj in range(L):
    if joined_opt_seq[jj]!=wildtype[jj]:
        IDseq.append(joined_opt_seq[jj].lower())
    if joined_opt_seq[jj]==wildtype[jj]:
        IDseq.append(joined_opt_seq[jj])
        

import matplotlib.pyplot as plt

plt.title('Rosetta E vs Model E for 100 conformations'+'\n'+str(''.join(IDseq))+'\n')
plt.scatter([sorted_state_e_est_list[i][0] for i in range(len(sorted_state_e_est_list))],[sorted_state_e_rosetta_list[i][0] for i in range(len(sorted_state_e_rosetta_list))])
plt.plot([-50,2000],[-50,2000])
plt.ylabel('Rosetta Energies')
plt.xlabel('SGD Model Energies')
plt.ylim([-50,700])
plt.xlim([-50,700])
plt.show()

