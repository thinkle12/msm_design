from __future__ import division, print_function
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 13:56:32 2016

@author: Trent
"""

import os
import math
from random import choice
from subprocess import Popen,PIPE 
import numpy
import random
from sys import stdout
import pickle
import time
import itertools
starttime = time.time()


import time

max_time = 252000





#File Imports and Identifiers
#Identify by #seq/E, type of learning, pair distance, alignment, and pdb state

namelist = []
for i in range(100):
    namelist.append('State'+str(i)+'.pdb')

for filename in os.listdir(os.curdir):
    if filename in namelist:
        pdbfile = filename

#print(pdbfile)
pdbtag = pdbfile[:-4]
#print(pdbtag)



hhh = open('w'+str(pdbtag)+'.txt','r+')
wtext = hhh.read()
wtext = wtext.split('\n')
del wtext[-1]
wfloatedlist = []
for numbers in wtext:
    wfloatedlist.append(float(numbers))
w = numpy.array(wfloatedlist)



#pdbfile = 'State0.pdb'
alignmentfile = 'smallalignment.txt'
pairfile = 'close_AA_pairs_ubiquitin.txt'

#Regression Parameters
#lambdaa = .0001
#ridgecoef = .01
lambdaa = 'RegularSGD.'
trainingsets = 50000
madcutoff = 2.0


#Identifier = str(trainingsets)+'_reg_6angpair_4positionmutations_top4AA_smallalign_State0_1200weight_2gamma_NEWTEST'
#
seqEfilename = 'ubiquitin_energies_'+str(pdbtag)+'.txt'
#
#pickle_dump_dir = '/home/thinkle2@ad.wisc.edu/pickle_output_files//'

AAs = ['A', 'C', 'D', 'E', 'F', 'G', 'H','I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W','Y'] #all 20 amino acids
L = 76 # the length of the protein

wildtype = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
wtlist = [wildtype[i] for i in range(len(wildtype))]

#################################################################################
#################################################################################
#################################################################################

h = open(alignmentfile,'r+')
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
for top4 in overall_solution_space:
    del top4[-1]
#for top3 in overall_solution_space:
#    del top3[-1]

print(overall_solution_space)
#reference_seq = [choice(overall_solution_space[i]) for i in range(len(overall_solution_space))]
reference_seq = wtlist
print(reference_seq)

#Close Pair Selection
#Here we read in the pairs we want to use from running the file coordinate_pair_selection.py
closepairs = open(pairfile).read()
closepairs = closepairs.split(',')
closepairs.pop()

newpairs = []
for i in range(len(closepairs)):
    newpairs.append(int(closepairs[i]))

newestpairs = []
for i in range(len(newpairs)//2):
    newestpairs.append((newpairs[2*i],newpairs[2*i+1]))

position_pairs = newestpairs

#Definitions
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

def pairs(aa1,aa2,i,h):
    AAchanges1 = [aa for aa in overall_solution_space[i] if aa!=reference_seq[i]]
    AAchanges2 = [bb for bb in overall_solution_space[h] if bb!=reference_seq[h]]
    Phi = []
    for aa in AAchanges1:
        for bb in AAchanges2:
            if (aa1==aa and aa2==bb):
                Phi.append(1)
            else:
                Phi.append(0)
    return Phi
    
def STF(a,z):
    if a>z:
        return a-z
    if a<-z:
        return a+z
    if -z<=a<=z:
        return 0

    

# Sample a random sequence, evaluate its energy according to our energy fucntion, and store the sequence and energy data


E_list = []
random_Seq_list = []
random_Seq_testinglist = []
CC = []
MAD = []
w = []
Philist = []
print('Stochastic Gradient Descent')
print('# of pair interactions = %s, %s training sequences, gamma = 2/((i+1)**0.5), lambda = %s, w start=zeros' % (len(position_pairs),trainingsets, lambdaa))


open(seqEfilename,'w')
open('CC'+str(pdbtag)+'.txt','w')
open('MAD'+str(pdbtag)+'.txt','w')


for n in range(trainingsets):
    if time.time()-starttime<max_time:
        #random_seq = ''.join([choice(overall_solution_space[i]) for i in range(L)]) # same a random sequence
        wildtype = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
        wtlist = [wildtype[i] for i in range(len(wildtype))]
        randAAlist = []
        for i in range(4):
            randAAlist.append(random.randint(0,75))
        for mutantpositions in randAAlist:
            wtlist[mutantpositions] = choice(overall_solution_space[mutantpositions])
        random_seq = ''.join(wtlist)
        cmd = 'python score_sequence_rosetta_submitnode.py %s' % random_seq #1. generate a terminal command as a string
        output = Popen(cmd,shell=True,stdout=PIPE).communicate() #2. send the command to the terminal shell
        seq_E = float(output[0]) #3. read the result and convert to a float
        
        print(random_seq, seq_E)
        open(seqEfilename,'a').write(random_seq+','+str(seq_E)+'\n')
        #random_Seq_list.append(sequences[n])
        E_list.append(seq_E)
        maxishvalue = sorted(E_list)[-int(len(E_list)*.01)]
        median_e_value = numpy.median(E_list)
        divisor_value = (median_e_value-maxishvalue)/math.log(.0189)
        #Convert Seq to Phi Here
        seq_Phi = [1] 
        for i in range(L):
            Phi_i = sigma2Phi(random_seq[i],i)
            seq_Phi.extend(Phi_i)
        #    for pair in position_pairs:
        #        firstpos = pair[0]
        #        secondpos = pair[1]
        #        pairs_j = pairs(random_seq[firstpos],random_seq[secondpos],firstpos,secondpos)
        #        seq_Phi.extend(pairs_j)
        if n<100:
            Philist.append(seq_Phi)
            if n==0:
                d = len(seq_Phi)
                #w = numpy.zeros(shape=(d,1)) #Starting vector for regular regression not L1
                #for f in range(1,len(Phi_i*L)):
                #    w[f][0] = .05
                print('Number of Variables in Model = '+str(d))
        if n>=100:
            Philist.append(seq_Phi)
            trainPhi = Philist[0]
            del Philist[0]
            E_testinglist = [E_list[-xx] for xx in range(1,101)]
            E_testinglist.reverse()
            trainE = E_list[n-100]
            weight = float(min(1,math.exp((median_e_value-trainE)/float(divisor_value))))
            
            gamma = 2/((n+1)**0.5) # have to do plus one cause first round i = 0
            x = numpy.array([trainPhi]).T # column vector
            #Regular
            w = w-gamma*x*weight*(numpy.dot(x.T,w)*weight - trainE*weight)
            #Ridge
            #w = w-gamma*(x*(numpy.dot(x.T,w) - trainE) + ridgecoef*w)
            #Lasso
            #for j in range(d):
            #    w[j][0]=STF(w[j][0],gamma*lambdaa)
        
            # analyze the fit of the current model w
            #Regular Matrix Multiplication
            #E_estimate = numpy.dot(newTestPhi,w)[:,0]
            #Iterative Matrix Multiplication
            #if not k%100:
            E_estimate=numpy.empty([len(Philist),1])
            for s in range(len(Philist)):
                E_estimate[s]=numpy.dot(Philist[s],w[:,0])
            E_estimate = E_estimate[:,0]
            cc=numpy.corrcoef(E_testinglist,E_estimate)[0,1]
            mad=numpy.mean(abs(E_estimate-E_testinglist))
            print('cc = '+str(cc))
            print('mad = '+str(mad))
            CC.append(cc) # calculate the correlation coeffcient, append into list
            MAD.append(mad)# calculate the mean absolute deviation, append into list
            open('CC'+str(pdbtag)+'.txt','a').write(str(cc)+'\n')
            open('MAD'+str(pdbtag)+'.txt','a').write(str(mad)+'\n')
            madlist100 = []
            if len(MAD)>=100:
                for g in range(1,101):
                    madlist100.append(MAD[-g])
                if all(mm<float(madcutoff) for mm in madlist100):
                    break

open('w'+str(pdbtag)+'.txt','w')
#open('wState0.txt','w')
for vv in range(len(w)):
    open('w'+str(pdbtag)+'.txt','a').write(str(w[vv][0])+'\n')


print('energy function fit:')
print('\tCC = %0.2f' % CC[(len(MAD))-1])
print('\tMAD = %0.2f' % MAD[(len(MAD))-1])

#plt.scatter(E_testinglist,E_estimate)
#plt.show()




#Test to see how many features are zero...
zerolist = []
for z in range(len(w)):
    if w[z][0]==0:
        zerolist.append(z)
print('Number of Features Equal to Zero in the Model = %s' % len(zerolist))

elapsedtime = time.time() - starttime 
print('time = %0.2f' % elapsedtime)

