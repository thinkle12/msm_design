from __future__ import division, print_function
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 13:56:32 2016

@author: Trent
"""


from random import choice
from subprocess import Popen,PIPE 
import numpy
import random
from sys import stdout
import pickle
import time
starttime = time.time()
from scipy.sparse import csr_matrix
from scipy.sparse import csc_matrix
import math
import time

#max_time = 252000
max_time = 50000
start_time = time.time()


AAs = ['A', 'C', 'D', 'E', 'F', 'G', 'H','I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W','Y'] #all 20 amino acids
L = 76 # the length of the protein

#overall_solution_space = []
#for i in range(L):
#    overall_solution_space.append(AAs)

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

#Temporary Text File Read in
h = open('ubiquitin_energies_State10.txt','r+')
#h = open('ubiquitin_State0_energies_500000_reg_noangpair_4positionmutations_top4AA_smallalign_State0_1200weight_2gamma_NEWTEST.txt','r+')
#h = open('ubiquitin_State0_energies_500000_reg_3angpair_4positionmutations_top4AA_smallalign_State0_1200weight_2gamma_NEWTEST.txt','r+')
#h = open('ubiquitin_State0_energies_500000_reg_6angpair_4positionmutations_top4AA_smallalign_State0_1200weight_2gamma_NEWTEST.txt','r+')
#h = open('ubiquitin_State0_energies_500000_reg_6angpair_4positionmutations_top4AA_smallalign_State0.txt','r+')
#h = open('10posmutations_ubiquitin_State0_energies_noErestriction.txt','r+')
#h = open('20posmutations_ubiquitin_State0_energies_noErestriction.txt','r+')
#h = open('35posmutations_ubiquitin_State0_energies_noErestriction.txt','r+')
#h = open('50posmutations_ubiquitin_State0_energies_noErestriction.txt','r+')
#h = open('65posmutations_ubiquitin_State0_energies_noErestriction.txt','r+')
#h = open('76posmutations_ubiquitin_State0_energies_noErestriction.txt','r+')
#h = open('ubiquitin_State0_energies_1000000_reg_6angpair_4positionmutations_top4aa_smallalign_State0_lowE7.txt','r+')
text = h.read()
text = text.split('\n')
#text = pickle.load(open('seq_energies_state0_ubiquitin.p'))
del text[-1]
sequences = [t.split(',')[0] for t in text]
Energies = [float(t.split(',')[1]) for t in text]



#hh = open('ubiquitin_State0_energies_200000_reg_30ang_smallalign_State0.txt','r+')
#textt = hh.read()
#textt = textt.split('\n')
##text = pickle.load(open('seq_energies_state0_ubiquitin.p'))
#del textt[-1]
#sequences.extend(t.split(',')[0] for t in textt)
#Energies.extend(float(t.split(',')[1]) for t in textt)
##print(len(Energies))
##print(len(sequences))
#
#
#hhh = open('ubiquitin_State0_energies_200000_reg_16ang_smallalign_State0_proper_pairs.txt','r+')
#texttt = hhh.read()
#texttt = texttt.split('\n')
##text = pickle.load(open('seq_energies_state0_ubiquitin.p'))
#del texttt[-1]
#sequences.extend(t.split(',')[0] for t in texttt)
#Energies.extend(float(t.split(',')[1]) for t in texttt)
#
#hhhh = open('ubiquitin_State0_energies_200000_01ridge_30ang_smallalign_State0.txt','r+')
#textttt = hhhh.read()
#textttt = textttt.split('\n')
##text = pickle.load(open('seq_energies_state0_ubiquitin.p'))
#del textttt[-1]
#sequences.extend(t.split(',')[0] for t in textttt)
#Energies.extend(float(t.split(',')[1]) for t in textttt)
#
#hhhhh = open('ubiquitin_State0_energies_200000_01ridge_16ang_smallalign_State0.txt','r+')
#texttttt = hhhhh.read()
#texttttt = texttttt.split('\n')
##text = pickle.load(open('seq_energies_state0_ubiquitin.p'))
#del texttttt[-1]
#sequences.extend(t.split(',')[0] for t in texttttt)
#Energies.extend(float(t.split(',')[1]) for t in texttttt)
#
#hhhhhh = open('ubiquitin_State0_energies_200000_001ridge_16ang_smallalign_State0.txt','r+')
#textttttt = hhhhhh.read()
#textttttt = textttttt.split('\n')
##text = pickle.load(open('seq_energies_state0_ubiquitin.p'))
#del textttttt[-1]
#sequences.extend(t.split(',')[0] for t in textttttt)
#Energies.extend(float(t.split(',')[1]) for t in textttttt)
#
#hhhhhhh = open('ubiquitin_State0_energies_500000_reg_6angpair_3angtriplet_smallalign_State0_proper_indicies_test.txt','r+')
#texttttttt = hhhhhhh.read()
#texttttttt = texttttttt.split('\n')
##text = pickle.load(open('seq_energies_state0_ubiquitin.p'))
#del texttttttt[-1]
#sequences.extend(t.split(',')[0] for t in texttttttt)
#Energies.extend(float(t.split(',')[1]) for t in texttttttt)
#


#hhh = open('w'+str('State0')+'.txt','r+')
#wtext = hhh.read()
#wtext = wtext.split('\n')
#del wtext[-1]
#wfloatedlist = []
#for numbers in wtext:
#    wfloatedlist.append([float(numbers)])
#w = numpy.array(wfloatedlist)


hf = open('ubiquitin_State0_energies_500000_reg_6angpair_4positionmutations_top4AA_smallalign_State0.txt','r+')
textf = hf.read()
textf = textf.split('\n')
#text = pickle.load(open('seq_energies_state0_ubiquitin.p'))
del textf[-1]
sequencesother = [t.split(',')[0] for t in textf]
Energiesother = [float(t.split(',')[1]) for t in textf]


lowEnergies = []
lowEsequences = []

for i in range(len(Energies)):
    lowEnergies.append(Energies[i])
    lowEsequences.append(sequences[i])
        
lowEtest = []
lowESeqtest = []
for i in range(50000):
    if Energiesother[i]<float(7) and len(lowEtest)<1000:
        lowEtest.append(Energiesother[i])
        lowESeqtest.append(sequencesother[i])
        
        
highEtest = []
highESeqtest = []
for i in range(50000):
    if Energiesother[i]>float(7) and len(highEtest)<1000:
        highEtest.append(Energiesother[i])
        highESeqtest.append(sequencesother[i])
#
#lowEnergies = Energies
#lowEsequences = sequences


#Start Close Pair Selection

#"Here we select the maximum distance for a pair of AAs to be considered important"
#########################################################################
closepairs = open('close_AA_pairs_ubiquitin.txt').read()
closepairs = closepairs.split(',')
closepairs.pop()

newpairs = []
for i in range(len(closepairs)):
    newpairs.append(int(closepairs[i]))

newestpairs = []
for i in range(len(newpairs)//2):
    newestpairs.append((newpairs[2*i],newpairs[2*i+1]))

position_pairs = newestpairs

#End of Close Pair Selection

#Start Close Triplet Selection

#"Here we select the maximum distance for a triplet of AAs to be considered important"
#########################################################################
newtrips = []
closetriplets = open('close_AA_triplets_ubiquitin.txt').read()
closetriplets = closetriplets.split('\n')
del closetriplets[-1]
for ctriplets in closetriplets:
    newtrips.append(ctriplets.split(','))
    
floatedtrips = []
for trips in newtrips:
    del trips[-1]
    singletrip = []
    for f in range(len(trips)):
        singletrip.append(int(trips[f]))
    floatedtrips.append(singletrip)
        
totaltriplets = [tuple(l) for l in floatedtrips]
        
#End of Close Triplet Selection


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
    
def triplets(aa1,aa2,aa3,i,h,r):
    AAchanges1 = [aa for aa in overall_solution_space[i] if aa!=reference_seq[i]]
    AAchanges2 = [bb for bb in overall_solution_space[h] if bb!=reference_seq[h]]
    AAchanges3 = [cc for cc in overall_solution_space[r] if cc!=reference_seq[r]]
    Phi = []
    for aa in AAchanges1:
        for bb in AAchanges2:
            for cc in AAchanges3:
                if (aa1==aa and aa2==bb and aa3==cc):
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
#lambdaa = 10
lambdaa = 'RegularSGD.'
trainingsets = len(lowEnergies)
state = 'state0'

E_list = []
random_Seq_list = []
random_Seq_testinglist = []
CC = []
MAD = []
#w = []
Philist = []
seq_list = []
trainE = []
trainPhi = []
sumterm = []
print('Stochastic Gradient Descent Multivariate Linear Regression')
print('# of pair interactions = %s, %s training sequences, gamma = 2/((i+1)**0.5), lambda = %s, w start[(0,.001)]' % (len(position_pairs),trainingsets, lambdaa))
print('/n')

print('4posmutations, 3 gamma 1200 weight 6 pairs,7E Threshold, lowE defined as less than 7')
print('textfile = ubiquitin_State0_energies_500000_reg_6angpair_4positionmutations_top4AA_smallalign_State0_1200weight_2gamma_NEWTEST.txt')


for n in range(trainingsets):
    if time.time()-start_time<max_time:
    #    random_seq = ''.join([choice(overall_solution_space[i]) for i in range(L)]) # same a random sequence
    #    cmd = 'python score_sequence_rosetta.py %s' % random_seq #1. generate a terminal command as a string
    #    output = Popen(cmd,shell=True,stdout=PIPE).communicate() #2. send the command to the terminal shell
    #    seq_E = float(output[0]) #3. read the result and convert to a float
    #    
    #    print(random_seq, seq_E)
    #    open(seqEfilename,'a').write(random_seq+','+str(seq_E)+'\n')
        random_Seq_list.append(lowEsequences[n])
        E_list.append(lowEnergies[n])
        #maxishvalue = sorted(E_list)[-int(len(E_list)*.01)]
        #median_e_value = numpy.median(E_list)
        #divisor_value = (median_e_value-maxishvalue)/math.log(.0189)
        #Convert Seq to Phi Here
        seq_Phi = [1] 
        for i in range(L):
            Phi_i = sigma2Phi(random_Seq_list[n][i],i)
            seq_Phi.extend(Phi_i)
    #    for pair in position_pairs:
    #        firstpos = pair[0]
    #        secondpos = pair[1]
    #        pairs_j = pairs(random_Seq_list[n][firstpos],random_Seq_list[n][secondpos],firstpos,secondpos)
    #        seq_Phi.extend(pairs_j)
    #    for triplet in totaltriplets:
    #        tfirstpos = triplet[0]
    #        tsecondpos = triplet[1]
    #        tthirdpos = triplet[2]
    #        triplets_j = triplets(lowEsequences[tfirstpos],lowEsequences[tsecondpos],lowEsequences[tthirdpos],tfirstpos,tsecondpos,tthirdpos)
    #        seq_Phi.extend(triplets_j)
        if n<100:
            Philist.append(seq_Phi)
            if n==0:
                d = len(seq_Phi)
                w = numpy.zeros(shape=(d,1)) #Starting vector for regular regression not L1
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
            #weight = float(min(1,math.exp((7-trainE)/float(1200))))
            #weight = float(min(1,math.exp((median_e_value-trainE)/float(divisor_value))))
            weight = 1
            
            gamma = 2/((n-100+1)**0.5)
    #        gamma = 6/((n-100+1)**0.5)
    #        #gamma = 1/(n+1)
    #        if n<50000:
    #            gamma = 1/((n-100+1)**0.5) # have to do plus one cause first round i = 0
    #        if n<100000:
    #            gamma = 2/((n-100+1)**0.5)
    #        if n<150000:
    #            gamma = 3/((n-100+1)**0.5)
    #        if n<200000:
    #            gamma = 4/((n-100+1)**0.5)
            #if n%100==0 and MAD[-1]/MAD[-100]<1:
                
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
            CC.append(cc) # calculate the correlation coeffcient, append into list
            MAD.append(float(mad))# calculate the mean absolute deviation, append into list
            stdout.write('\rTraining set = %s' % str(n+1)) #This prints the current training set index on same line as an update
            stdout.flush()
            madlist100 = []
            if len(MAD)>=100:
                for g in range(1,101):
                    madlist100.append(MAD[-g])
                if all(mm<float(2) for mm in madlist100):
                    break

#open('w'+str(pdbfile)+'.txt','w')
#open('wState0.txt','w')
#for vv in range(len(w)):
#    open('wState0.txt','a').write(str(w[vv][0])+'\n')
#        
import matplotlib.pyplot as plt
plt.plot(range(len(MAD)),MAD)
plt.show()


print('energy function fit:')
print('\tCC = %0.2f' % CC[(len(MAD))-1])
print('\tMAD = %0.2f' % MAD[(len(MAD))-1])

plt.scatter(E_testinglist,E_estimate)
plt.show()

#removepointlambdaa = str('83ktrainedregSGD')
#if removepointlambdaa[0]=='0':
#    removepointlambdaa=removepointlambdaa[1:]
#path="C:\Users\Trent\Desktop\PythonFiles\Project\PickleFilesGraphComparison\\"
#import pickle
#pickle.dump(w,open(path+'L1_'+str(removepointlambdaa)+"_w.p", "wb"))
#pickle.dump(E_estimate,open(path+'L1_'+str(removepointlambdaa)+"_E_estimate.p", "wb"))
#pickle.dump(testE,open(path+'L1_'+str(removepointlambdaa)+"_testE.p", "wb"))
#pickle.dump(CC,open(path+'L1_'+str(removepointlambdaa)+"_CC.p", "wb"))
#pickle.dump(MAD,open(path+'L1_'+str(removepointlambdaa)+"_MAD.p", "wb"))

#Test to see how many features are zero...
zerolist = []
for z in range(len(w)):
    if w[z][0]==0:
        zerolist.append(z)
print('Number of Features Equal to Zero in the Model = %s' % len(zerolist))

elapsedtime = time.time() - starttime 
print('time = %0.2f' % elapsedtime)

sumlist = []
for i in range(1000):
    sumlist.append(MAD[-i])
print('Avg MAD of last 1000 terms = '+str(sum(sumlist)/1000))
            
lowEPhilist = []
for n in range(len(lowESeqtest)):
    seq_Phi = [1] 
    for i in range(L):
        Phi_i = sigma2Phi(lowESeqtest[n][i],i)
        seq_Phi.extend(Phi_i)
#    for pair in position_pairs:
#        firstpos = pair[0]
#        secondpos = pair[1]
#        pairs_j = pairs(lowESeqtest[n][firstpos],lowESeqtest[n][secondpos],firstpos,secondpos)
#        seq_Phi.extend(pairs_j)
    lowEPhilist.append(seq_Phi)
    
lowE_estimatelist = []
lowEmadlist = []
for s in range(len(lowEPhilist)):
    lowE_estimate=numpy.dot(lowEPhilist[s],w)
    lowE_estimatelist.append(lowE_estimate)
    lowEmad=numpy.mean(abs(lowE_estimatelist[s]-lowEtest[s]))
    lowEmadlist.append(lowEmad)

mad1000lowE = float(sum(lowEmadlist))/len(lowEtest)
    
highEPhilist = []
for n in range(len(highESeqtest)):
    seq_Phi = [1] 
    for i in range(L):
        Phi_i = sigma2Phi(highESeqtest[n][i],i)
        seq_Phi.extend(Phi_i)
#    for pair in position_pairs:
#        firstpos = pair[0]
#        secondpos = pair[1]
#        pairs_j = pairs(lowESeqtest[n][firstpos],lowESeqtest[n][secondpos],firstpos,secondpos)
#        seq_Phi.extend(pairs_j)
    highEPhilist.append(seq_Phi)
    
highE_estimatelist = []
highEmadlist = []
for s in range(len(highEPhilist)):
    highE_estimate=numpy.dot(highEPhilist[s],w)
    highE_estimatelist.append(highE_estimate)
    highEmad=numpy.mean(abs(highE_estimatelist[s]-highEtest[s]))
    highEmadlist.append(highEmad)
    
mad1000highE = float(sum(highEmadlist))/len(highEtest)

under100above7Etest = []
under100above7ESeqtest = []
for i in range(50000):
    if Energiesother[i]>float(7) and Energiesother[i]<100 and len(under100above7Etest)<1000:
        under100above7Etest.append(Energiesother[i])
        under100above7ESeqtest.append(sequencesother[i])
        
under100above7EPhilist = []
for n in range(len(under100above7ESeqtest)):
    seq_Phi = [1] 
    for i in range(L):
        Phi_i = sigma2Phi(under100above7ESeqtest[n][i],i)
        seq_Phi.extend(Phi_i)
#    for pair in position_pairs:
#        firstpos = pair[0]
#        secondpos = pair[1]
#        pairs_j = pairs(under100above7ESeqtest[n][firstpos],under100above7ESeqtest[n][secondpos],firstpos,secondpos)
#        seq_Phi.extend(pairs_j)
    under100above7EPhilist.append(seq_Phi)
    
under100above7E_estimatelist = []
under100above7Emadlist = []
for s in range(len(under100above7EPhilist)):
    under100above7E_estimate=numpy.dot(under100above7EPhilist[s],w)
    under100above7E_estimatelist.append(under100above7E_estimate)
    under100above7Emad=numpy.mean(abs(under100above7E_estimatelist[s]-under100above7Etest[s]))
    under100above7Emadlist.append(under100above7Emad)
    
mad1000under100above7E = float(sum(under100above7Emadlist))/len(under100above7Etest)

print('Low E MAD = '+str(mad1000lowE))
print('High E MAD = '+str(mad1000highE))
print('7<E<100 MAD = '+str(mad1000under100above7E))
