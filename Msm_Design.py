from __future__ import division, print_function
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 13:56:32 2016

@author: Trent
"""

class SGD_Static:
    
    def __init__(self,pdbfile,sequence_energy_file,sequence_alignment_file,wt_seq,pair_select_list=False,gamma_multiplier=2,regular_linear_regression=True,lasso=False,lambda_lasso_coef=.01,ridge=False,ridge_coef=.01):
        self.pdbfile = pdbfile
        self.sequence_energy_file = sequence_energy_file
        self.sequence_alignment_file = sequence_alignment_file
        self.pair_select_list = pair_select_list
        self.gamma_multiplier = gamma_multiplier
        self.wt_seq = wt_seq
        self.regular_linear_regression = regular_linear_regression
        self.lasso = lasso
        self.ridge = ridge
        self.lambda_lasso_coef = lambda_lasso_coef
        self.ridge_coef = ridge_coef
	

    def static_model(self,):
        if not self.regular_linear_regression and not self.lasso and not self.ridge:
            raise ValueError('No model type selected, one of regular_linear_regression, lasso, or ridge must be true')
        if self.regular_linear_regression and self.lasso:
            raise ValueError('Cannot use both regular and lasso regression, select one or the other')
        if self.regular_linear_regression and self.ridge:
            raise ValueError('Cannot use both regular and ridge regression, select one or the other')
        if self.ridge and self.lasso:
            raise ValueError('Cannot use both ridge and lasso regression, select one or the other')
        if self.regular_linear_regression and self.lasso and self.ridge:
            raise ValueError('Cannot use regular, ridge, and lasso regression at the same time, select one or the other')
        
        
        import numpy
        from sys import stdout
        import time
        starttime = time.time()

        import time
        
        #max_time = 252000
        max_time = 5000000
        start_time = time.time()
        
        
        AAs = ['A', 'C', 'D', 'E', 'F', 'G', 'H','I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W','Y'] #all 20 amino acids
        L = len(self.wt_seq) # the length of the protein
        
        
        
        # this is the amino acids allowed at each postion
        h = open(self.sequence_alignment_file,'r+')
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
        
        
        #wildtype = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
        wildtype = self.wt_seq
        wtlist = [wildtype[i] for i in range(len(wildtype))]
        reference_seq = wtlist
        print(reference_seq)
        
        
        h = open(self.sequence_energy_file,'r+')
        text = h.read()
        text = text.split('\n')
        del text[-1]
        sequences = [t.split(',')[0] for t in text]
        Energies = [float(t.split(',')[1]) for t in text]
        
        

        if self.pair_select_list:
            position_pairs = self.pair_select_list
        else:
            position_pairs = []
        
        #Definitions
        def sigma2Phi(sigma,i):
            """This takes in an amino acid and a position and returns the 1x19 binary indicator vector Phi"""
            AAchanges = [aa for aa in overall_solution_space[i] if aa!=reference_seq[i]] # these are all 19 possible AA changes (i.e. AAs that are not the reference sequence AA)
            #AAchanges = [aa for aa in AAs if aa!=reference_seq[i]]
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
        #lambdaa = .01
        trainingsets = len(Energies)
        print('Number of training sets = '+str(len(Energies)))
        
        #ridgecoef = .000001
        E_list = []
        random_Seq_list = []
        CC = []
        MAD = []
        #w = []
        Philist = []
        trainE = []
        trainPhi = []
        print('Stochastic Gradient Descent Multivariate Linear Regression')
        print('# of pair interactions = %s, gamma = %s/((i+1)**0.5)' % (len(position_pairs), self.gamma_multiplier))
        print('\n')
        
        
        
        for n in range(trainingsets):
            if time.time()-start_time<max_time:
            #    random_seq = ''.join([choice(overall_solution_space[i]) for i in range(L)]) # same a random sequence
            #    cmd = 'python score_sequence_rosetta.py %s' % random_seq #1. generate a terminal command as a string
            #    output = Popen(cmd,shell=True,stdout=PIPE).communicate() #2. send the command to the terminal shell
            #    seq_E = float(output[0]) #3. read the result and convert to a float
            #    
            #    print(random_seq, seq_E)
            #    open(seqEfilename,'a').write(random_seq+','+str(seq_E)+'\n')
                random_Seq_list.append(sequences[n])
                E_list.append(Energies[n])

                #Convert Seq to Phi Here
                seq_Phi = [1] 
                for i in range(L):
                    Phi_i = sigma2Phi(random_Seq_list[n][i],i)
                    seq_Phi.extend(Phi_i)
                if self.pair_select_list:
                    for pair in position_pairs:
                        firstpos = pair[0]
                        secondpos = pair[1]
                        pairs_j = pairs(random_Seq_list[n][firstpos],random_Seq_list[n][secondpos],firstpos,secondpos)
                        seq_Phi.extend(pairs_j)
                else:
                    pass
                if n<100:
                    Philist.append(seq_Phi)
                    if n==0:
                        d = len(seq_Phi)
                        w = numpy.zeros(shape=(d,1)) #Starting vector for regular regression not L1
                        if self.lasso:
                            for f in range(1,len(Phi_i*L)):
                                w[f][0] = .01
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
                    
                    gamma = self.gamma_multiplier/((n-100+1)**0.5)

                        
                    x = numpy.array([trainPhi]).T # column vector
                    if self.regular_linear_regression:
                        w = w-gamma*x*weight*(numpy.dot(x.T,w)*weight - trainE*weight)
                    #Ridge
                    if self.ridge:
                        w = w-gamma*(x*weight*(numpy.dot(x.T,w)*weight - trainE) + self.ridge_coef*w)
                    #Lasso
                    if self.lasso:
                        w = w-gamma*x*weight*(numpy.dot(x.T,w)*weight - trainE*weight)
                        for j in range(d):
                            w[j][0]=STF(w[j][0],gamma*self.lambda_lasso_coef)
            
                    # analyze the fit of the current model w
                    #Regular Matrix Multiplication
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
                        if all(mm<float(.1) for mm in madlist100):
                            break
                        
                        
          
        import matplotlib.pyplot as plt
        plt.plot(range(len(MAD)),MAD)
        plt.show()
        
        plt.plot(range(len(CC)),CC)
        plt.show()
        
        print('energy function fit:')
        print('\tCC = %0.2f' % CC[(len(MAD))-1])
        print('\tMAD = %0.2f' % MAD[(len(MAD))-1])
        
        plt.scatter(E_testinglist,E_estimate)
        plt.show()
        
        
        #Test to see how many features are zero...
        zerolist = []
        for z in range(len(w)):
            if w[z][0]==0:
                zerolist.append(z)
        print('Number of Features Equal to Zero in the Model = %s' % len(zerolist))
        
        sumlist = []
        for i in range(1000):
            sumlist.append(MAD[-i])
        print('Avg MAD of last 1000 terms = '+str(sum(sumlist)/1000))
        
        elapsedtime = time.time() - starttime 
        print('time = %0.2f' % elapsedtime)
        return w



class SGD_Online:
    
    def __init__(self,pdbfile,sequence_alignment_file,wt_seq,energy_scoring_script='score_sequence_amber.py',output_energy_file=False,custom_tag=False,pair_select_list=False,gamma_multiplier=2,regular_linear_regression=True,lasso=False,lambda_lasso_coef=.01,ridge=False,ridge_coef=.01,output_cc=False,output_mad=False):
        'Takes in a pdb file, sequence alignment file, wild type sequence, and energy scoring script'
        'We provide two scoring scripts score_sequence_amber.py and score_sequence_rosetta.py'
        'The rosetta script requires rosetta protein structure software'
        'The amber script requires openmm molecular modeling software as well as rosetta for energy minimization'
        'this class must run in command line to make use of energy scoring scripts which output energies to command line'        
        self.pdbfile = pdbfile
        self.energy_scoring_script = energy_scoring_script
        self.sequence_alignment_file = sequence_alignment_file
        self.pair_select_list = pair_select_list
        self.gamma_multiplier = gamma_multiplier
        self.wt_seq = wt_seq
        self.regular_linear_regression = regular_linear_regression
        self.lasso = lasso
        self.ridge = ridge
        self.lambda_lasso_coef = lambda_lasso_coef
        self.ridge_coef = ridge_coef
        self.output_energy_file = output_energy_file
        self.output_cc = output_cc
        self.output_mad = output_mad
        self.custom_tag = custom_tag

    def online_model(self,num_mutations=4,max_computation_time=252000,max_training_sets=50000,mad_cutoff=2):
        if not self.regular_linear_regression and not self.lasso and not self.ridge:
            raise ValueError('No model type selected, one of regular_linear_regression, lasso, or ridge must be true')
        if self.regular_linear_regression and self.lasso:
            raise ValueError('Cannot use both regular and lasso regression, select one or the other')
        if self.regular_linear_regression and self.ridge:
            raise ValueError('Cannot use both regular and ridge regression, select one or the other')
        if self.ridge and self.lasso:
            raise ValueError('Cannot use both ridge and lasso regression, select one or the other')
        if self.regular_linear_regression and self.lasso and self.ridge:
            raise ValueError('Cannot use regular, ridge, and lasso regression at the same time, select one or the other')
        
        
        from random import choice
        from subprocess import Popen,PIPE 
        import numpy
        import random
        import time
        starttime = time.time()
        
        
        import time
        
        max_time = max_computation_time
        

        
        alignmentfile = self.sequence_alignment_file
        

        trainingsets = max_training_sets
        madcutoff = mad_cutoff
        if not self.custom_tag:
            if '/' in self.pdbfile:
                pdbtag = self.pdbfile.split('/')[-1].split('.')[0]
            else:
                pdbtag = self.pdbfile
            if '\\' in self.pdbfile:
                pdbtag = self.pdbfile.split('\\')[-1].split('.')[0]        
        else: 
            pdbtag = self.pdbfile
        if self.custom_tag:
             pdbtag = self.custom_tag
            
        if self.output_energy_file:
            seqEfilename = 'ubiquitin_energies_'+str(pdbtag)+'.txt'

        
        AAs = ['A', 'C', 'D', 'E', 'F', 'G', 'H','I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W','Y'] #all 20 amino acids
        
        
        
        wildtype = self.wt_seq
        wtlist = [wildtype[i] for i in range(len(wildtype))]
        
        L = len(wildtype) # the length of the protein
        

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
        reference_seq = wtlist
        print(reference_seq)
        

        
        if self.pair_select_list:
            position_pairs = self.pair_select_list
        else:
            position_pairs = []
        
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
        CC = []
        MAD = []
        w = []
        Philist = []
        print('Stochastic Gradient Descent')
        print('# of pair interactions = %s, %s maximum training sequences, gamma = %s/((i+1)**0.5)' % (len(position_pairs),trainingsets, self.gamma_multiplier))
        
	if self.output_energy_file:    
            open(seqEfilename,'w')
        if self.output_cc:
            open('CC'+str(pdbtag)+'.txt','w')
        if self.output_mad:
            open('MAD'+str(pdbtag)+'.txt','w')
        
        
        for n in range(trainingsets):
            if time.time()-starttime<max_time:
                wildtype = self.wt_seq
                wtlist = [wildtype[i] for i in range(len(wildtype))]
                randAAlist = []
                for i in range(num_mutations):
                    randAAlist.append(random.randint(0,75))
                for mutantpositions in randAAlist:
                    wtlist[mutantpositions] = choice(overall_solution_space[mutantpositions])
                random_seq = ''.join(wtlist)
                cmd = 'python '+str(self.energy_scoring_script)+' %s' % random_seq #1. generate a terminal command as a string
                output = Popen(cmd,shell=True,stdout=PIPE).communicate() #2. send the command to the terminal shell
                seq_E = float(output[0]) #3. read the result and convert to a float
                
                print(random_seq, seq_E)
                if self.output_energy_file:
		    open(seqEfilename,'a').write(random_seq+','+str(seq_E)+'\n')
                E_list.append(seq_E)
                #Convert Seq to Phi Here
                seq_Phi = [1] 
                for i in range(L):
                    Phi_i = sigma2Phi(random_seq[i],i)
                    seq_Phi.extend(Phi_i)
                if self.pair_select_list:
                    for pair in position_pairs:
                        firstpos = pair[0]
                        secondpos = pair[1]
                        pairs_j = pairs(random_seq[firstpos],random_seq[secondpos],firstpos,secondpos)
                        seq_Phi.extend(pairs_j)
                else:
                    pass
                if n<100:
                    Philist.append(seq_Phi)
                    if n==0:
                        d = len(seq_Phi)
                        w = numpy.zeros(shape=(d,1)) #Starting vector for regular regression not L1
                        if self.lasso:
                            for f in range(1,len(Phi_i*L)):
                                w[f][0] = .01
                        print('Number of Variables in Model = '+str(d))
                if n>=100:
                    Philist.append(seq_Phi)
                    trainPhi = Philist[0]
                    del Philist[0]
                    E_testinglist = [E_list[-xx] for xx in range(1,101)]
                    E_testinglist.reverse()
                    trainE = E_list[n-100]
                    weight = 1
                    
                    gamma = self.gamma_multiplier/((n+1)**0.5) # have to do plus one cause first round i = 0
                    
                    x = numpy.array([trainPhi]).T # column vector
                    if self.regular_linear_regression:
                        w = w-gamma*x*weight*(numpy.dot(x.T,w)*weight - trainE*weight)
                    #Ridge
                    if self.ridge:
                        w = w-gamma*(x*weight*(numpy.dot(x.T,w)*weight - trainE) + self.ridge_coef*w)
                    #Lasso
                    if self.lasso:
                        w = w-gamma*x*weight*(numpy.dot(x.T,w)*weight - trainE*weight)
                        for j in range(d):
                            w[j][0]=STF(w[j][0],gamma*self.lambda_lasso_coef)
            

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
                    if self.output_cc:
                        open('CC'+str(pdbtag)+'.txt','a').write(str(cc)+'\n')
                    if self.output_mad:
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
        return w

trial = SGD_Static(pdbfile='C:/Users/Trent/Desktop/PythonFiles/Project/Package/State0.pdb',
            sequence_energy_file='C:/Users/Trent/Desktop/PythonFiles/Project/ubiquitin_State0_energies_amber_test.txt',
            sequence_alignment_file='C:/Users/Trent/Desktop/PythonFiles/Project/smallalignment.txt',
            wt_seq='MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG',
            gamma_multiplier=1,regular_linear_regression=False,lasso=True,
            ridge=False,ridge_coef=.01,lambda_lasso_coef=.0001)

trial.static_model()

#trial = SGD_Online('State0.pdb','smallalignment.txt','MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG',energy_scoring_script='score_sequence_amber.py')
#trial.online_model(num_mutations=4,max_computation_time=252000,max_training_sets=500,mad_cutoff=2)
