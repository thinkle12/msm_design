from __future__ import division, print_function
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 13:56:32 2016

@author: Trent
"""

from random import choice
from subprocess import Popen
import numpy
import random
import time
from sys import stdout

class SGD_Static:
    
    def __init__(self,pdbfile,sequence_energy_file,wt_seq,sequence_alignment_file=False,reduce_alignment=1,pair_select_list=False,pair_dist=4,gamma_multiplier=2,regular_linear_regression=True,lasso=False,lambda_lasso_coef=.01,ridge=False,ridge_coef=.01,custom_tag=False):
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
        self.custom_tag = custom_tag
        self.overall_solution_space=None
        self.energies=None
        self.sequences=None
        self.CC=None
        self.MAD=None
        self.model_parameters=None
        self.pair_dist=pair_dist
        self.reduce_alignment=reduce_alignment
	


    def static_model(self,mad_thresh=.01,model_output_text_file=False,display_plots=True):
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
        
        

        starttime = time.time()

        
        #max_time = 252000
        max_time = 5000000
        start_time = time.time()
        
        
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
             
        print('Filename Tag will be = '+str(pdbtag))
            
        
        L = len(self.wt_seq) # the length of the protein
        
        
        overall_solution_space = import_sequence_alignment(self.sequence_alignment_file,self.reduce_alignment,L)
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
        
        self.sequences=sequences
        self.energies=Energies
        

        if self.pair_select_list:
            position_pairs = pair_select(self.pdbfile,self.pair_dist,self.wt_seq)
        else:
            position_pairs = []
        

        
            
        
        
        
        # Sample a random sequence, evaluate its energy according to our energy fucntion, and store the sequence and energy data
        trainingsets = len(Energies)
        print('Number of training sets = '+str(len(Energies)))
        
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
                random_Seq_list.append(sequences[n])
                E_list.append(Energies[n])

                #Convert Seq to Phi Here
                seq_Phi = [1] 
                for i in range(L):
                    Phi_i = sigma2Phi(random_Seq_list[n][i],i,overall_solution_space,reference_seq)
                    seq_Phi.extend(Phi_i)
                if self.pair_select_list:
                    for pair in position_pairs:
                        firstpos = pair[0]
                        secondpos = pair[1]
                        pairs_j = pairs(random_Seq_list[n][firstpos],random_Seq_list[n][secondpos],firstpos,secondpos,overall_solution_space,reference_seq)
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
                        if all(mm<float(mad_thresh) for mm in madlist100):
                            break
           
        if model_output_text_file:
            open('w'+str(pdbtag)+'.txt','w')
            for vv in range(len(w)):
                open('w'+str(pdbtag)+'.txt','a').write(str(w[vv][0])+'\n')             
                  
                  
        self.CC=CC
        self.MAD=MAD
        self.model_parameters=w
        if display_plots:  
            import matplotlib.pyplot as plt
            plt.plot(range(len(MAD)),MAD)
            plt.show()
            
            plt.plot(range(len(CC)),CC)
            plt.show()
            
            
            
            plt.scatter(E_testinglist,E_estimate)
            plt.show()
        
        print('energy function fit:')
        print('\tCC = %0.2f' % CC[(len(MAD))-1])
        print('\tMAD = %0.2f' % MAD[(len(MAD))-1])
        
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
    
    def __init__(self,wt_seq,pdbfile,sequence_alignment_file=False,reduce_alignment=2,output_energy_file=False,
                 custom_tag=False,pair_select_list=False,pair_dist=4,gamma_multiplier=2,
                 regular_linear_regression=True,lasso=False,lambda_lasso_coef=.01,ridge=False,
                 ridge_coef=.01,output_cc=True,output_mad=True,
                 rosetta_score_path='/home/romeroroot/code/rosetta_src_2016.17.58663_bundle/main/source/bin/',
                 rosetta_database_path='/home/romeroroot/code/rosetta_src_2016.17.58663_bundle/main/database'):
                     
                     
        'Takes in a pdb file, sequence alignment file, wild type sequence, and energy scoring script'
        'We provide two scoring scripts score_sequence_amber.py and score_sequence_rosetta.py'
        'The rosetta script requires rosetta protein structure software'
        'The amber script requires openmm molecular modeling software as well as rosetta for energy minimization'
        'this class must run in command line to make use of energy scoring scripts which output energies to command line'        
        self.pdbfile = pdbfile
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
        self.overall_solution_space=None
        self.CC=None
        self.MAD=None
        self.model_parameters=None
        self.reduce_alignment = reduce_alignment
        self.pair_dist=pair_dist
        self.sequence_alignment_file=sequence_alignment_file
        self.rosetta_score_path = rosetta_score_path
        self.rosetta_database_path = rosetta_database_path


    def online_model(self,energy_scoring_function='rosetta',num_mutations=4,max_computation_time=252000,max_training_sets=50000,mad_cutoff=1,w_start_file=False):
        #Scoring function is either 'rosetta' or 'amber'
        
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
        
        

        starttime = time.time()
        
        
        
        max_time = max_computation_time
        
        if w_start_file:
            hh = open(w_start_file,'r+')
            newtext = hh.read()
            newtext = newtext.split('\n')
            del newtext[-1]
            w = []
            for numbers in newtext:
                w.append(float(numbers))


        wildtype = self.wt_seq
        wtlist = [wildtype[i] for i in range(len(wildtype))]
        
        L = len(wildtype)
        
        overall_solution_space = import_sequence_alignment(self.sequence_alignment_file,self.reduce_alignment,L)
        print(overall_solution_space)

        trainingsets = max_training_sets
        madcutoff = mad_cutoff
        
        
        #Here Creating the tag that follows to output files...
        #This is skipped if a custom_tag variable is defined
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
             
        print('Filename Tag will be = '+str(pdbtag))
            
        if self.output_energy_file:
            seqEfilename = 'ubiquitin_energies_'+str(pdbtag)+'.txt'

        
        
        
        
         # the length of the protein
        


        
        reference_seq = wtlist
        print(reference_seq)
        

        
        if self.pair_select_list:
            position_pairs = pair_select(self.pdbfile,self.pair_dist,self.wt_seq)
        else:
            position_pairs = []
        
        #Definitions
        
        # Sample a random sequence, evaluate its energy according to our energy fucntion, and store the sequence and energy data
        
        
        E_list = []
        CC = []
        MAD = []
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
                
                #cmd = 'python '+str(self.energy_scoring_script)+' %s' % random_seq #1. generate a terminal command as a string
                #output = Popen(cmd,shell=True,stdout=PIPE).communicate() #2. send the command to the terminal shell
                #seq_E = float(output[0]) #3. read the result and convert to a float
                
                'Add a new energy scoring function here if wanted'
                'Input should be random_seq and an associated pdbfile'
                'Should create a definition similar to the two energy function definitions provided'
                'if there is a need or want for a different energy function' 
                if energy_scoring_function=='amber':
                    seq_E = score_sequence_amber(random_seq,pdbfile=self.pdbfile,
                                                   rosetta_path = self.rosetta_score_path,
                                                   rosetta_db = self.rosetta_database_path)          
                
                if energy_scoring_function=='rosetta':
                    seq_E = score_sequence_rosetta(random_seq,pdbfile=self.pdbfile,
                                                   rosetta_path = self.rosetta_score_path,
                                                   rosetta_db = self.rosetta_database_path)          
                
                
                #This bypasses Amber not returning scores for certain sequences
                if seq_E=='Nan':
                    print('Sequence did not generate a score')
                    pass                
                
                if seq_E!='Nan':
                    #print(seq_E)
                    seq_E = float(seq_E) #3. read the result and convert to a float
                
                
    
    
                    print(random_seq, seq_E)
                    if self.output_energy_file:
                        open(seqEfilename,'a').write(random_seq+','+str(seq_E)+'\n')
                    E_list.append(seq_E)
                    #Convert Seq to Phi Here
                    seq_Phi = [1] 
                    for i in range(L):
                        Phi_i = sigma2Phi(random_seq[i],i,overall_solution_space,reference_seq)
                        seq_Phi.extend(Phi_i)
                    if self.pair_select_list:
                        for pair in position_pairs:
                            firstpos = pair[0]
                            secondpos = pair[1]
                            pairs_j = pairs(random_seq[firstpos],random_seq[secondpos],firstpos,secondpos,overall_solution_space,reference_seq)
                            seq_Phi.extend(pairs_j)
                    else:
                        pass
                    if len(E_list)-1<100:
                        Philist.append(seq_Phi)
                        if n==0:
                            d = len(seq_Phi)
                            if not w_start_file:
                                w = numpy.zeros(shape=(d,1)) #Starting vector for regular regression not L1
                                if self.lasso:
                                    for f in range(1,len(Phi_i*L)):
                                        w[f][0] = .01
                            print('Number of Variables in Model = '+str(d))
                    if len(E_list)-1>=100:
                        Philist.append(seq_Phi)
                        trainPhi = Philist[0]
                        del Philist[0]
                        E_testinglist = [E_list[-xx] for xx in range(1,101)]
                        E_testinglist.reverse()
                        trainE = E_list[len(E_list)-1-100]
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
        
        self.CC=CC
        self.MAD=MAD
        self.model_parameters=w
        
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
        
        
    
    
    
def score_sequence_amber(newseq=None,pdbfile=None,rosetta_path = '/home/romeroroot/code/rosetta_src_2016.17.58663_bundle/main/source/bin/',rosetta_db = '/home/romeroroot/code/rosetta_src_2016.17.58663_bundle/main/database'):   
    'Provide rosetta path to fixbb.linuxgccrelease'
    from random import choice
    from os import remove
    import simtk.openmm.app as app
    import simtk.openmm as op
    import simtk.unit as unit

    def rand_tag():
        '''creates a unique tag that can be used to identify the current files a script is working with.  necessary for running the same program multiple times in the same directory'''
        alpha = 'abcdefghijklmnopqrstuvwxyz0123456789'
        tag = ''.join([choice(alpha) for i in range(15)])
        return tag
    
    
    
    def thread_repack_score_amber(pdbfile,newseq):
    
        tag = rand_tag() # generate a random tag to associate with this run
    
        # generate a resfile to mutate the pdb
        resfile = """
        #header
        NATAA # keep this so all sites aren't designed
        USE_INPUT_SC # this to also consider the WT rotamer, if we're mutating WT to WT
    
        start
    
        #body\n"""
    
        for i in range(len(newseq)):
            resfile += '%i    A    PIKAA %s \n'%(i+1,newseq[i])
    
        open('resfile_'+tag+'.res','w').write(resfile)
    
    
        # the options for the run
        options = ['nice',
                   rosetta_path+'fixbb.linuxgccrelease', # fixbb is the program used for threading a repacking
                   '-s '+pdbfile, 
                   '-database '+rosetta_db, # good to explicitly define location of rosetta DB
                   '-resfile resfile_'+tag+'.res', # a resfile to specify mutations
                   '-out:suffix _'+tag] # this adds a suffix to the score output file: score_[RNDTAG].sc
    
        # run fixbb
        Popen(options,stdout=open('/dev/null','w')).wait() # send stdout to the trash, report stderr
    
    
        # read the output
        scorefile = open('score_%s.sc'%tag).read().split('\n')
        names = scorefile[1].split(':')[1].split()
        values = scorefile[2].split(':')[1].split()
        score = dict((names[i],float(values[i])) for i in range(len(names)-1)) # use i-1 because the last entry is the filename
        pdbname = (pdbfile[:-4]+'_'+tag+'_0001.pdb').split('/')[-1]
    
    
        # now load into openMM and score with amber
        pdb = app.PDBFile(pdbname)
        forcefield = app.ForceField('amber03.xml','amber03_obc.xml')
        system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff, constraints=None)
        integrator = op.LangevinIntegrator(300*unit.kelvin, 1/unit.picosecond, 1e-9*unit.picoseconds)
        simulation = app.Simulation(pdb.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)
    
    
        #Must do try/except here because Amber occasionally throws an Exception error and does not return a score
        try:
            simulation.minimizeEnergy()
            state = simulation.context.getState(getPositions=True, getEnergy=True)
            score = state.getPotentialEnergy()/unit.kilojoule_per_mole
    
        except Exception:
            score = 'Nan'
            pass
    
    
        #delete the evidence
        remove('score_%s.sc'%tag)
        remove('resfile_'+tag+'.res')
        remove(pdbname)
    
        return score
    
    
    score = thread_repack_score_amber(pdbfile,newseq)
    return score
    
    
    
def score_sequence_rosetta(newseq,pdbfile,rosetta_path = '/home/romeroroot/code/rosetta_src_2016.17.58663_bundle/main/source/bin/',rosetta_db = '/home/romeroroot/code/rosetta_src_2016.17.58663_bundle/main/database'):
    from random import choice
    from subprocess import Popen
    from os import remove

    def rand_tag():
        '''creates a unique tag that can be used to identify the current files a script is working with.  necessary for running the same program multiple times in the same directory'''
        alpha = 'abcdefghijklmnopqrstuvwxyz0123456789'
        tag = ''.join([choice(alpha) for i in range(15)])
        return tag
    
    
    
    def fast_thread_repack_score(pdbfile,newseq,savepdb=False):
        """this function inputs a pdbfile name and a sequence, mutates the pdb file, repacks, and returns the scores as a dict
        REQUIRES: pdbfile to be consecutivly numbered from 1 to N, and len(newseq) = N. Used when we want to score seqence variants many times
        """
    
        tag = rand_tag() # generate a random tag to associate with this run
    
    
        # generate a resfile to mutate the pdb
        resfile = """
        #header
        NATAA # keep this so all sites aren't designed
        USE_INPUT_SC # this to also consider the WT rotamer, if we're mutating WT to WT
    
        start
    
        #body\n"""
    
        for i in range(len(newseq)):
            resfile += '%i    A    PIKAA %s \n'%(i+1,newseq[i])
    
        open('resfile_'+tag+'.res','w').write(resfile)
    
    
        # the options for the run
        options = ['nice',
                   rosetta_path+'fixbb.linuxgccrelease', # fixbb is the program used for threading a repacking
                   '-s '+pdbfile, 
                   '-database '+rosetta_db, # good to explicitly define location of rosetta DB
                   '-resfile resfile_'+tag+'.res', # a resfile to specify mutations
                   '-out:suffix _'+tag] # this adds a suffix to the score output file: score_[RNDTAG].sc
    
        # run fixbb
        Popen(options,stdout=open('/dev/null','w')).wait() # send stdout to the trash, report stderr
    
    
        # read the output
        scorefile = open('score_%s.sc'%tag).read().split('\n')
        names = scorefile[1].split(':')[1].split()
        values = scorefile[2].split(':')[1].split()
        score = dict((names[i],float(values[i])) for i in range(len(names)-1)) # use i-1 because the last entry is the filename
    
    
        #delete the evidence
        remove('score_%s.sc'%tag)
        remove('resfile_'+tag+'.res')
    
        pdbname = (pdbfile[:-4]+'_'+tag+'_0001.pdb').split('/')[-1]
        if savepdb:
            return score,pdbname
        else:
            remove(pdbname)
            return score

    score = fast_thread_repack_score(pdbfile,newseq)
    return score['total_score']


def pair_select(pdbfile,dist,wt_seq):
    "Here we select the maximum distance for a pair of AAs to be considered important"
    angstrom_distance = dist
    #########################################################################
    #this should really use the 'state' pdb files....
    L = len(wt_seq)
    
    
    
    
    
    #AAs = ['A', 'C', 'D', 'E', 'F', 'G', 'H','I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W','Y'] #all 20 amino acids
    
    

    pdbdata = open(pdbfile).read()
    
    pdblines = pdbdata.split('\n')
    
    coordinates = []
    for line in pdblines:
        if line[:4]=='ATOM':
            coordinates.append(line)
    
    #"Why is there a lot of extra sequences in the pdb file"
    #for i in range(4123):
    #    coordinates.pop()
    
    "Here we make lists of all atoms x,y,z coordinates and pair them by which AA they are in into lists, each list corresponds to the atoms in each AA"
    xcoors = {}
    for q in range(1,L+1):
        xcoors['AAxcoor'+str(q)] = []
    
    for j in range(L+1):
        for line in coordinates:
            if int(line[22:26])==j:
                x = line[30:38]
                xnum = float(x)
                xcoors['AAxcoor'+str(j)].append(xnum)
    
    ycoors = {}
    for q in range(1,L+1):
        ycoors['AAycoor'+str(q)] = []
    
    for k in range(L+1):
        for line in coordinates:
            if int(line[22:26])==k:
                y = line[38:46]
                ynum = float(y)
                ycoors['AAycoor'+str(k)].append(ynum)
    
    zcoors = {}
    for q in range(1,L+1):
        zcoors['AAzcoor'+str(q)] = []
    
    for s in range(L+1):
        for line in coordinates:
            if int(line[22:26])==s:
                z = line[46:54]
                znum = float(z)
                zcoors['AAzcoor'+str(s)].append(znum)
                
    
    "Pair all xyz together"
    listofxyz = {}
    for q in range(1,L+1):
        listofxyz['listofxyz'+str(q)] = []
    
    "XYZ points of all atoms sorted into lists of corresponding AAs"
    for q in range(1,L+1):
        for j in range(len(xcoors['AAxcoor'+str(q)])):
            listofxyz['listofxyz'+str(q)].append((xcoors['AAxcoor'+str(q)][j],ycoors['AAycoor'+str(q)][j],zcoors['AAzcoor'+str(q)][j]))
        
    "Here we find the distance between every atom in the protein separated by the atoms that belong to specific AAs"
    atomtoatom = {}
    for i in range(1,L+1):
        for j in range(1,L+1):
            if i<j:
                atomtoatom['atomtoatom'+str(i)+'-'+str(j)] = []
                for k in range(len(listofxyz['listofxyz'+str(i)])):
                    for f in range(len(listofxyz['listofxyz'+str(j)])):
                        atomtoatom['atomtoatom'+str(i)+'-'+str(j)].append(numpy.linalg.norm(numpy.array(listofxyz['listofxyz'+str(i)])[k]-numpy.array(listofxyz['listofxyz'+str(j)])[f]))
    
    "Here we take the minimum distance between atoms between all AAs"
    minimumdistances = {}
    for i in range(1,L+1):
        for j in range(1,L+1):
            if i<j and atomtoatom['atomtoatom'+str(i)+'-'+str(j)] != []: #New shit here...
                minimumdistances['minimumdistances'+str(i)+'-'+str(j)] = []
                minimumdistances['minimumdistances'+str(i)+'-'+str(j)] = min(atomtoatom['atomtoatom'+str(i)+'-'+str(j)])
    
    "Here we determine which pairs are close enough by setting an angstrom_distance constraint"
    closeatoms = []
    for i in range(1,L+1):
        for j in range(1,L+1):
            if all([i<j and minimumdistances['minimumdistances'+str(i)+'-'+str(j)]<float(angstrom_distance)]):
            #if all([i<j and minimumdistances['minimumdistances'+str(i)+'-'+str(j)]<float(angstrom_distance) and abs(i-j)!=1]):
                closeatoms.append((i,j))
    
    "Here we subtract of (1,1) from each pair to match the indecies of the regression format"  
    

              
    newcloseatoms = numpy.array(closeatoms) - (1,1)
    newcloseatoms = map(tuple, newcloseatoms)
    
    print(len(newcloseatoms))
    return newcloseatoms


def import_sequence_alignment(seq_align_file,reduce_alignment,seq_len):
    if seq_align_file:
        h = open(seq_align_file,'r+')
        text = h.read()
        text = text.split('\n')
        del text[-1]
        
        overall_solution_space = []
        for topaminoacids in text:
            tempsolspace = []
            for i in range(len(topaminoacids)):
                tempsolspace.append(topaminoacids[i])
            overall_solution_space.append(tempsolspace)
            
        #This part deletes the last AA at each position to have one less AA per position
        #Take this out Eventually
        if reduce_alignment:
            for i in range(reduce_alignment):
                 for top5 in overall_solution_space:
                     del top5[-1]  
            return overall_solution_space
        else:
            return overall_solution_space
    else:
        AAs = ['A', 'C', 'D', 'E', 'F', 'G', 'H','I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W','Y'] #all 20 amino acids
        for i in range(len(seq_len)):
            overall_solution_space.append(AAs)
        return overall_solution_space
        
#Definitions
def sigma2Phi(sigma,i,overall_solution_space,reference_seq):
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

def pairs(aa1,aa2,i,h,overall_solution_space,reference_seq):
    #This generates the pairs associated with currently solution space 
    #Similar to sigma2Phi but for the pairs
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
    


#Used in the lasso function
def STF(a,z):
    if a>z:
        return a-z
    if a<-z:
        return a+z
    if -z<=a<=z:
        return 0




