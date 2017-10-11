# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 23:26:11 2017

@author: Trent
"""



import numpy
import os
import re
import time
from sys import stdout
import random
import math

import msm_model


#Optimization for pairs is not yet supported
#Hold off on pairs for now...


class Sequence_Optimization:
    
    def __init__(self,model_file_dir,wt_seq,sequence_alignment_file=False,reduce_alignment=1,mad_file_dir=False,cc_file_dir=False,display_plot=True,pair_interaction_dist=False,pdb_dir=False):
        self.model_file_dir = model_file_dir
        self.mad_file_dir = mad_file_dir
        self.cc_file_dir = cc_file_dir
        self.wt_seq = wt_seq
        self.sequence_alignment_file = sequence_alignment_file
        self.display_plot = display_plot
        self.best_results = None
        self.reduce_alignment = reduce_alignment
#        if pair_interaction_dist and pdb_dir:
#            
#            pair_interaction_list=Msm_Design.pair_select(pdbfile,pair_interaction_dist,wt_seq)
#            self.pair_interaction_list = pair_interaction_list
        
    def import_model_files(self):

        newfilelist = os.listdir(self.model_file_dir)
        
        newnamelist = []
        for newnames in newfilelist:
            newnamelist.append(re.sub("[^0-9]",'',newnames))
        
        wStates = {}
        for q in range(len(newfilelist)):
            wStates['wState'+str(newnamelist[q])] = []
        
        
        for j in range(len(newfilelist)):
            hh = open(self.model_file_dir+str(newfilelist[j]),'r+')
            newtext = hh.read()
            newtext = newtext.split('\n')
            del newtext[-1]
            newfloatedlist = []
            for numbers in newtext:
                newfloatedlist.append(float(numbers))
            wStates['wState'+str(newnamelist[j])] = numpy.array(newfloatedlist)
        return wStates
        

        
    def test_model_performance(self):
        if self.mad_file_dir:
            filelist = os.listdir(self.mad_file_dir)
            
            #print filelist
            
            namelist = []
            for names in filelist:
                namelist.append(re.sub("[^0-9]",'',names))
            
            MADStates = {}
            for q in range(len(filelist)):
                MADStates['MADState'+str(namelist[q])] = []
                
                
            sumlist = []
            for j in range(len(filelist)):
                h = open(self.mad_file_dir+str(filelist[j]),'r+')
                text = h.read()
                text = text.split('\n')
                #text = pickle.load(open('seq_energies_state0_ubiquitin.p'))
                del text[-1]
                floatedlist = []
                for items in text:
                    try:
                        floatedlist.append(float(items))
                    except ValueError:
                        pass
                if floatedlist != []:
                    MADStates['MADState'+str(namelist[j])] = floatedlist
                    thousandlist = []
                    if len(floatedlist)<1000:
                        for i in range(len(floatedlist)):
                            thousandlist.append(floatedlist[i])
                        sumlist.append(sum(thousandlist)/len(thousandlist))
                    if len(floatedlist)>=1000:
                        for i in range(1,1001):
                            thousandlist.append(floatedlist[-i])
                        sumlist.append([sum(thousandlist)/1000])
                sumlist[-1].append(str(filelist[j]))
            print sumlist
            return sumlist
        else:
            print 'No MAD Directory Set'
            
    def generate_seq_energy(self,seq,state_number):
        'Pair models not yet supported'
        wStates=self.import_model_files()
        overall_solution_space = msm_model.import_sequence_alignment(self.sequence_alignment_file,self.reduce_alignment,len(self.wt_seq))
        
        wildtype = self.wt_seq
        wtlist = [wildtype[i] for i in range(len(wildtype))]
        reference_seq = wtlist        
        
                
        
        def binary_conversion(sequence):
            seq_Phi = [1] 
            for i in range(len(self.wt_seq)):
                Phi_i = msm_model.sigma2Phi(sequence[i],i,overall_solution_space,reference_seq)
                seq_Phi.extend(Phi_i)
            #    for pair in position_pairs:
            #        firstpos = pair[0]
            #        secondpos = pair[1]
            #        pairs_j = pairs(lowESeqtest[n][firstpos],lowESeqtest[n][secondpos],firstpos,secondpos)
            #        seq_Phi.extend(pairs_j)
            return seq_Phi
        
        
        
        rnd_list = []
        for j in range(len(seq)):
            if seq[j] not in overall_solution_space[j]:
                print 'Residue '+str(seq[j])+' at position '+str(j)+' not in amino acid space '+str(overall_solution_space[j])
                rnd_list.append(j)                
                break

        if len(rnd_list)<=0:      
            binary_seq=binary_conversion(seq)
            energy = numpy.dot(binary_seq,wStates['wState'+str(state_number)])
            return energy

    def single_state_design(self,states_to_optimize_for_by_state_number,t=298.15,k=8.31447e-3,mut_thresh=10,num_mutation_pathway_passes=500,inner_loop_trials=1000,mean_shift=True,plot_separate=True,plot_together=False,pdf_output_plot=False,write_opt_file_name=False):
        'Pair distances not supported'


        mutation_threshold = len(self.wt_seq)-mut_thresh
        
        

        
        starttime = time.time()
        
        
        if self.mad_file_dir:
            print self.test_model_performance()
        else:
            pass
        
        wStates = self.import_model_files()
                
        modelsfiles = os.listdir(self.model_file_dir)
        
        modelsfilelist = []
        for newnames in modelsfiles:
            modelsfilelist.append(re.sub("[^0-9]",'',newnames))
            
        
        
        L = len(self.wt_seq)        
        # this is the amino acids allowed at each postion

        overall_solution_space = msm_model.import_sequence_alignment(self.sequence_alignment_file,self.reduce_alignment,L)
        print overall_solution_space        
        wildtype = self.wt_seq
        wtlist = [wildtype[i] for i in range(len(wildtype))]
        reference_seq = wtlist
        print(reference_seq)
        
            
            
        
        
        
        def binary_conversion(sequence):
            seq_Phi = [1] 
            for i in range(L):
                Phi_i = Msm_Design.sigma2Phi(sequence[i],i,overall_solution_space,reference_seq)
                seq_Phi.extend(Phi_i)
#                for pair in position_pairs:
#                    firstpos = pair[0]
#                    secondpos = pair[1]
#                    pairs_j = pairs(sequence[n][firstpos],sequence[n][secondpos],firstpos,secondpos)
#                    seq_Phi.extend(pairs_j)
            return seq_Phi
        
        
        
        
        
        
        def single_probability_estimation(binary, State_to_optimize, mean):
                E_estimate = math.exp(-(numpy.dot(binary,wStates['w'+str(State_to_optimize)])-mean)/(k*t))
                return E_estimate
                
        def multi_probability_estimation(binary, mean):
            State_E_list = []
            for a in range(len(modelsfilelist)):
                if modelsfilelist[a]!='98' and modelsfilelist[a]!='45':
                    E_estimate_mult = math.exp(-(numpy.dot(binary,wStates['wState'+str(modelsfilelist[a])])-mean)/(k*t))
                    State_E_list.append(E_estimate_mult)
            return sum(State_E_list)
        
        def energy_mean_calc(binary):
            State_E_list = []
            for a in range(len(modelsfilelist)):
                if modelsfilelist[a]!='98' and modelsfilelist[a]!='45':
                    E_mult = numpy.dot(binary,wStates['wState'+str(modelsfilelist[a])])
                    State_E_list.append(E_mult)
            return numpy.mean(State_E_list)
        
        
        if write_opt_file_name:
            opt_filename = write_opt_file_name
            
            open(opt_filename,'w')
            
            open(opt_filename,'a').write('Probability'+','+str('Sequence')+','+str('chosen_state')+'\n')
            
        #state_opt = ['State12']
        
        #Might need to take out WT mutations from overall_solution_space to avoid mutating to wt
        best_seq_perstate_permutnum = {}
        if pdf_output_plot:
            from matplotlib.backends.backend_pdf import PdfPages
            pp = PdfPages(pdf_output_plot)
        for state_num in states_to_optimize_for_by_state_number:
            chosen_state = 'State'+str(state_num)
            print chosen_state
            rand_seq_list = []
            boltz_list = []
            objective_func_sequence_dict = {}
            for muts in range(L-mutation_threshold+1):
                objective_func_sequence_dict['Mutations_'+str(muts)] = []
            #Test Optimization 
            overall = []
            for nn in range(num_mutation_pathway_passes):
                trial = []
                for n in range(inner_loop_trials):
                    if len(trial)<=mut_thresh:
                        if n==0:
                            aminoacidchain = [wildtype[i] for i in range(len(wildtype))]
                            random_seq = ''.join(aminoacidchain)
                            rand_seq_list.append(aminoacidchain)
                            Phi = binary_conversion(random_seq)
                            if mean_shift:
                                mean_energy = energy_mean_calc(Phi)
                            else:
                                mean_energy = 0
                            single_E = single_probability_estimation(Phi,chosen_state,mean_energy)
                            multi_E = multi_probability_estimation(Phi,mean_energy)
                            boltz_dist = single_E/multi_E
                            boltz_list.append(boltz_dist)
                            Mut_num = L-len([a for a, b in zip(rand_seq_list[-1], reference_seq) if a == b])
                            objective_func_sequence_dict['Mutations_'+str(Mut_num)].append([boltz_dist,random_seq])
                            trial.append(boltz_dist)
                        if n>0:
                            randAAlist = []
                            for i in range(1):
                                randAAlist.append(random.randint(0,L-1))
                            if len([a for a, b in zip(rand_seq_list[-1], reference_seq) if a == b])<=mutation_threshold:
                                aachaintobemutated = list(aminoacidchain)
                            if len([a for a, b in zip(rand_seq_list[-1], reference_seq) if a == b])>mutation_threshold:
                                aachaintobemutated = list(rand_seq_list[-1])
                            for mutantpositions in randAAlist:
                                aachaintobemutated[mutantpositions] = random.choice(overall_solution_space[mutantpositions])
                            random_seq = ''.join(aachaintobemutated)
                            listseq = list(random_seq)
                            if len([a for a, b in zip(rand_seq_list[-1], reference_seq) if a == b])<=mutation_threshold:
                                rand_seq_list.append(aminoacidchain)
                                Phi = binary_conversion(rand_seq_list[-1])
                                if mean_shift:
                                    mean_energy = energy_mean_calc(Phi)
                                else:
                                    mean_energy = 0
                                single_E = single_probability_estimation(Phi,chosen_state,mean_energy)
                                multi_E = multi_probability_estimation(Phi,mean_energy)
                                boltz_dist = single_E/multi_E
                                boltz_list.append(boltz_dist)
                
                            if len([a for a, b in zip(rand_seq_list[-1], reference_seq) if a == b])>mutation_threshold:
                                rand_seq_list.append(listseq)
                                Phi = binary_conversion(rand_seq_list[-1])
                                if mean_shift:
                                    mean_energy = energy_mean_calc(Phi)
                                else:
                                    mean_energy = 0
                                single_E = single_probability_estimation(Phi,chosen_state,mean_energy)
                                multi_E = multi_probability_estimation(Phi,mean_energy)
                                boltz_dist = single_E/multi_E
                                boltz_list.append(boltz_dist)
                                Mut_num = L-len([a for a, b in zip(rand_seq_list[-1], reference_seq) if a == b])
    
                                if boltz_list[-1]>boltz_list[-2]:
                                    trial.append(boltz_dist)
                                    objective_func_sequence_dict['Mutations_'+str(Mut_num)].append([boltz_dist,''.join(listseq)])
                                    pass
                                if boltz_list[-1]<=boltz_list[-2]:
                                    del rand_seq_list[-1]
                                    del boltz_list[-1]
                                
                                
                                
                        
                        
                        stdout.write('\rMutation Pathway Number = %s' % str(nn+1)) #This prints the current training set index on same line as an update
                        stdout.flush()
                        
                overall.append(trial)
            best_seq_perstate_permutnum[str(chosen_state)] = []
            for muts in range(L-mutation_threshold+1):
                try:
                    best_seq_perstate_permutnum[str(chosen_state)].append(max(objective_func_sequence_dict['Mutations_'+str(muts)]))
                except ValueError:
                    pass
            print '\n'

            

            
	    
            if self.display_plot and plot_separate:
                import matplotlib.pyplot as plt
                listt = []
                for i in range(len(best_seq_perstate_permutnum[str(chosen_state)])):
                    listt.append([i,best_seq_perstate_permutnum[str(chosen_state)][i][0]])
                states=[x[0] for x in listt]
                probs=[math.log(x[1]) for x in listt]
                plt.scatter(states,probs)
                plt.title('Single State Optimization for State '+str(state_num))
                if pdf_output_plot:                
                    pp.savefig()
                #plt.show()                  
                plt.close()
                

            
            if write_opt_file_name:
                for i in range(len(best_seq_perstate_permutnum[str(chosen_state)])):
                    open(opt_filename,'a').write(str(best_seq_perstate_permutnum[str(chosen_state)][i][0])+','+str(best_seq_perstate_permutnum[str(chosen_state)][i][1])+','+str(chosen_state)+','+str('num_Mutations = ')+str(i)+'\n')
           
        if pdf_output_plot:
            pp.close()    
        self.best_results=best_seq_perstate_permutnum 
        


        
        if self.display_plot and plot_together:
            import matplotlib.pyplot as plt
            for stuff in best_seq_perstate_permutnum:
                listt = []
                for i in range(len(best_seq_perstate_permutnum[stuff])):
                    listt.append([i,best_seq_perstate_permutnum[stuff][i][0]])
                states=[x[0] for x in listt]
                probs=[math.log(x[1]) for x in listt]
                plt.title('Log Shifted Probability Plot'+'\n'+'For Most Probable Sequence Per State Per Mutation Number')
                plt.xlabel('Mutation Number')
                plt.ylabel('Log of Probability')
                plt.scatter(states,probs)
                

        

        elapsedtime = time.time() - starttime 
        print('time = %0.2f' % elapsedtime)
        return best_seq_perstate_permutnum
        
        
        
    def two_state_design(self,states_to_optimize_for_by_state_number,t=298.15,k=8.31447e-3,mut_thresh=10,num_mutation_pathway_passes=100,inner_loop_trials=1000,mean_shift=True,plot_separate=True,plot_together=False,pdf_output_plot=False,write_opt_file_name=False):
        
        
        mutation_threshold = len(self.wt_seq)-mut_thresh
        

        
        starttime = time.time()
        
        
        if self.mad_file_dir:
            print self.test_model_performance()
        else:
            pass
        
        wStates = self.import_model_files()
                
        modelsfiles = os.listdir(self.model_file_dir)
        
        modelsfilelist = []
        for newnames in modelsfiles:
            modelsfilelist.append(re.sub("[^0-9]",'',newnames))
            
        
        
        L = len(self.wt_seq)        
        # this is the amino acids allowed at each postion

        overall_solution_space = msm_model.import_sequence_alignment(self.sequence_alignment_file,self.reduce_alignment,len(self.wt_seq))
        
        wildtype = self.wt_seq
        wtlist = [wildtype[i] for i in range(len(wildtype))]
        reference_seq = wtlist
        print(reference_seq)        
        

            
      
        
        def binary_conversion(sequence):
            seq_Phi = [1] 
            for i in range(L):
                Phi_i = msm_model.sigma2Phi(sequence[i],i,overall_solution_space,reference_seq)
                seq_Phi.extend(Phi_i)
            #    for pair in position_pairs:
            #        firstpos = pair[0]
            #        secondpos = pair[1]
            #        pairs_j = pairs(lowESeqtest[n][firstpos],lowESeqtest[n][secondpos],firstpos,secondpos)
            #        seq_Phi.extend(pairs_j)
            return seq_Phi
        
        
        
        
        
        
        def single_probability_estimation(binary, State_to_optimize, mean):
                E_estimate = math.exp(-(numpy.dot(binary,wStates['w'+str(State_to_optimize)])-mean)/(k*t))
                return E_estimate
                
        def multi_probability_estimation(binary, mean):
            State_E_list = []
            for a in range(len(modelsfilelist)):
                if modelsfilelist[a]!='98' and modelsfilelist[a]!='45':
                    E_estimate_mult = math.exp(-(numpy.dot(binary,wStates['wState'+str(modelsfilelist[a])])-mean)/(k*t))
                    State_E_list.append(E_estimate_mult)
            return sum(State_E_list)
        
        def energy_mean_calc(binary):
            State_E_list = []
            for a in range(len(modelsfilelist)):
                if modelsfilelist[a]!='98' and modelsfilelist[a]!='45':
                    E_mult = numpy.dot(binary,wStates['wState'+str(modelsfilelist[a])])
                    State_E_list.append(E_mult)
            return numpy.mean(State_E_list)
        
        



        if write_opt_file_name:
            opt_filename = write_opt_file_name
        
            open(opt_filename,'w')
        
            open(opt_filename,'a').write('Objective_Function_Value'+','+str('Sequence')+','+str('chosen_state1')+','+str('chosen_state2')+'\n')
        

        best_seq_perstate_permutnum = {}
        
        if pdf_output_plot:
            from matplotlib.backends.backend_pdf import PdfPages
            pp = PdfPages(pdf_output_plot)
        #Might need to take out WT mutations from overall_solution_space to avoid mutating to wt
        #Test Optimization
        for jj in range(len(states_to_optimize_for_by_state_number)):
            state_num1 = states_to_optimize_for_by_state_number[jj][0]
            state_num2 = states_to_optimize_for_by_state_number[jj][1]
            chosen_state1 = 'State'+str(state_num1)
            chosen_state2 = 'State'+str(state_num2)
            print chosen_state1
            print chosen_state2
            two_state_opt_list = []
            boltz_list1 = []
            boltz_list2 = []
            rand_seq_list = []
            objective_func_sequence_dict = {}
            for muts in range(L-mutation_threshold+1):
                objective_func_sequence_dict['Mutations_'+str(muts)] = []
    
            overall = []
            for nn in range(num_mutation_pathway_passes):
                trial = []
                for n in range(inner_loop_trials):
                    if len(trial)<=mut_thresh:
                        if n==0:
                            aminoacidchain = [wildtype[i] for i in range(len(wildtype))]
                            random_seq = ''.join(aminoacidchain)
                            rand_seq_list.append(aminoacidchain)
                            Phi = binary_conversion(random_seq)
                            mean_energy = energy_mean_calc(Phi)
                            single_E1 = single_probability_estimation(Phi,chosen_state1,mean_energy)
                            single_E2 = single_probability_estimation(Phi,chosen_state2,mean_energy)
                            multi_E = multi_probability_estimation(Phi,mean_energy)
                            boltz_dist1 = single_E1/multi_E
                            boltz_dist2 = single_E2/multi_E
                            boltz_list1.append(boltz_dist1)
                            boltz_list2.append(boltz_dist2)
                            opt_func = (boltz_dist1-.5)**2+(boltz_dist2-.5)**2
                            two_state_opt_list.append(opt_func)
                            Mut_num = L-len([a for a, b in zip(rand_seq_list[-1], reference_seq) if a == b])
                            objective_func_sequence_dict['Mutations_'+str(Mut_num)].append([opt_func,''.join(aminoacidchain),boltz_dist1,boltz_dist2])
                            trial.append(opt_func)                    
                        if n>0:
                            randAAlist = []
                            for i in range(1):
                                randAAlist.append(random.randint(0,L-1))
                            if len([a for a, b in zip(rand_seq_list[-1], reference_seq) if a == b])<=mutation_threshold:
                                aachaintobemutated = list(aminoacidchain)
                            if len([a for a, b in zip(rand_seq_list[-1], reference_seq) if a == b])>mutation_threshold:
                                aachaintobemutated = list(rand_seq_list[-1])
                            for mutantpositions in randAAlist:
                                aachaintobemutated[mutantpositions] = random.choice(overall_solution_space[mutantpositions])
                            random_seq = ''.join(aachaintobemutated)
                            listseq = list(random_seq)
                            if len([a for a, b in zip(rand_seq_list[-1], reference_seq) if a == b])<=mutation_threshold:
                                rand_seq_list.append(aminoacidchain)
                                Phi = binary_conversion(rand_seq_list[-1])
                                mean_energy = energy_mean_calc(Phi)
                                single_E1 = single_probability_estimation(Phi,chosen_state1,mean_energy)
                                single_E2 = single_probability_estimation(Phi,chosen_state2,mean_energy)
                                multi_E = multi_probability_estimation(Phi,mean_energy)
                                boltz_dist1 = single_E1/multi_E
                                boltz_dist2 = single_E2/multi_E
                                boltz_list1.append(boltz_dist1)
                                boltz_list2.append(boltz_dist2)
                                opt_func = (boltz_dist1-.5)**2+(boltz_dist2-.5)**2
                                two_state_opt_list.append(opt_func) 
                            if len([a for a, b in zip(rand_seq_list[-1], reference_seq) if a == b])>mutation_threshold:
                                rand_seq_list.append(listseq)
                                Phi = binary_conversion(rand_seq_list[-1])
                                mean_energy = energy_mean_calc(Phi)
                                single_E1 = single_probability_estimation(Phi,chosen_state1,mean_energy)
                                single_E2 = single_probability_estimation(Phi,chosen_state2,mean_energy)
                                multi_E = multi_probability_estimation(Phi,mean_energy)
                                boltz_dist1 = single_E1/multi_E
                                boltz_dist2 = single_E2/multi_E
                                boltz_list1.append(boltz_dist1)
                                boltz_list2.append(boltz_dist2)
                                opt_func = (boltz_dist1-.5)**2+(boltz_dist2-.5)**2
                                two_state_opt_list.append(opt_func)
                                Mut_num = 76-len([a for a, b in zip(rand_seq_list[-1], reference_seq) if a == b])
                                if (boltz_list1[-1]>boltz_list1[-2] and boltz_list2[-1]>boltz_list2[-2]):
                                #if two_state_opt_list[-1]<=two_state_opt_list[-2]:
                                    trial.append(opt_func)
                                    objective_func_sequence_dict['Mutations_'+str(Mut_num)].append([opt_func,''.join(listseq),boltz_dist1,boltz_dist2])
                                    pass
                                if (boltz_list1[-1]<=boltz_list1[-2] and boltz_list2[-1]<=boltz_list2[-2]) or (boltz_list1[-1]>boltz_list1[-2] and boltz_list2[-1]<=boltz_list2[-2]) or (boltz_list1[-1]<=boltz_list1[-2] and boltz_list2[-1]>boltz_list2[-2]):
                                #if two_state_opt_list[-1]>two_state_opt_list[-2]:                                
                                    del rand_seq_list[-1]
                                    del two_state_opt_list[-1]
                                    del boltz_list1[-1]
                                    del boltz_list2[-1]
        
        
                        stdout.write('\rMutation Pathway Number = %s' % str(nn+1)) #This prints the current training set index on same line as an update
                        stdout.flush()
            
            print '\n'            
            overall.append(trial)
            best_seq_perstate_permutnum[str(chosen_state1)+'_'+str(chosen_state2)] = []
            for muts in range(L-mutation_threshold+1):
                best_seq_perstate_permutnum[str(chosen_state1)+'_'+str(chosen_state2)].append(min(objective_func_sequence_dict['Mutations_'+str(muts)]))
        

            
            if self.display_plot and plot_separate:
                import matplotlib.pyplot as plt
                listt = []
                for i in range(len(best_seq_perstate_permutnum[str(chosen_state1)+'_'+str(chosen_state2)])):
                    listt.append([i,best_seq_perstate_permutnum[str(chosen_state1)+'_'+str(chosen_state2)][i][0]])
                states=[x[0] for x in listt]
                probs=[x[1] for x in listt]
                plt.scatter(states,probs)
                plt.title('Two State Optimization for State'+str(state_num1)+' and State'+str(state_num2))
                plt.xlabel('Mutation Number')
                plt.ylabel('Objective Function Value')
                if pdf_output_plot:                
                    pp.savefig()
                #plt.show()                  
                plt.close()
                
                
            
            if write_opt_file_name:
                for i in range(len(best_seq_perstate_permutnum[str(chosen_state1)+'_'+str(chosen_state2)])):
                    open(opt_filename,'a').write(str(best_seq_perstate_permutnum[str(chosen_state1)+'_'+str(chosen_state2)][i][0])+','+str(best_seq_perstate_permutnum[str(chosen_state1)+'_'+str(chosen_state2)][i][1])+','+str(chosen_state1)+'_probability='+str(best_seq_perstate_permutnum[str(chosen_state1)+'_'+str(chosen_state2)][i][2])+','+str(chosen_state2)+'_probability='+str(best_seq_perstate_permutnum[str(chosen_state1)+'_'+str(chosen_state2)][i][3])+','+str('num_Mutations=')+str(i)+'\n')
        
        if pdf_output_plot:
            pp.close() 
         
        if self.display_plot and plot_together:
            import matplotlib.pyplot as plt
            for stuff in best_seq_perstate_permutnum:
                listt = []
                for i in range(len(best_seq_perstate_permutnum[stuff])):
                    listt.append([i,best_seq_perstate_permutnum[stuff][i][0]])
                states=[x[0] for x in listt]
                probs=[x[1] for x in listt]
                plt.scatter(states,probs)
                plt.xlabel('Mutation Number')
                plt.ylabel('Objective Function Value')
                plt.title('Two State Optimization for All Defined Pairs of States')
        

        self.best_results=best_seq_perstate_permutnum
        
        print '\n' 

        
        
        
        elapsedtime = time.time() - starttime 
        print('time = %0.2f' % elapsedtime)

  
        

