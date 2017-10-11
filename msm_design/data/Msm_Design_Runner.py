# -*- coding: utf-8 -*-
"""
Created on Thu Oct 05 22:09:36 2017

@author: Trent
"""

#import matplotlib
#matplotlib.use('Agg')  

from msm_design import msm_model
from msm_design import msm_optimization


#trial = msm_model.SGD_Static(pdbfile='C:/Users/Trent/Desktop/PythonFiles/Project/Package/State0.pdb',
#            sequence_energy_file='C:/Users/Trent/Desktop/PythonFiles/Project/ubiquitin_State0_energies_amber_test.txt',
#            sequence_alignment_file='C:/Users/Trent/Desktop/PythonFiles/Project/smallalignment.txt',
#            wt_seq='MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG',
#            gamma_multiplier=1,regular_linear_regression=True,lasso=False,
#            ridge=False,ridge_coef=.0001,lambda_lasso_coef=.0001,pair_select_list=False,pair_dist=3)
##
#trial.static_model(model_output_text_file=False)

#trial = SGD_Online(pdbfile='/home/thinkle2@ad.wisc.edu/ubiquitin_structures/State76.pdb',sequence_alignment_file='smallalignment.txt',
#                   wt_seq='MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG',
#                  pair_select_list=True,pair_dist=2,output_energy_file=True)
                  
#print(score_sequence_rosetta(pdbfile='/home/thinkle2@ad.wisc.edu/ubiquitin_structures/State6.pdb',newseq='MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'))
#print(score_sequence_rosetta(pdbfile=trial.pdbfile,newseq=trial.wt_seq))
#print(score_sequence_amber(pdbfile='/home/thinkle2@ad.wisc.edu/ubiquitin_structures/State6.pdb',newseq='MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'))
#print(score_sequence_amber(pdbfile=trial.pdbfile,newseq=trial.wt_seq))
#trial.online_model(num_mutations=20,energy_scoring_function='rosetta',max_computation_time=252000,max_training_sets=200,mad_cutoff=2)
#


opt=msm_optimization.Sequence_Optimization(model_file_dir='C:/Users/Trent/Desktop/PythonFiles/Project/w_amber_nopair//',
    wt_seq='MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG',
    sequence_alignment_file='C:/Users/Trent/Desktop/PythonFiles/Project/smallalignment.txt',reduce_alignment=1)
    
So = opt.two_state_design([[7,34],[5,19]],plot_separate=True,num_mutation_pathway_passes=100,
    write_opt_file_name='two_state_output.txt',
    pdf_output_plot='two_state_plot_output.pdf')