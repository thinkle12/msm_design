# -*- coding: utf-8 -*-
"""
Created on Thu Oct 05 22:09:36 2017

@author: Trent
"""

import matplotlib
matplotlib.use('Agg')  

#from msm_design import msm_model
#from msm_design import msm_optimization
import msm_model
import msm_optimization


#We start with a simple SGD_Static model from msm_model
#Define the pdbfile with its path
pdbfilepath = '/home/thinkle2@ad.wisc.edu/mod_test/data/pdb_files/State0.pdb'
#Define the sequence energy file
#For information on how to generate the sequence energy file see 'Energy_Generate_script.py' in 'data'
sequence_energy_file_path = '/home/thinkle2@ad.wisc.edu/mod_test/data/ubiquitin_State0_energies_rosetta_4positionmutations_top4AA.txt'
#path to sequence alignment file, If none let sequence_alignment_file=False
sequence_alignment_file_path = '/home/thinkle2@ad.wisc.edu/mod_test/data/smallalignment.txt'
#sequence_alignment_file_path=False

wildtypeseq = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'

#Other parameters are variable... if pair_select_list=False then pair_dist does nothing
#standard parameters are shown below

trial = msm_model.SGD_Static(pdbfile=pdbfilepath,
            sequence_energy_file=sequence_energy_file_path,
            reduce_alignment=2,
	    sequence_alignment_file=sequence_alignment_file_path,
            wt_seq=wildtypeseq,
            gamma_multiplier=1,regular_linear_regression=True,lasso=False,
            ridge=False,ridge_coef=.0001,lambda_lasso_coef=.0001,pair_select_list=False,pair_dist=3,custom_tag='Test_Run')

#Model trains here using the data stored in trial
#Model will output a text file if model_output_text_file=True with the custom_tag from trial embeded in the name
#If no custom tag set the algorithm will attempt to mimic the name of the given pdb file
trial.static_model(model_output_text_file=True)


#Next we train an online model which will train the model at the same time the energies are being generated
#Again define pdbfile, sequence_alignment_file and wt_seq
#A custom tag is good here to keep track of what models we are training
newpdbfile = '/home/thinkle2@ad.wisc.edu/ubiquitin_structures/State76.pdb'
new_seq_align='/home/thinkle2@ad.wisc.edu/mod_test/data/smallalignment.txt'
wt='MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'

#Here you need to put the path to the rosetta function fixbb.linuxgccrelease
#This will be different for you so make sure you can locate this on your machine
#Without this path being proper the code will not execute
path_to_rosetta_fixbb='/home/romeroroot/code/rosetta_src_2016.17.58663_bundle/main/source/bin/'
#We also need the path to the rosetta database
path_to_rosetta_db='/home/romeroroot/code/rosetta_src_2016.17.58663_bundle/main/database'
trial2 = msm_model.SGD_Online(pdbfile=newpdbfile,
                sequence_alignment_file=new_seq_align,
                reduce_alignment=2,
                wt_seq=wt,
                pair_select_list=False,pair_dist=2,output_energy_file=True,
                custom_tag='Online_Test_20muts_State76',
                rosetta_score_path=path_to_rosetta_fixbb,
                rosetta_database_path=path_to_rosetta_db)

#trial2 will store some of the pre model training information
#Next we run an online model and select important parameters including:
#The number of mutations we want to attempt every time we sample a sequence
#The energy scoring function which will either be rosetta or amber
#The max computation time - can set to a high number if wanted
#The max number of training sets
#The MAD cutoff (Model stops if the MAD is below the cutoff for 100 straight model updates)
trial2.online_model(num_mutations=20,energy_scoring_function='rosetta',max_computation_time=252000,max_training_sets=200,mad_cutoff=2)
#running this will automatically output a text file of the model parameters in the same directory as the script


opt=msm_optimization.Sequence_Optimization(model_file_dir='/home/thinkle2@ad.wisc.edu/mod_test/data/w_files//',
    wt_seq='MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG',
    sequence_alignment_file='/home/thinkle2@ad.wisc.edu/mod_test/data/smallalignment.txt',reduce_alignment=2)
    

iterate = [2, 4, 5, 15, 17, 24, 27, 30, 35, 41, 42, 47]
ll = []
for i in iterate:
    for j in iterate:
        if i<j:
            ll.append([i,j])


So = opt.two_state_design(ll,plot_separate=True,num_mutation_pathway_passes=10000,
    write_opt_file_name='two_state_trial.txt',
    pdf_output_plot='two_state_trial.pdf')
