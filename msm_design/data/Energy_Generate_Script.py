# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 23:14:41 2017

@author: Trent
"""

import msm_model
import random

#Define the number of mutations you want
num_mutations=35
#Define the total number of energies you want to generate
num_energies=250
#Name the output file
seqEfilename = 'Test_Writing_Energy_File.txt'
#Path to sequence alignment file, if no file change seq_align_file=False
seq_align_file = '/home/thinkle2@ad.wisc.edu/mod_test/data/smallalignment.txt'
#Path to pdbfile
pdbfile='home/thinkle2@ad.wisc.edu/mod_test/data/smallalignment.txt'

wildtype = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
L=len(wildtype)
wtlist = [wildtype[i] for i in range(len(wildtype))]

#Generate the solution space
overall_solution_space = msm_model.import_sequence_alignment(
                        seq_align_file=seq_align_file,
                        reduce_alignment=1,seq_len=L)

open(seqEfilename,'w')
#Iterate 'num_energies' times
for j in range(num_energies):
    randAAlist = []
    #Iterate over the maximum possible mutations
    for i in range(num_mutations):
        randAAlist.append(random.randint(0,L-1))
    #Find mutant positions and mutate
    for mutantpositions in randAAlist:
        wtlist[mutantpositions] = random.choice(overall_solution_space[mutantpositions])
    random_seq = ''.join(wtlist)
    #Score the new sequence
    energy=msm_model.score_sequence_rosetta(random_seq,
                                     pdbfile=pdbfile,
                                     rosetta_path = '/home/romeroroot/code/rosetta_src_2016.17.58663_bundle/main/source/bin/',
                                     rosetta_db = '/home/romeroroot/code/rosetta_src_2016.17.58663_bundle/main/database')
                                
    #Show and write the output
    print(random_seq, energy)
    open(seqEfilename,'a').write(random_seq+','+str(energy)+'\n')
