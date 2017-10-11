# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 22:52:07 2017

@author: Trent
"""

from setuptools import setup

setup(
    name='msm_design',
    version='0.9',
    description='Markov State Modeling Protein Design Package utilizing Machine Learning and Brute Force Optimization of a specified conformational state of a protein',
    author='Trent Hinkle',
    author_email='trenth12@gmail.com',
    url='https://github.com/thinkle12/msm_design/',
    packages = ['msm_design'],
    package_data = {
    # If any package contains *.txt or *.rst files, include them:
    'msm_design/data': ['Msm_Design_Runner.py','smallalignment.txt',
                      'ubiquitin_State0_energies_amber_4positionmutations_top4AA.txt',
                      'ubiquitin_State0_energies_rosetta_4positionmutations_top4AA.txt'],
        'msm_design/data/model_files':['w_nopairnoweight_State0.txt','w_nopairnoweight_State1.txt',
                                      'w_nopairnoweight_State2.txt','w_nopairnoweight_State3.txt',
                                      'w_nopairnoweight_State4.txt','w_nopairnoweight_State5.txt,
                                      'w_nopairnoweight_State6.txt','w_nopairnoweight_State7.txt',
                                       'w_nopairnoweight_State8.txt','w_nopairnoweight_State9.txt'],
        'msm_design/data/pdb_files':['State0.pdb','State1.pdb','State2.pdb','State3.pdb','State4.pdb,
                                    'State5.pdb','State6.pdb','State7.pdb','State8.pdb','State9.pdb]
        }
    )
