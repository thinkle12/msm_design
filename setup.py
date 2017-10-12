# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 22:52:07 2017

@author: Trent
"""

from setuptools import setup

setup(
    name='msm_design',
    version='0.12',
    description='Markov State Modeling Protein Design Package utilizing Machine Learning and Brute Force Optimization of a specified conformational state of a protein',
    author='Trent Hinkle',
    author_email='trenth12@gmail.com',
    url='https://github.com/thinkle12/msm_design/',
    packages = ['msm_design'],
    include_package_data=True,
    package_data = {
    # If any package contains *.txt or *.rst files, include them:
    'msm_design/data': ['msm_design/data/Msm_Design_Runner.py','msm_design/data/smallalignment.txt',
                      'msm_design/data/ubiquitin_State0_energies_amber_4positionmutations_top4AA.txt',
                      'msm_design/data/ubiquitin_State0_energies_rosetta_4positionmutations_top4AA.txt'],
        'msm_design/data/model_files':['msm_design/data/model_files/w_nopairnoweight_State0.txt','msm_design/data/model_files/w_nopairnoweight_State1.txt',
                                      'msm_design/data/model_files/w_nopairnoweight_State2.txt','msm_design/data/model_files/w_nopairnoweight_State3.txt',
                                      'msm_design/data/model_files/w_nopairnoweight_State4.txt','msm_design/data/model_files/w_nopairnoweight_State5.txt',
                                      'msm_design/data/model_files/w_nopairnoweight_State6.txt','msm_design/data/model_files/w_nopairnoweight_State7.txt',
                                       'msm_design/data/model_files/w_nopairnoweight_State8.txt','msm_design/data/model_files/w_nopairnoweight_State9.txt'],
        'msm_design/data/pdb_files':['msm_design/data/pdb_files/State0.pdb','msm_design/data/pdb_files/State1.pdb',
                                     'msm_design/data/pdb_files/State2.pdb','msm_design/data/pdb_files/State3.pdb',
                                     'msm_design/data/pdb_files/State4.pdb','msm_design/data/pdb_files/State5.pdb',
                                     'msm_design/data/pdb_files/State6.pdb','msm_design/data/pdb_files/State7.pdb',
                                     'msm_design/data/pdb_files/State8.pdb','msm_design/data/pdb_files/State9.pdb']
        }
    )
