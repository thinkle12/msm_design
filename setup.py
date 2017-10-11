# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 22:52:07 2017

@author: Trent
"""

from setuptools import setup

setup(
    name='msm_design',
    version='0.6',
    description='Markov State Modeling Protein Design Package utilizing Machine Learning and Brute Force Optimization of a specified conformational state of a protein',
    author='Trent Hinkle',
    author_email='trenth12@gmail.com',
    url='https://github.com/thinkle12/msm_design/',
    packages = ['msm_design'],
    package_data = {
    # If any package contains *.txt or *.rst files, include them:
    '': ['*.txt', '*.pdb']
        }
    )
