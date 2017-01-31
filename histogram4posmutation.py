# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 19:35:31 2016

@author: Trent
"""

h = open('ubiquitin_State0_energies_500000_reg_6angpair_4positionmutations_top4AA_smallalign_State0.txt'
,'r+')
text = h.read()
text = text.split('\n')
#text = pickle.load(open('seq_energies_state0_ubiquitin.p'))
del text[-1]
sequences = [t.split(',')[0] for t in text]
Energies = [float(t.split(',')[1]) for t in text]

import matplotlib.pyplot as plt

plt.hist(Energies, bins=600)
plt.axis([-50,150,0,80000])
plt.show()

