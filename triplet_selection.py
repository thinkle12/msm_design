# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 21:37:14 2016

@author: Trent
"""

import numpy
import itertools

"Here we select the maximum distance for a pair of AAs to be considered important"
angstrom_distance = 2.3
#########################################################################
#this should really use the 'state' pdb files....
L = 76
pdbfile = 'State0.pdb'
AAs = ['A', 'C', 'D', 'E', 'F', 'G', 'H','I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W','Y'] #all 20 amino acids



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

tridistances = {}
for i in range(1,L+1):
    for j in range(1,L+1):
        for k in range(1,L+1):
            if i<j and j<k:
                point1 = j
                try:
                    minimumdistances['minimumdistances'+str(j)+'-'+str(k)]
                    point2 = k
                except KeyError:
                    pass
                tridistance = minimumdistances['minimumdistances'+str(j)+'-'+str(k)]
                tridistances['tridist'+str(i)+'-'+str(j)+'-'+str(k)] = []
                tridistances['tridist'+str(i)+'-'+str(j)+'-'+str(k)] = [minimumdistances['minimumdistances'+str(i)+'-'+str(j)],minimumdistances['minimumdistances'+str(i)+'-'+str(k)],minimumdistances['minimumdistances'+str(j)+'-'+str(k)]]

"Here we determine which pairs are close enough by setting an angstrom_distance constraint"
closetripletatoms = []
for i in range(1,L+1):
    for j in range(1,L+1):
        for k in range(1,L+1):
            if all([i<j and j<k and all(n<float(angstrom_distance) for n in tridistances['tridist'+str(i)+'-'+str(j)+'-'+str(k)])]):
                closetripletatoms.append((i,j,k))

"Here we subtract of (1,1) from each pair to match the indecies of the regression format"  

file1='close_AA_triplets_ubiquitin.txt'
open(file1,'w')
          
newclosetripletatoms = numpy.array(closetripletatoms) - (1,1,1)
newclosetripletatoms = map(tuple, newclosetripletatoms)

yyy = []
for triplets in newclosetripletatoms:
    for d in range(len(triplets)):
        open(file1,'a').write(str(triplets[d])+',')
    open(file1,'a').write('\n')

print newclosetripletatoms
print str('Angstrom Distance = ')+str(angstrom_distance)
