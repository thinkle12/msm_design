import numpy
import itertools
import os
"Here we select the maximum distance for a pair of AAs to be considered important"
angstrom_distance = 4
#########################################################################
#this should really use the 'state' pdb files....
L = 76

namelist = []
for i in range(100):
    namelist.append('State'+str(i)+'.pdb')

for filename in os.listdir(os.curdir):
    if filename in namelist:
        pdbfile = filename

pdbfile = 'State0.pdb'

AAs = ['A', 'C', 'D', 'E', 'F', 'G', 'H','I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W','Y'] #all 20 amino acids



#WTFile = open('1ubq.fasta.txt').read()
#WTFileLines = WTFile.split('\n')
#WildType = WTFileLines[1]
##Replace reference_seq with the wild type import from pdb fasta txt file
#AAoccuranceslist = []
#for aas in AAs:
#    print(str(aas)+str(    )+str(''.join(WildType).count(aas)))
#    AAoccuranceslist.append((aas,''.join(WildType).count(aas)))
#AAoccuranceslist.sort(key=lambda tup: tup[1])
#AAoccuranceslist.reverse()
#tempsolspace = []
#for i in range(5):
#    tempsolspace.append(AAoccuranceslist[i][0])
#
#solspace = list(tempsolspace)
#overall_solution_space = []
#for j in range(len(WildType)):
#    if WildType[j] not in solspace:
#        tempsolspace.append(WildType[j])
#        solution = list(tempsolspace)
#        overall_solution_space.append(solution)
#        del tempsolspace[-1]
#    if WildType[j] in solspace:
#        overall_solution_space.append(solspace)
    
#Select top 5 most occuring AA's to be able to change....
#Not sure what to do with AA's in WT that arent top 5 AA
#Do we include the WT AA at each position to be able to vary as well?
#Or do we not allow those AA's to change...?
#How do we define sequence space...


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

file1='close_AA_pairs_ubiquitin.txt'
open(file1,'w')
          
newcloseatoms = numpy.array(closeatoms) - (1,1)
newcloseatoms = map(tuple, newcloseatoms)

yyy = []
for i in range(len(newcloseatoms)):
    yyy = str(newcloseatoms[i][0])+','+str(newcloseatoms[i][1])
    open(file1,'a').write(yyy+',')

print newcloseatoms
print str('Angstrom Distance = ')+str(angstrom_distance)
