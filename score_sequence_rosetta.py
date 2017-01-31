from sys import argv
from random import choice
from subprocess import Popen
from os import remove



rosetta_path = '/home/romeroroot/code/rosetta_src_2016.17.58663_bundle/main/source/bin/'
rosetta_db = '/home/romeroroot/code/rosetta_src_2016.17.58663_bundle/main/database'



def rand_tag():
    '''creates a unique tag that can be used to identify the current files a script is working with.  necessary for running the same program multiple times in the same directory'''
    alpha = 'abcdefghijklmnopqrstuvwxyz0123456789'
    tag = ''.join([choice(alpha) for i in range(15)])
    return tag



def fast_thread_repack_score(pdbfile,newseq,savepdb=False):
    """this function inputs a pdbfile name and a sequence, mutates the pdb file, repacks, and returns the scores as a dict
    REQUIRES: pdbfile to be consecutivly numbered from 1 to N, and len(newseq) = N. Used when we want to score seqence variants many times
    """

    tag = rand_tag() # generate a random tag to associate with this run


    # generate a resfile to mutate the pdb
    resfile = """
    #header
    NATAA # keep this so all sites aren't designed
    USE_INPUT_SC # this to also consider the WT rotamer, if we're mutating WT to WT

    start

    #body\n"""

    for i in range(len(newseq)):
        resfile += '%i    A    PIKAA %s \n'%(i+1,newseq[i])

    open('resfile_'+tag+'.res','w').write(resfile)


    # the options for the run
    options = [rosetta_path+'fixbb.linuxgccrelease', # fixbb is the program used for threading a repacking
               '-s '+pdbfile, 
               '-database '+rosetta_db, # good to explicitly define location of rosetta DB
               '-resfile resfile_'+tag+'.res', # a resfile to specify mutations
               '-out:suffix _'+tag] # this adds a suffix to the score output file: score_[RNDTAG].sc

    # run fixbb
    Popen(options,stdout=open('/dev/null','w')).wait() # send stdout to the trash, report stderr


    # read the output
    scorefile = open('score_%s.sc'%tag).read().split('\n')
    names = scorefile[1].split(':')[1].split()
    values = scorefile[2].split(':')[1].split()
    score = dict((names[i],float(values[i])) for i in range(len(names)-1)) # use i-1 because the last entry is the filename


    #delete the evidence
    remove('score_%s.sc'%tag)
    remove('resfile_'+tag+'.res')

    pdbname = (pdbfile[:-4]+'_'+tag+'_0001.pdb').split('/')[-1]
    if savepdb:
        return score,pdbname
    else:
        remove(pdbname)
        return score

pdbfile = '/home/romeroroot/temp_promero2/projects/trent/ubiquitin_structures_from_greg/State0.pdb'
newseq = argv[1]
score = fast_thread_repack_score(pdbfile,newseq)
print score['total_score']


