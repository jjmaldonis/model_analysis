import sys
import os
import gr
import numpy as np
from collections import defaultdict
from model import Model

def allcna(m):
    """ Common neighbor analysis """

    cna = defaultdict(list)
    for atom in m.atoms:
        #print(atom)
        #print(atom.neighs)
        if(atom.cn == 0):
            raise Exception("m.generate_neighbors() must be run before calling allcna")
        for n in atom.neighs:
            type = one_cna(atom,n)
            cna[type] += (atom,n)
    return cna

def one_cna(atom1,atom2):
    cna = [0,0,0]
    # Calculate the number of neighbors they have in common
    # Store the common neighbors in 'common_neighs.keys()'
    common_neighs = defaultdict(list)
    for n in atom1.neighs:
        if n in atom2.neighs:
            common_neighs[n] = []
            cna[0] += 1
    if(cna == 0):
        return cna
    cnl = list(common_neighs)
    #print(atom1,atom2)
    #print(common_neighs)

    for x in common_neighs:
        for y in common_neighs:
            if x in y.neighs:
                common_neighs[x] += [y]
                common_neighs[y] += [x]
    for x in common_neighs:
        common_neighs[x] = list(set(common_neighs[x]))

    #print(common_neighs)
    #return cna

    # Calculate the number of bonds between the common neighbors
    for i,atomi in enumerate(cnl):
        for atomj in cnl[i+1:]:
            if(atomi in atomj.neighs):
                cna[1] += 1

    # Calculate the number of bonds in the longest continuous
    # chain formed by the k bonds between common neighbors
    longestpath = 0
    for cn in common_neighs: #cn = common neighbor
        lpath,lp = temp(common_neighs, common_neighs[cn], [cn], 0, 0, [cn])
        if( lpath[0] in lpath[-1].neighs and lp > 1): lp += 1 # found a cycle
        if(lp > longestpath): longestpath = lp
    #print(longestpath)
    cna[2] = longestpath
    #print(cna)
    #if(cna == [5,4,4]):
    #    print(' ')
    #    #print("WOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO")
    #    print(atom1,atom2)
    #    print(common_neighs)
    #    print( sorted([len(common_neighs[x]) for x in common_neighs]) )
    #    if( sorted([len(common_neighs[x]) for x in common_neighs]) == [0,2,2,2,2]): g = f + 1
    #if(atom1.id == 1244 and atom2.id == 261): g = f+1
    #if( sorted([len(common_neighs[x]) for x in common_neighs]) == [1,1,1,3]): g = f + 1
    return tuple(cna)

def temp(cns,l,v, pd, lp, lpath):
    if( pd > lp):
        lp = pd
        lpath.append(v[-1])
    #print(l,v,pd,lp, lpath)
    for x in l:
        #if(x == v[0]):
        #    pd += 1
        if x not in v:
            lpath, lp = temp(cns,cns[x],v+[x], pd+1, lp, lpath)
            if( pd > lp):
                lp = pd
                lpath.append(v[-1])
    return lpath, lp
            
def cna_print(cna):
    pass


def main():
    #modelfile = sys.argv[1]
    #m = Model(modelfile)
    #m.generate_neighbors(3.5)
    models_prefix = sys.argv[1]
    print("Collecting paths and models...")
    end = models_prefix.rfind('/')+1
    if(end == 0):
        path = './'
    else:
        path = models_prefix[:models_prefix.rfind('/')+1]
    print(path)
    files = os.listdir(path)
    modelfiles = []
    for file in files:
        if models_prefix[models_prefix.rfind('/')+1:] in file:
            modelfiles.append(path+file)
    print("Sorting by filename...")
    modelfiles.sort()
    print('Last model to be counted is {0}'.format(modelfiles[-1]))

    i = 0
    for model in modelfiles:
        if(i%10 == 0):
            print(' ')
            print(model)
            m = Model(model)
            m.generate_neighbors(3.5)
            cnas = allcna(m)
            for x in sorted(cnas.keys()):
                print('{0} {1}'.format(x,len(cnas[x])/2))
        i += 1

if __name__ == "__main__":
    main()
