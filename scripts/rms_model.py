
import os
import sys
from model import Model
import math


def rms_closest(m1,m2):
    if m1.natoms != m2.natoms: raise Exception("Error! The two models don't have the same number of atoms!")

    r = 0.0
    for atomi in m1.atoms:
        #m2.add_atom(atomi) # Need to add atomi to m2 for next function to work. Actually we don't!
        closest = m2.nearest_neigh_of_same_type(atomi) # This works correctly
        d = m1.dist(atomi,closest)
        r += d**2
        #m2.remove_atom(atomi)
    r = math.sqrt(r/float(m1.natoms))
    return r


def rms(m1,m2):
    if m1.natoms != m2.natoms: raise Exception("Error! The two models don't have the same number of atoms!")

    # This won't always be true, but assume atom 1 in m1
    # is atom 1 in m2.
    r = 0.0
    for i,atomi in enumerate(m1.atoms):
        d = m1.dist(atomi,m2.atoms[i])
        r += d**2
    r = math.sqrt(r/float(m1.natoms))
    return r


def find_nearest_atom(atom,m):
    #return m.nearest_neigh(atom) # Not the same models! Actually that might not matter.
    mind = 100000.0
    for atomi in m.atoms:
        if(atomi != atom):# and atom.z == atomi.z):
            d = m.dist(atom,atomi)
            if(d < mind):
                mind = d
                closest_atom = atomi
    return closest_atom


def main():
    m1 = Model(sys.argv[1])
    m2 = Model(sys.argv[2])
    print("RMS {0}".format(rms_closest(m1,m2)))
    return

    models_prefix = sys.argv[1]

    print("Collecting paths and models...")
    if( '/' in models_prefix):
        path = models_prefix[:models_prefix.rfind('/')+1]
    else:
        path = './'
    files = os.listdir(path)
    modelfiles = []
    for file in files:
        if models_prefix[models_prefix.rfind('/')+1:] in file:
            modelfiles.append(path+file)
    print("Sorting by filename...")
    modelfiles.sort()
    
    m1 = Model(modelfiles[0])

    print("Running rms analysis on all models.")
    i = 1
    for model in modelfiles:
        m2 = Model(model)
        print("{0} {1}".format(model.strip().split()[-1][:-4],rms_closest(m1,m2)))
        i += 1
        



if __name__ == "__main__":
    main()
