import sys
from model import  Model
import numpy as np
from recenter_model import recenter_model
from atom import Atom
import copy

def periodic_continue_model(m,ex):
    atoms = m.atoms
    pa = np.zeros(m.natoms*((ex)**3),dtype = Atom)
    n = 0
    for i in range(0,ex):
        #print('i={0}'.format(i))
        for j in range(0,ex):
            #print('j={0}'.format(j))
            #print(range(0,ex+1))
            for k in range(0,ex):
                #print('k={0}'.format(k))
                ta = [i*m.lx,j*m.ly,k*m.lz]
                for l in xrange(m.natoms):
                    x = m.atoms[l].coord[0] + ta[0]
                    y = m.atoms[l].coord[1] + ta[1]
                    z = m.atoms[l].coord[2] + ta[2]
                    pa[n] = Atom(n, m.atoms[l].z, x, y, z)
                    #print(pa[n])
                    n += 1
    mp = Model(m.comment, m.lx*ex, m.ly*ex, m.lz*ex, pa)
    recenter_model(mp)
    return mp


def main():
    modelfile = sys.argv[1]
    expansion = int(sys.argv[2])
    m = Model(modelfile)
    mp = periodic_continue_model(m,expansion)
    #mp.write_real_xyz()
    mp.write_our_xyz()


if __name__ == "__main__":
    main()
