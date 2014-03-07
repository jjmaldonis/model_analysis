

import sys
from model import Model
import math


def rms(m1,m2):
    if m1.natoms != m2.natoms: raise Exception("Error! The two models don't have the same number of atoms!")

    # This won't always be true, but assume atom 1 in m1
    # is the closest atom to atom 1 in m2.
    r = 0.0
    for atomi in m1.atoms:
        closest = find_nearest_atom(atomi,m2)
        d = m1.dist(atomi,closest)
        r += d**2
    r = math.sqrt(r/float(m1.natoms))
    return r


def find_nearest_atom(atom,m):
    mind = 100000.0
    for atomi in m.atoms:
        d = m.dist(atom,atomi)
        if(d < mind):
            mind = d
            closest_atom = atomi
    return closest_atom


def main():
    m1 = Model(sys.argv[1])
    m2 = Model(sys.argv[2])

    print(rms(m1,m2))


if __name__ == "__main__":
    main()
