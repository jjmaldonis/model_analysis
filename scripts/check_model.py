

import sys
from model import Model
import znum2sym
from pprint import pprint


#class Checks(object):
#    def __init__(self):
#        super(Checks,self).__init__()

def check_dists(model):
    dists = model.get_all_dists()
    dists = [dist[2] for dist in dists]
    print("Minimimum atomic spacing: {0}\n".format(min(dists)))

def atoms_in_box(model):
    for atom in model.atoms:
        if( atom.coord[0] < -model.lx/2): raise Exception('ERROR! Atom {0} is out of your box!'.format(atom))
        if( atom.coord[0] > model.lx/2): raise Exception('ERROR! Atom {0} is out of your box!'.format(atom))
        if( atom.coord[1] < -model.ly/2): raise Exception('ERROR! Atom {0} is out of your box!'.format(atom))
        if( atom.coord[1] > model.ly/2): raise Exception('ERROR! Atom {0} is out of your box!'.format(atom))
        if( atom.coord[2] < -model.lz/2): raise Exception('ERROR! Atom {0} is out of your box!'.format(atom))
        if( atom.coord[2] > model.lz/2): raise Exception('ERROR! Atom {0} is out of your box!'.format(atom))
    print("\nYour atom coords are all inside the box.\n")

def atom_density(model):
    ad = model.natoms / (model.lx * model.ly * model.lz)
    print("\nAtom density: {0}".format(ad))
    print(' ')

def composition(model):
    d = model.composition()
    for key in d:
        print('{0}: {1}'.format(key,d[key]))
    print(' ')

def positions(model):
    print('\nBox sizes: {0} x {1} x {2}'.format(model.lx, model.ly, model.lz))
    print('Half box sizes: {0}, {1}, {2}'.format(model.lx/2, model.ly/2, model.lz/2))
    xx = [atom.coord[0] for atom in model.atoms]
    yy = [atom.coord[1] for atom in model.atoms]
    zz = [atom.coord[2] for atom in model.atoms]
    print('x-min: {0}\t x-max: {1}'.format(min(xx),max(xx)))
    print('y-min: {0}\t y-max: {1}'.format(min(yy),max(yy)))
    print('z-min: {0}\t z-max: {1}\n'.format(min(zz),max(zz)))

def number_of_atoms(model):
    print("\nNumber of atoms {0}:".format(sum(model.atomtypes.values())))
    pprint(model.atomtypes)
    print

def main():
    modelfile = sys.argv[1]
    m = Model(modelfile)
    number_of_atoms(m)
    composition(m)
    atom_density(m)
    atoms_in_box(m)
    positions(m)
    check_dists(m)

    print("You should also check the RDFs! Partials included!\n")


if __name__ == "__main__":
    main()

