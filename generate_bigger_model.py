

import sys
from model import Model
from atom import Atom


def glm(modelfile,mult):
    model = Model(modelfile)

    outfile = modelfile+'.bigger'

    model = Model(modelfile)
    atoms = model.atoms[:]
    natoms = model.natoms
    for atom in atoms:
        for i in range(0,mult):
            for j in range(0,mult):
                for k in range(0,mult):
                    if(i == j == k == 0): continue
                    natoms += 1
                    x = atom.coord[0] + i*model.lx
                    y = atom.coord[1] + j*model.ly
                    z = atom.coord[2] + k*model.lz
                    model.add_atom(Atom(natoms,atom.z,x,y,z))

    # Shift right 1/2 world and left mult/2 worlds
    # which is equivalent to left (mult-1)/2 worlds
    for i,atom in enumerate(model.atoms):
        model.atoms[i].set_coord(model.atoms[i].coord[0] - (mult-1)/2.0*model.lx, model.atoms[i].coord[1] - (mult-1)/2.0*model.ly, model.atoms[i].coord[2] - (mult-1)/2.0*model.lz)

    model.lx *= mult
    model.ly *= mult
    model.lz *= mult

    model.write_our_xyz(outfile)


def main():
    modelfile = sys.argv[1]
    size_multiplier = int(sys.argv[2])
    glm(modelfile,size_multiplier)



if __name__ == "__main__":
    main()
