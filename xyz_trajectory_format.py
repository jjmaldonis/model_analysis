import sys
import os
from model import Model
from rot_3d import rot


def main():
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

    print("Atom species transform:")
    m = Model(modelfiles[0])
    types = []
    for atom in m.atoms:
        if atom.z not in types:
            types.append(atom.z)
    types.sort()
    types_dict = {}
    i = 1
    for type in types:
        if type not in types_dict:
            types_dict[type] = i
            i += 1
    print types_dict

    print('Getting trajectories...')
    of = open('traj.xyz','w')

    rot_arr = [ 0.984104, 0.007169, -0.177450, 0.008311, 0.996231, 0.086340, 0.177400, -0.086442, 0.980335]
    for model in modelfiles:
        print(model)
        m = Model(model)
        rot(m,rot_arr)
        of.write('{0}\nAtoms\n'.format(m.natoms))
        for atom in m.atoms:
            of.write('{0} {1} {2} {3}\n'.format(types_dict[atom.z], atom.coord[0], atom.coord[1], atom.coord[2]))

    of.close()

    print('Output filename is traj.xyz')


if __name__ == "__main__":
    main()
