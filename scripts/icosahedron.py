""" This file generates a perfect <0,0,12,0> icosahedron.
    It is actually called a dodecahedron (platonic solid).
    You can specify the bond length between nearest neighbor atoms in the code """
import sys, random, copy
from math import sqrt
from fractions import Fraction
from model import Model
from atom import Atom

def dodecahedron(b,save=False,filename=None):
    # http://en.wikipedia.org/wiki/Dodecahedron#Regular_dodecahedron
    if(filename == None): filename = 'dodecahedron.xyz'
    # b is the nearest-neighbor interatomic bond distance
    p = 1.61803398875 # golden ratio
    b = b * 0.5 * p
    coords = [[0,0,0]]

    coords.append([p,0,1./p])
    coords.append([-p,0,1./p])
    coords.append([-p,0,-1./p])
    coords.append([p,0,-1./p])

    coords.append([1./p, p, 0])
    coords.append([1./p, -p, 0])
    coords.append([-1./p, -p, 0])
    coords.append([-1./p, p, 0])

    coords.append([0, 1./p, p])
    coords.append([0, 1./p, -p])
    coords.append([0, -1./p, -p])
    coords.append([0, -1./p, p])

    coords.append([1,1,1])
    coords.append([1,-1,1])
    coords.append([-1,-1,1])
    coords.append([-1,1,1])

    coords.append([-1,1,-1])
    coords.append([1,1,-1])
    coords.append([1,-1,-1])
    coords.append([-1,-1,-1])

    coords = [ [b*x for x in c] for c in coords]
    m = 0
    for coord in coords:
        for x in coord:
            if(abs(x) > m): m = abs(x)

    atoms = [Atom(i,14,c[0],c[1],c[2]) for i,c in enumerate(coords)]
    model = Model('dodecahedron',m,m,m,atoms)
    model.filename = filename

    if(save):
        model.write_real_xyz(model.filename)
        #f = open(filename,'w')
        #f.write(str(len(coords))+'\n')
        #f.write('{0} {0} {0} comment\n'.format(m))
        #for c in coords:
        #    f.write('Si ' + ' '.join([str(x) for x in c]) + '\n')
        #f.close()
    return model

def icosahedron(b,save=False,filename=None):
    # http://en.wikipedia.org/wiki/Regular_icosahedron
    if(filename == None): filename = 'icosahedron.xyz'
    # b is the nearest-neighbor interatomic bond distance
    b = b * sqrt(2)* 1.113516364 / 1.84375
    coords = [atom.coord for atom in dodecahedron(b).atoms]
    vertices = range(1,21)
    faces = [
        [15, 2, 12, 16, 9],
        [11, 19, 10, 4, 18],
        [6, 7, 14, 15, 12],
        [10, 18, 17, 5, 8],
        [9, 16, 13, 8, 5],
        [7, 6, 20, 19, 11],
        [5, 1, 13, 4, 18],
        [4, 19, 6, 14, 1],
        [1, 9, 13, 12, 14],
        [2, 3, 20, 7, 15],
        [2, 16, 8, 17, 3],
        [3, 20, 11, 10, 17]]
    face_coords = [(0,0,0)]
    for f in faces:
        center = (sum(coords[i][0] for i in f)/5., sum(coords[i][1] for i in f)/5., sum(coords[i][2] for i in f)/5.)
        face_coords.append(center)
    m = 0
    for coord in face_coords:
        for x in coord:
            if(abs(x) > m): m = abs(x)

    atoms = [Atom(i,14,c[0],c[1],c[2]) for i,c in enumerate(face_coords)]
    model = Model('icosahedron',m,m,m,atoms)
    model.filename = filename
    if(save):
        model.write_real_xyz(model.filename)
        #f = open(filename,'w')
        #f.write(str(len(face_coords))+'\n')
        #f.write('{0} {0} {0} comment\n'.format(m))
        #for c in face_coords:
        #    f.write('Si ' + ' '.join([str(x) for x in c]) + '\n')
        #f.close()
    return model

def perturb_atom(m):
    i = random.randint(0,m.natoms-1)
    mindist = sorted(m.get_all_dists())[0][-1]
    amount = random.uniform(-mindist/10,mindist/10)
    axis = random.randint(0,2)
    c = list(m.atoms[i].coord)
    c[axis] = c[axis] + amount
    m.atoms[i].coord = tuple(c)
    #print("Perturbed atom {0} by {1} in direction {2}".format(i,amount,axis))
    return m

def swap_atom(m):
    i = random.randint(0,m.natoms-1)
    j = random.randint(0,m.natoms-1)
    while(i == j): j = random.randint(0,m.natoms-1)
    atom = copy.deepcopy(m.atoms[i])
    m.atoms[i] = m.atoms[j]
    m.atoms[j] = atom
    m.atoms[i].id = i
    m.atoms[j].id = j
    #print("Swapped atoms {0} and {1}".format(i,j))
    return m

def main():
    lattice_param = 1
    perfect = icosahedron(lattice_param,save=True,filename='gams/icosahedron.xyz')
    #for atom in perfect.atoms:
    #    print(atom.coord)
    #print('')
    for i in range(10):
        m = copy.deepcopy(perfect)
        perturb_atom(m)
        swap_atom(m)
        #for atom in m.atoms:
        #    print(atom.coord)
        #print('')


if __name__ == '__main__':
    main()
