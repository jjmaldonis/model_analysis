""" This file generates a perfect <0,0,12,0> icosahedron.
    It is actually called a dodecahedron (platonic solid).
    You can specify the bond length between nearest neighbor atoms in the code """
import sys, random, copy
from math import sqrt
from fractions import Fraction
from model import Model
from atom import Atom
from rot_3d import rot as rotate
from rot_3d import calc_rot_array
from model_to_gams_input import convert

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

def perturb_atom(m,i=None,amount=None,axis=None):
    if(i == None):
        i = random.randint(0,m.natoms-1)
    if(axis == None):
        axis = random.randint(0,2)
    elif(axis == 'all'):
        for axis in [0,1,2]:
            perturb_atom(m,i=i,axis=axis)
        return None
    if(amount == None):
        lx = m.lx
        ly = m.ly
        lz = m.lz
        m.lx = 10000
        m.ly = 10000
        m.lz = 10000
        mindist = sorted(m.get_all_dists())
        m.lx = lx
        m.ly = ly
        m.lz = lz
        mindist = mindist[0][-1]
        frac = 10
        amount = random.uniform(-mindist/frac,mindist/frac)
    c = list(m.atoms[i].coord)
    c[axis] = c[axis] + amount
    m.atoms[i].coord = tuple(c)
    #print("Perturbed atom {0} by {1} in direction {2}".format(i,amount,axis))
    return("Perturbed atom {0} by {1} in direction {2}".format(i,amount,axis))

def swap_atom(m,i=None,j=None):
    if(i == None): i = random.randint(0,m.natoms-1)
    if(j == None): j = random.randint(0,m.natoms-1)
    while(i == j): j = random.randint(0,m.natoms-1)
    atom = copy.deepcopy(m.atoms[i])
    m.atoms[i] = m.atoms[j]
    m.atoms[j] = atom
    m.atoms[i].id = i
    m.atoms[j].id = j
    #print("Swapped atoms {0} and {1}".format(i,j))
    return("Swapped atoms {0} and {1}".format(i,j))

def create_basic_models():
    lattice_param = 1
    a = 30
    b = 135
    c = 45
    rot_arr = calc_rot_array(a,b,c,deg=True)

    perfect = icosahedron(lattice_param,save=True,filename='icosahedron.perfect.xyz')
    convert(perfect,'prototype','icosahedron.perfect.txt')

    sm = copy.deepcopy(perfect)
    swap_atom(sm,2,3)
    sm.write_real_xyz('icosahedron.perfect.swap.xyz')
    convert(sm,'polyhedron','icosahedron.perfect.swap.txt')

    rm = copy.deepcopy(perfect)
    rotate(rm,rot_arr)
    rm.write_real_xyz('icosahedron.perfect.rot.xyz')
    convert(rm,'polyhedron','icosahedron.perfect.rot.txt')

    rm = copy.deepcopy(sm)
    rotate(rm,rot_arr)
    rm.write_real_xyz('icosahedron.perfect.swap.rot.xyz')
    convert(rm,'polyhedron','icosahedron.perfect.swap.rot.txt')

    imperfect = copy.deepcopy(perfect)
    perturb_atom(imperfect,1,0.1,0)
    imperfect.write_real_xyz('icosahedron.imperfect.xyz')
    convert(imperfect,'prototype','icosahedron.imperfect.txt')

    sim = copy.deepcopy(imperfect)
    swap_atom(sim,2,3)
    sim.write_real_xyz('icosahedron.imperfect.swap.xyz')
    convert(sim,'polyhedron','icosahedron.imperfect.swap.txt')

    rim = copy.deepcopy(imperfect)
    rotate(rim,rot_arr)
    rim.write_real_xyz('icosahedron.imperfect.rot.xyz')
    convert(rim,'polyhedron','icosahedron.imperfect.rot.txt')

    rim = copy.deepcopy(sim)
    rotate(rim,rot_arr)
    rim.write_real_xyz('icosahedron.imperfect.swap.rot.xyz')
    convert(rim,'polyhedron','icosahedron.imperfect.swap.rot.txt')

def create_randomized_model(num,dir=dir):
    lattice_param = 1
    perfect = icosahedron(lattice_param,save=False,filename='icosahedron.perfect.xyz')
    for i in range(num):
        print(i)
        a = random.uniform(0,360)
        b = random.uniform(0,360)
        c = random.uniform(0,360)
        #a = random.randrange(0,345,15)
        #b = random.randrange(0,345,15)
        #c = random.randrange(0,345,15)
        rot_arr = calc_rot_array(a,b,c,deg=True)

        imperfect = copy.deepcopy(perfect)
        perturb_atom(imperfect,axis='all')
        #swap_atom(imperfect)
        imperfect.write_real_xyz(dir+'icosahedron.{0}.random.xyz'.format(i))
        convert(imperfect,'polyhedron',dir+'icosahedron.{0}.random.txt'.format(i))
        rotate(imperfect,rot_arr)
        imperfect.write_real_xyz(dir+'icosahedron.{0}.random.rot.xyz'.format(i))
        convert(imperfect,'polyhedron',dir+'icosahedron.{0}.random.rot.txt'.format(i))

def main():
    #perfect = Model('icosahedron.perfect.xyz')
    perfect = icosahedron(2.85)
    print(perfect)
    #imperfect = Model('icosahedron.imperfect.xyz')
    #create_basic_models()
    #create_randomized_model(1000,dir='random_models/')
    

if __name__ == '__main__':
    main()
