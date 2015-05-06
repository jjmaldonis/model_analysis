""" This file generates a perfect <0,0,12,0> icosahedron.
    It is actually called a dodecahedron (platonic solid).
    You can specify the bond length between nearest neighbor atoms in the code """
import sys
from math import sqrt
from fractions import Fraction

def dodecahedron(b,save=False):
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
    if(save):
        m = 0
        for coord in coords:
            for x in coord:
                if(abs(x) > m): m = abs(x)
        f = open('dodecahedron.xyz','w')
        f.write(str(len(coords))+'\n')
        f.write('{0} {0} {0} comment\n'.format(m))
        for c in coords:
            f.write('Si ' + ' '.join([str(x) for x in c]) + '\n')
        f.close()
    return coords

def icosahedron(b,save=False):
    # b is the nearest-neighbor interatomic bond distance
    b = b * sqrt(2)* 1.113516364 / 1.84375
    coords = dodecahedron(b)
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
    if(save):
        m = 0
        for coord in face_coords:
            for x in coord:
                if(abs(x) > m): m = abs(x)
        f = open('icosahedron.xyz','w')
        f.write(str(len(face_coords))+'\n')
        f.write('{0} {0} {0} comment\n'.format(m))
        for c in face_coords:
            f.write('Si ' + ' '.join([str(x) for x in c]) + '\n')
        f.close()
    return face_coords


def main():
    lattice_param = 2.
    icosahedron(lattice_param,save=True)

if __name__ == '__main__':
    main()
