
import sys
import math
import numpy as np
from model import Model

def rxrydist(atomi,atomj,m):
    x = (atomj.coord[0] - atomi.coord[0])
    y = (atomj.coord[1] - atomi.coord[1])
    while(x > m.lx/2): x = -m.lx + x
    while(y > m.ly/2): y = -m.ly + y
    while(x < -m.lx/2): x = m.lx + x
    while(y < -m.ly/2): y = m.ly + y
    return x,y


def rdf_2d(m,dr):
    # dr is the bin size in atom coord units (A probably)
    # it should be large enough so that the intensity isn't 1 in every bin
    # but small enough to not overlook important information (ie peaks)

    # Algorithm:
    # Calculation the distance between all pairs of atoms, and calculate
    # the rx and ry components of them. This will results in two huge
    # lists: rx and ry. Then create a square matrix of size boxlen/dr,
    # and go through rx and ry simultaneously, incrementing the matrix
    # index if an rx,ry pair falls in that spot. Contour plot the matrix.

    rx = []
    ry = []

    #print("Calculating all distances...")
    for i,atomi in enumerate(m.atoms):
        for j,atomj in enumerate(m.atoms):
            if(j != i):
                x,y = rxrydist(atomi,atomj,m)
                rx.append(x)
                ry.append(y)

    size = int(np.ceil(m.lx/dr))

    # This does by hand what the np.histogram does automattically.
    #print("Initializing matrix: {0}x{1}".format(size,size))
    #mat = np.zeros((size,size))
    #print(mat)
    #print("Making histogram: {0} {1}".format(len(rx),len(ry)))
    #for i in range(0,len(rx)):
    #    x = int(math.floor((rx[i]+m.lx/2.0)/dr))
    #    y = int(math.floor((ry[i]+m.lx/2.0)/dr))
    #    mat[x][y] += 1
    #print("Printing final matrix")
    #print(mat.tolist())

    hist2d = np.histogram2d(rx,ry,size)
    print(hist2d[0].tolist())


def main():
    m = Model(sys.argv[1])
    rdf_2d(m,0.1)

if __name__ == '__main__':
    main()

