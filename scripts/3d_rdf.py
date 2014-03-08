
import sys
import math
import numpy as np
from model import Model

def rxryrzdist(atomi,atomj,m):
    x = (atomj.coord[0] - atomi.coord[0])
    y = (atomj.coord[1] - atomi.coord[1])
    z = (atomj.coord[2] - atomi.coord[2])
    while(x > m.lx/2): x = -m.lx + x
    while(y > m.ly/2): y = -m.ly + y
    while(z > m.lz/2): z = -m.lz + z
    while(x < -m.lx/2): x = m.lx + x
    while(y < -m.ly/2): y = m.ly + y
    while(z < -m.lz/2): z = m.lz + z
    return x,y,z


def GetAtomsInBin(m,bx,by,rx,ry,dr):
    xmin = dr*bx - m.lx/2.0
    xmax = dr*(bx+1) - m.lx/2.0
    ymin = dr*by - m.ly/2.0
    ymax = dr*(by+1) - m.ly/2.0
    atoms = []
    for i in range(0,len(rx)):
        if( xmin < rx[i] < xmax and ymin < ry[i] < ymax ):
            atoms.append(int(i/m.natoms))
            atoms.append(i%m.natoms)
    atoms = list(set(atoms)) # Remove dups
    return atoms # List


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

    size = int(np.ceil(m.lx/dr))
    #print("Initializing matrix: {0}x{1}".format(size,size))
    hist3d = np.zeros((size,size,size),dtype=np.int)
    rx = []
    ry = []
    rz = []
    bindict = {}

    #print("Calculating all distances...")
    for i,atomi in enumerate(m.atoms):
        for j,atomj in enumerate(m.atoms):
            if(j != i):
                x,y,z = rxryrzdist(atomi,atomj,m)
                xbin = int(math.floor((x+m.lx/2.0)/dr))
                ybin = int(math.floor((y+m.ly/2.0)/dr))
                zbin = int(math.floor((z+m.lz/2.0)/dr))
                try:
                    bindict[(xbin,ybin,zbin)].append( [i,j] )
                except KeyError:
                    bindict[(xbin,ybin,zbin)] = [[i,j]]
                rx.append(x)
                ry.append(y)
                rz.append(z)
                hist3d[xbin][ybin][zbin] += 1
    print(hist3d.tolist())

    ## This does by hand what the np.histogram does automattically.
    ##print(mat)
    ##print("Making histogram: {0} {1}".format(len(rx),len(ry)))
    ##for i in range(0,len(rx)):
    ##    x = int(math.floor((rx[i]+m.lx/2.0)/dr))
    ##    y = int(math.floor((ry[i]+m.lx/2.0)/dr))
    ##    mat[x][y] += 1
    ##print("Printing final matrix")
    ##print(mat.tolist())

    ##hist2d = np.histogram2d(rx,ry,size)
    ##print(hist2d[0].tolist())
    ##l = hist2d[0].tolist()

    #std = np.std(hist2d)
    #mean = np.mean(hist2d)
    #sm = 3*std+mean
    ##lsm = [[0 if x < sm else x for x in list] for list in l]
    ##print(lsm)

    ## Get atoms that are in spots:
    ## Algorithm: Create another matrix of size sizeXsize
    ## that will be used as a holder to identify all positions
    ## in the histogram that have intensities higher than
    ## sm = 3*stddev + mean. ie for every entry in mat that is
    ## greater than sm, increment the index in this new matrix
    ## as well as the indexes around that index. Then, go back
    ## through this matrix and set every index > 1 to True, and
    ## and everything else to False. These are the 'pixels' that
    ## we want to extract atoms from. Using GetAtomsInBin(...),
    ## put all the atoms in the True bins into a new model.
    ## Save that model!

    ## This says the spot size radius is 10 pixels.
    ## Need to look at image in Igor to figure this out.
    #tol = 10
    #tfmat = np.zeros((size,size))
    #for i in range(0,len(hist2d)):
    #    for j in range(0,len(hist2d[i])):
    #        kmin = i-tol
    #        kmax = i+tol
    #        if(kmin < 0): kmin = 0
    #        if(kmax > size-1): kmax = size-1
    #        for k in range(kmin,kmax):
    #            nmin = i-tol
    #            nmax = i+tol
    #            if(nmin < 0): nmin = 0
    #            if(nmax > size-1): nmax = size-1
    #            for n in range(nmin,nmax):
    #                if(hist2d[k][n] >= sm): tfmat[k][n] += 1
    #tfmat = [ [1 if x > 1 else 0 for x in ilist] for ilist in tfmat]

    #atoms = []
    #for i in range(0,size):
    #    for j in range(0,size):
    #        if(hist2d[i][j]):
    #            print(bindict[(i,j)])
    #            for k in range(0,len(bindict[(i,j)])):
    #                atoms.append(bindict[(i,j)][k][0])
    #                atoms.append(bindict[(i,j)][k][1])
    #atoms = list(set(atoms))
    #atoms = [m.atoms[i] for i in atoms]
    #new_model = Model("2D Hist extracted atoms\n", m.lx,m.ly,m.lz, atoms)

    #new_model.write_real_xyz('temp.xyz')
    #new_model.write_real_xyz()




def main():
    m = Model(sys.argv[1])
    rdf_2d(m,0.1)

if __name__ == '__main__':
    main()

