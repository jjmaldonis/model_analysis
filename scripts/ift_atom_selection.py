import sys
#import numpy as np
from basic_model import Model
import math

def atom_selection(modelfile, intensityfile, outbase=None):
    npix = 256
    #intensities = np.zeros((npix,npix,npix),dtype=float)
    intensities = [[[0.0 for i in range(npix)] for i in range(npix)] for i in range(npix)]
    m = Model(modelfile)
    maxx = 0.0
    with open(intensityfile) as f:
        for l,line in enumerate(f):
            line = line.strip().split()
            for i,x in enumerate(line):
                x = float(x)
                if(x > maxx): maxx = x
                intensities[i%npix][int(i/npix)][l] = x
            print("Read in line {0}.".format(l))

    #ints = np.zeros((m.natoms),dtype=float)
    ints = [0.0 for i in range(m.natoms)]
    for i,atom in enumerate(m.atoms):
        xbin = int((m.lx/2+atom.coord[0])/m.lx*npix)
        ybin = int((m.ly/2+atom.coord[1])/m.ly*npix)
        zbin = int((m.ly/2+atom.coord[2])/m.lz*npix)
        atom.intensity = intensities[xbin][ybin][zbin]/maxx
        ints[i] = intensities[xbin][ybin][zbin]/maxx

    #print("np.mean: {0}".format(np.mean(ints)))
    #print("np.stdev: {0}".format(np.std(ints)))
    #print("Variance: {0}".format(np.var(ints)))
    #print("Median: {0}".format(np.median(ints)))
    mean = sum(ints)/len(ints)
    print("Mean: {0}".format(mean))
    stdev = 0.0
    for x in ints:
        stdev += (x-mean)**2
    stdev = math.sqrt(stdev/len(ints))
    print("Stdev: {0}".format(stdev))
    mini = mean + stdev/2.0
    print("Accepting atoms with intensity above {0}".format(mini))
    atoms = [atom for atom in m.atoms if atom.intensity > mini]
    print("Found {0} atoms".format(len(atoms)))
    if(outbase == None):
        outbase = 'temp'
    Model(m.comment,m.lx,m.ly,m.lz,atoms).write_cif(outbase+'.cif')
    Model(m.comment,m.lx,m.ly,m.lz,atoms).write_our_xyz(outbase+'.xyz')

def main():
    modelfile = sys.argv[1]
    intensityfile = sys.argv[2]
    atom_selection(modelfile, intensityfile)


if __name__ == '__main__':
    main()
