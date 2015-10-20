import sys
#import numpy as np
from basic_model import Model
import math

def atom_selection(modelfile, intensityfile, npix, outbase=None):
    intensities = [[[0.0 for i in range(npix)] for i in range(npix)] for i in range(npix)]
    m = Model(modelfile)
    maxx = 0.0
    mean = 0.0
    with open(intensityfile) as f:
        for l,line in enumerate(f):
            line = line.strip().split()
            if(l == 0):
                npixx, npixy, npixz = tuple([int(x) for x in line])
                if(npixx == npixy == npixz):
                    npix = npixx
                else:
                    pass
            else:
                for i,x in enumerate(line):
                    x = float(x)
                    if(x > maxx): maxx = x
                    intensities[i%npix][int(i/npix)][l-1] = x
                    mean += x
                print("Read in line {0} (python).".format(l-1))
    if(maxx == 0.0):
        print("Something went wrong somewhere. Check the stdev and mgrid gfx files first.")
        return 1
    mean /= npix**3
    mean /= maxx
    print("Mean of entire box is {0}".format(mean))
    stdev = 0.0
    with open(intensityfile) as f:
        for l,line in enumerate(f):
            line = line.strip().split()
            if(l == 0):
                npixx, npixy, npixz = tuple([int(x) for x in line])
                if(npixx == npixy == npixz):
                    npix = npixx
                else:
                    pass
            else:
                for i,x in enumerate(line):
                    x = float(x)
                    stdev += (x/maxx-mean)**2
                print("Read in line {0} again (python).".format(l-1))
    stdev = math.sqrt(stdev/npix**3)

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
    #mean = sum(ints)/len(ints)
    print("Mean: {0}".format(mean))
    #stdev = 0.0
    #for x in ints:
    #    stdev += (x-mean)**2
    #stdev = math.sqrt(stdev/len(ints))
    print("Stdev: {0}".format(stdev))
    #mini = mean + stdev/2.0
    mini = mean + stdev/1.0
    print("Accepting atoms with intensity above {0}".format(mini))
    atoms = [atom for atom in m.atoms if atom.intensity > mini]
    print("Found {0} atoms".format(len(atoms)))
    if(outbase == None):
        outbase = 'temp'
    Model(m.comment,m.lx,m.ly,m.lz,atoms).write_cif(outbase+'.cif')
    Model(m.comment,m.lx,m.ly,m.lz,atoms).write_our_xyz(outbase+'.xyz')
    Model(m.comment,m.lx,m.ly,m.lz,atoms).write_real_xyz(outbase+'.real.xyz')
    return 0

def main():
    modelfile = sys.argv[1]
    intensityfile = sys.argv[2]
    atom_selection(modelfile, intensityfile, 256)


if __name__ == '__main__':
    main()
