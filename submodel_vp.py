import sys
from model import Model
from voronoi_3d import voronoi_3d

def icofrac(m):
    """ This modifies the model. """
    nbins = 6
    del_bin = 100.0/nbins
    fracs = []
    for atom in m.atoms:
        fracs.append((float(atom.vp.index[2])/float(sum(atom.vp.index))*100.0))
        bin = int( (float(atom.vp.index[2])/float(sum(atom.vp.index))*100.0) /(100.0/(nbins-1)))
        atom.z = bin+1
    fracs.sort()
    #print('Min %: {0}. Max %: {1}'.format(min(fracs),max(fracs)))


def main():
    # NOTE: Cutoff can either be a single integer or it
    # can be a dictionary where the keys are two-tuples
    # of atomic numbers (e.g. (40,13)=3.5 for Zr,Al).
    modelfile = sys.argv[1]
    submodelfile = sys.argv[2]
    rotatedsubmodelfile = sys.argv[3]
    m = Model(modelfile)
    try:
        cut = 3.5 #float(sys.argv[2])
        cutoff = {}
        for z1 in m.atomtypes:
            for z2 in m.atomtypes:
                cutoff[(z1,z2)] = cut
                cutoff[(z2,z1)] = cut
    except:
        print("You didn't input a cutoff so you much define it in the code.")
    voronoi_3d(m,cutoff)

    subm = Model(submodelfile)
    rotsubm = Model(rotatedsubmodelfile)
    for atom in subm.atoms:
        if atom in m.atoms:
            rotsubm.atoms[atom.id].vp = m.atoms[m.atoms.index(atom)].vp
        else:
            print("Couldn't find atom {0} in full model!".format(atom))
    icofrac(rotsubm)
    #rotsubm.write_real_xyz()
    rotsubm.write_real_xyz(rotatedsubmodelfile[:-4]+'.icofrac.real.xyz')

if __name__ == '__main__':
    main()
