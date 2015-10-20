import sys
from vor import Vor
from model import Model


def icofrac(m, modelfile, cutoff):
    nbins = 6
    del_bin = 100.0/nbins
    fracs = []
    for atom in m.atoms:
        fracs.append((float(atom.vp.index[2])/float(sum(atom.vp.index))*100.0))
        bin = int( (float(atom.vp.index[2])/float(sum(atom.vp.index))*100.0) /(100.0/(nbins-1)))
        atom.z = bin+1
    fracs.sort()
    print('Min %: {0}. Max %: {1}'.format(min(fracs),max(fracs)))

def main():
    modelfile = sys.argv[1]
    cutoff = 3.5
    m = Model(modelfile)
    vor_instance = Vor()
    vor_instance.runall(modelfile,cutoff)
    vor_instance.set_atom_vp_indexes(m)
    icofrac(m,modelfile,cutoff)
    m.write_real_xyz()

if __name__ == "__main__":
    main()
