import sys
from vor import Vor
from model import Model


def generate_map(cutoff,modelfile):
    """ cutoff is for vor.f90 """
    vor_instance = Vor()
    vor_instance.runall(modelfile,cutoff)
    index = vor_instance.get_index()

    icofrac = []
    index = [index[i].strip().split() for i in range(0,len(index))]
    for i in range(0,len(index)):
        for j in range(0,len(index[i])):
            try:
                index[i][j] = int(index[i][j])
            except:
                try:
                    index[i][j] = float(index[i][j])
                except:
                    pass
        index[i] = index[i][6:11]
        icofrac.append(int( float(index[i][2])/float(sum(index[i]))*100 ))

    model = Model(modelfile)
    atoms = model.get_atoms()
    
    atoms = [atom.set_znum(icofrac[i]) for i,atom in enumerate(atoms)]
    #for atom in atoms:
    #    if atom.get_znum() == 0:
    #        atom.set_znum(1)
    nbins = 6
    del_bin = 100.0/nbins
    for atom in atoms:
        for i in range(1,nbins+1):
            if( atom.z <= i*del_bin):
                atom.z = i
                break

    atoms = [atom.convert_to_sym() for i,atom in enumerate(atoms)]

    print(model.natoms)
    print('{0} {1} {2}'.format(model.lx,model.ly,model.lz))
    for atom in atoms:
        print atom.vesta()

def icofrac(m, modelfile, cutoff):
    vor_instance = Vor()
    vor_instance.runall(modelfile,cutoff)
    vor_instance.set_atom_vp_indexes(m)
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
    modelfile = sys.argv[1]
    m = Model(modelfile)
    icofrac(m,modelfile,3.5)
    m.write_real_xyz()
    #generate_map(3.5,sys.argv[1])

if __name__ == "__main__":
    main()
