import sys
from vor import Vor
from model import Model


class GradientGraph(object):
    """ """

    def __init__(self):
        super(GradientGraph,self).__init__()

    def generate_map(self,cutoff,modelfile):
        vor_instance = Vor()
        vor_instance.runall(cutoff,modelfile)
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
        nbins = 5
        del_bin = 100.0/nbins
        for atom in atoms:
            for i in range(1,nbins+1):
                if( atom.z <= i*del_bin):
                    atom.z = i
                    break

        atoms = [atom.convert_to_sym() for i,atom in enumerate(atoms)]

        for atom in atoms:
            print atom

def main():
    gg = GradientGraph()
    gg.generate_map(3.5,sys.argv[1])

if __name__ == "__main__":
    main()
