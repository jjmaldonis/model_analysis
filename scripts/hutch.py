import math
import sys
import znum2sym
from atom import Atom

class Hutch(object):
    """ implements a hutch for a 3D atom model """

    def __init__(self, modelfile=None):
        """ constructor """
        if modelfile is None:
            raise Exception("No input model was given.")

        super(Hutch, self).__init__()
        self.hutchs = {}
        with open(modelfile) as f:
            content = f.readlines()
        content = [item.strip() for item in content]
        self.lx = eval(content[1].strip().split()[0])
        self.ly = eval(content[1].strip().split()[1])
        self.lz = eval(content[1].strip().split()[2])
        content.pop(0)
        content.pop(0)
        content.pop(-1)

        self.natoms = len(content)
        self.nhutchs = int(round(self.natoms**(1.0/3.0)))
        self.hutchsize = self.lx/self.nhutchs

        # Create hutch dictionary
        for i in range(0,self.nhutchs):
            for j in range(0,self.nhutchs):
                for k in range(0,self.nhutchs):
                    self.add_hutch((i,j,k))
        #print(self.hutchs)

        # content now contains only the atoms, as strings
        for i, item in enumerate(content):
            content[i] = content[i].split()
            id = i
            try:
                znum = eval(content[i][0])
            except:
                znum = znum2sym.sym2z(content[i][0])
            x = eval(content[i][1])
            y = eval(content[i][2])
            z = eval(content[i][3])
            content[i] = Atom(id, znum, x,y,z)
            self.add_atom(content[i])
        self.atoms = content

    def check_hutches(self):
        for atom in self.atoms:
            hutch = self.hutch_position(atom)
            hutch = str(hutch[0]) + ',' + str(hutch[1]) + ',' + str(hutch[2])
            if atom not in self.hutchs[hutch]:
                raise Exception("You have an error in your hutches!")
    
    def add_hutch(self,hutch):
        hutch = str(hutch[0]) + ',' + str(hutch[1]) + ',' + str(hutch[2])
        if hutch in self.hutchs:
            raise Exception("Trying to create a hutch that already exists: {0}".format(hutch))
        self.hutchs[hutch] =  []

    def hutch_position(self,atom):
        """ Returns the hutch that the atom should be located in using a tuple """
        x = int(round( (atom.coord[0] + 0.5*self.lx) / self.hutchsize) % self.nhutchs)
        y = int(round( (atom.coord[1] + 0.5*self.ly) / self.hutchsize) % self.nhutchs)
        z = int(round( (atom.coord[2] + 0.5*self.lz) / self.hutchsize) % self.nhutchs)
        if (x == 0): x = 1
        if (y == 0): y = 1
        if (z == 0): z = 1
        return (x,y,z)

    def add_atom(self,atom):
        hutch = self.hutch_position(atom)
        hutch = str(hutch[0]) + ',' + str(hutch[1]) + ',' + str(hutch[2])
        if hutch not in self.hutchs:
            #print(self.hutchs)
            raise Exception("You gave an improper hutch address! {0}".format(hutch))
        if atom in self.hutchs[hutch]:
            raise Exception("Atom already exists in that hutch.")
        self.hutchs[hutch].append(atom)

    def remove_atom(self,atom):
        hutch = self.hutch_position(atom)
        hutch = str(hutch[0]) + ',' + str(hutch[1]) + ',' + str(hutch[2])
        if atom not in self.hutchs[hutch]:
            raise Exception("Trying to remove an atom that doesn't exist in the hutch found!")
        self.hutchs[hutch].remove(atom)

    def get_atoms_in_hutch(self,hutch):
        hutch = str(hutch[0]) + ',' + str(hutch[1]) + ',' + str(hutch[2])
        return self.hutchs[hutch]

    def get_atoms_in_cutoff(self,atom,cutoff):
        x_start = atom.coord[0] - cutoff
        y_start = atom.coord[1] - cutoff
        z_start = atom.coord[2] - cutoff
        x_end   = atom.coord[0] + cutoff
        y_end   = atom.coord[1] + cutoff
        z_end   = atom.coord[2] + cutoff
        if(x_start < -self.lx/2.0): x_start = x_start + self.lx #PBC
        if(y_start < -self.ly/2.0): y_start = y_start + self.ly #PBC
        if(z_start < -self.lz/2.0): z_start = z_start + self.lz #PBC
        if(x_end > self.lx/2.0): x_end = x_end - self.lx #PBC
        if(y_end > self.ly/2.0): y_end = y_end - self.ly #PBC
        if(z_end > self.lz/2.0): z_end = z_end - self.lz #PBC
        atom_s = Atom(0,0,x_start,y_start,z_start)
        hutch_s = self.hutch_position(atom_s)
        i_start = hutch_s[0]
        j_start = hutch_s[1]
        k_start = hutch_s[2]
        atom_e = Atom(0,0,x_end,y_end,z_end)
        hutch_e = self.hutch_position(atom_e)
        i_end = hutch_e[0]
        j_end = hutch_e[1]
        k_end = hutch_e[2]
        list = []
        for i in range(0,self.nhutchs):
            if(i_start <= i_end):
                if(i < i_start or i > i_end): continue
            else:
                if(i < i_start and i > i_end): continue
            for j in range(0,self.nhutchs):
                if(j_start <= j_end):
                    if(j < j_start or j > j_end): continue
                else:
                    if(j < j_start and j > j_end): continue
                for k in range(0,self.nhutchs):
                    if(k_start <= k_end):
                        if(k < k_start or k > k_end): continue
                    else:
                        if(k < k_start and k > k_end): continue

                    list = list + self.get_atoms_in_hutch((i,j,k))
        list.remove(atom)
        list = [atomi for atomi in list if self.dist(atom,atomi) <= cutoff]
        #for atomi in list:
        #    if(self.dist(atom,atomi) > cutoff):
        #        list.remove(atomi)
        #    else:
        #        print("dist={0}".format(self.dist(atom,atomi)))
        return list

    def dist(self,atom1,atom2):
        x2 = (atom1.coord[0] - atom2.coord[0])**2
        y2 = (atom1.coord[1] - atom2.coord[1])**2
        z2 = (atom1.coord[2] - atom2.coord[2])**2
        return math.sqrt(x2+y2+z2)
    def get_all_dists(self):
        dists = []
        for atomi in self.atoms:
            for atomj in self.atoms[self.atoms.index(atomi)+1:]:
                dists.append([atomi,atomj,self.dist(atomi,atomj)])
        return dists

        

    def get_all_atoms(self):
        return self.atoms
    def get_natoms(self):
        return self.natoms

def main():
    hutch = Hutch(sys.argv[1])
    hutch.check_hutches()
    #print(hutch.hutchs)
    print(hutch.get_atoms_in_hutch((1,2,3)))
    print("Atoms near atoms[0]:")
    print(hutch.get_atoms_in_cutoff(hutch.atoms[0],1.5))

if __name__ == "__main__":
        main()
