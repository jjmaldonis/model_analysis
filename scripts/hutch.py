import math
import sys
import znum2sym
from atom import Atom

class Hutch(object):
    """ implements a hutch for a 3D atom model """

    def __init__(self, model=None):
        """ constructor """
        if model is None:
            raise Exception("No input model was given to Hutch().")
        if( model.lx != model.ly != model.lz):
            raise Exception("The model must be a cube!")
        super(Hutch, self).__init__()
        # Set basic self variables
        self.nhutchs = max(int(round(model.natoms**(1.0/3.0))),1)
        self.hutchsize = model.lx/self.nhutchs
        self.lx = model.lx
        self.ly = model.ly
        self.lz = model.lz
        # Create hutch dictionary
        self.hutchs = {}
        for i in range(0,self.nhutchs):
            for j in range(0,self.nhutchs):
                for k in range(0,self.nhutchs):
                    self.add_hutch((i,j,k))
        # Put atoms into their correct hutch
        for atom in model.atoms:
            self.add_atom(atom)
        self.check_hutches(model)

    def check_hutches(self,model):
        for atom in model.atoms:
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
        x = int(round( (atom.coord[0] + 0.5*self.lx) / self.hutchsize)) % self.nhutchs
        y = int(round( (atom.coord[1] + 0.5*self.ly) / self.hutchsize)) % self.nhutchs
        z = int(round( (atom.coord[2] + 0.5*self.lz) / self.hutchsize)) % self.nhutchs
        #if (x == 0): x = 1
        #if (y == 0): y = 1
        #if (z == 0): z = 1
        return (x,y,z)

    def add_atom(self,atom):
        hutch = self.hutch_position(atom)
        hutch = str(hutch[0]) + ',' + str(hutch[1]) + ',' + str(hutch[2])
        if hutch not in self.hutchs:
            #print(self.hutchs)
            raise Exception("You gave an improper hutch address! {0}".format(hutch))
        #if atom in self.hutchs[hutch]:
        #    raise Exception("Atom already exists in that hutch: {0}".format(atom))
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

def main():
    model = Model(sys.argv[1])
    hutch = Hutch(model)
    #print(hutch.hutchs)
    print(hutch.get_atoms_in_hutch((1,2,3)))
    print("Atoms near atoms[0]:")
    print(hutch.get_atoms_in_cutoff(hutch.atoms[0],1.5))

if __name__ == "__main__":
        main()
