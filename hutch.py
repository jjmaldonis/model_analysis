from math import ceil,floor,sqrt
import sys
import znum2sym
from atom import Atom
from collections import defaultdict

def intify(func):
    def inner(*args,**kwargs):
        x = func(*args,**kwargs)
        x = int(x)
        return x
    return inner
ceil = intify(ceil)
floor = intify(floor)

class Hutch(object):
    """ Implements a hutch for a 3D atom model """

    def __init__(self, model=None, check=False):
        if model is None:
            raise Exception("No input model was given to Hutch().")
        if(model.xsize != model.ysize != model.zsize):
            raise Exception("The model must be a cube!")
        self.model = model # Keep a pointer to the parent model
        self.nhutchs = max(int(round(model.natoms**(1.0/3.0))),1)
        self.hutchsize = self.model.xsize/self.nhutchs
        # Create hutch dictionary
        self.hutchs = {}
        for i in range(0,self.nhutchs):
            for j in range(0,self.nhutchs):
                for k in range(0,self.nhutchs):
                    self.add_hutch((i,j,k))
        # Put atoms into their correct hutch
        for atom in model.atoms:
            self.add_atom(atom)
        if check:
            self.check_hutches(model)

    @property
    def xsize(self):
        return self.model.xsize
    @property
    def ysize(self):
        return self.model.ysize
    @property
    def zsize(self):
        return self.model.zsize

    def check_hutches(self,model):
        for atom in model.atoms:
            hutch = self._get_hutch(atom)
            if atom not in self.hutchs[hutch]:
                raise Exception("You have an error in your hutches!")

    def add_hutch(self,hutch):
        if hutch in self.hutchs:
            raise Exception("Trying to create a hutch that already exists: {0}".format(hutch))
        self.hutchs[hutch] =  []

    def _get_hutch(self,atom):
        """ Returns the hutch that the atom should be located in using a tuple """
        #x = int(round( (atom.coord[0] + 0.5*self.xsize) / self.hutchsize)) % self.nhutchs
        #y = int(round( (atom.coord[1] + 0.5*self.ysize) / self.hutchsize)) % self.nhutchs
        #z = int(round( (atom.coord[2] + 0.5*self.zsize) / self.hutchsize)) % self.nhutchs
        x = int(floor( (atom.coord[0] + 0.5*self.xsize) / self.hutchsize)) % self.nhutchs
        y = int(floor( (atom.coord[1] + 0.5*self.ysize) / self.hutchsize)) % self.nhutchs
        z = int(floor( (atom.coord[2] + 0.5*self.zsize) / self.hutchsize)) % self.nhutchs
        return (x,y,z)

    def add_atom(self,atom):
        hutch = self._get_hutch(atom)
        if hutch not in self.hutchs:
            raise Exception("You gave an improper hutch address! {0}".format(hutch))
        #if atom in self.hutchs[hutch]:
        #    raise Exception("Atom already exists in that hutch: {0}".format(atom))
        self.hutchs[hutch].append(atom)

    def remove_atom(self,atom):
        hutch = self._get_hutch(atom)
        try:
            self.hutchs[hutch].remove(atom)
        except IndexError:
            raise Exception("Trying to remove an atom that doesn't exist in the hutch found!")

    def get_atoms_in_same_hutch(self,atom):
        hutch = self._get_hutch(atom)
        return self.hutchs[hutch]

    def get_atoms_in_radius(self,theatom,cutoff):
        if(isinstance(cutoff,dict)):
            r2 = {key:r*r for key,r in cutoff.iteritems()}
            radius = max(r for _,r in cutoff.iteritems())
        else:
            radius = cutoff
            r2 = defaultdict(lambda: radius**2)
            cutoff = defaultdict(lambda: cutoff)
        x,y,z = theatom.coord
        hx = (x + self.xsize*0.5) /self.hutchsize
        hy = (y + self.ysize*0.5) /self.hutchsize
        hz = (z + self.zsize*0.5) /self.hutchsize
        hr = radius/self.hutchsize
        lhx = floor(hx-hr)
        lhy = floor(hy-hr)
        lhz = floor(hz-hr)
        uhx = ceil(hx+hr)
        uhy = ceil(hy+hr)
        uhz = ceil(hz+hr)
        hutches  = [(i,j,k,) for i in range(lhx,uhx) for j in range(lhy,uhy) for k in range(lhz,uhz)]
        hutches = [tuple([q if 0<=q<self.nhutchs else q%self.nhutchs for q in h]) for h in hutches] #pbc
        atoms = [ atom for hutch in hutches for atom in self.hutchs[hutch]
            if(self.model.dist2(atom,theatom) < r2[(atom.z,theatom.z)]) ]
        if(theatom in atoms): atoms.remove(theatom)
        return atoms

def main():
    model = Model(sys.argv[1])
    #hutch = Hutch(model)
    hutch = model.hutch
    #print(hutch.hutchs)
    print("Atoms near atoms[0]:")
    #print(hutch.get_atoms_in_radius(model.atoms[0],3.5))
    #print(model.get_atoms_in_cutoff(model.atoms[0],3.5))
    cutoff = {}
    cutoff[(46,46)] = 3.5
    cutoff[(14,14)] = 3.5
    cutoff[(46,14)] = 3.5
    cutoff[(14,46)] = 3.5
    atomid = 0
    for atomid in range(model.natoms):
    #for atomid in [2]:
        h = hutch.get_atoms_in_radius(model.atoms[atomid],cutoff)
        m = model.get_atoms_in_cutoff(model.atoms[atomid],cutoff)
        h.sort(key=lambda atom: atom.id)
        m.sort(key=lambda atom: atom.id)
        #print(h)
        #print(m)
        #print("In hutch routine")
        #for atom in h:
        #    if(atom not in m):
        #        print(model.dist(atom,model.atoms[atomid]))
        #        print(hutch._get_hutch(atom))
        #print("In model routine")
        #for atom in m:
        #    if(atom not in h):
        #        print(model.dist(atom,model.atoms[atomid]))
        #        print(hutch._get_hutch(atom))
        #        print(hutch.hutchs[hutch._get_hutch(atom)])
        assert h == m

if __name__ == "__main__":
    main()
