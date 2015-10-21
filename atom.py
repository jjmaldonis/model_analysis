import znum2sym
import copy

class VoronoiPoly(object):
    """ A structure for holding the Voronoi Polyhedron surrounding an atom. """
    def __init__(self, index=None, type=None, nnabsp=None, neighs=None, volume=None):
        self.index = index # List
        self.type = type # User defined string, e.g. 'Crystal-like'
        self.nnabsp = nnabsp # Dictionary
        self.neighs = neighs # List of this atoms neighbors
        self.volume = volume # Float

    def copy(self):
        new = VoronoiPoly()
        new.index = copy.copy(self.index)
        new.type = self.type
        new.nnabsp = copy.copy(self.nnabsp)
        new.neighs = copy.copy(self.neighs)
        new.volume = self.volume
        return new

    def compute_type(self, vp_dict):
        """ Computes, sets, and returns the vp type based on the input dictionary.
            self.index must already be set. """
        for key in vp_dict:
            for vps in vp_dict[key]:
                found = True
                for i in range(4):
                    if(vps[i] != '*' and vps[i] != self.index[i]):
                        found = False
                if(found):
                    self.type = key[:-1]
                    return key[:-1]
        return "Undef"

    def __repr__(self):
        return '<{0}>'.format(' '.join(str(x) for x in self.index))
    def __str__(self):
        return '<{0}>'.format(' '.join(str(x) for x in self.index))

class Atom(object):
    def __init__(self, id, znum, x, y, z):
        self.id = id
        if(isinstance(znum,int)):
            self.z = znum
        elif(isinstance(znum,str)):
            self.z = znum2sym.sym2z(znum)
        self.coord = (x,y,z)
        self.vp = VoronoiPoly()
        self.neighs = None # Neighbors. List when set
        self.cn = None # Coordination number. int when set
        self.sym = znum2sym.z2sym(self.z) #atomic symbol

    def __eq__(self, other):
        x,y,z = [round(xx,6) for xx in self.coord]
        a,b,c = [round(xx,6) for xx in other.coord]
        return self.z == other.z and x==a and y==b and z==c
        #return (self.z == other.z and self.coord == other.coord)

    def copy(self):
        new = Atom(self.id, self.z, self.coord[0], self.coord[1], self.coord[2])
        new.vp = self.vp.copy()
        new.neighs = copy.copy(self.neighs)
        new.cn = self.cn
        return new

    def __repr__(self):
        return 'atomID={0}'.format(self.id)
        return "{0}\t{1}\t({2}, {3}, {4})".format(self.id, self.z, round(self.coord[0],3), round(self.coord[1],3), round(self.coord[2],3))

    def realxyz(self):
        """ Prints the atom in xyz file format. """
        return "{0}\t{1}\t{2}\t{3}".format(self.sym, self.coord[0], self.coord[1], self.coord[2])

    @property
    def frac11(self, lx, ly, lz):
        """ Returns the coordinatese in fractional form between -1 and 1. Inputs are the world sizes. """
        return (self.coord[0]/lx*2., self.coord[1]/ly*2., self.coord[2]/lz*2.)

    @property
    def frac01(self, lx, ly, lz):
        """ Returns the coordinatese in fractional form between 0 and 1. Inputs are the world sizes. """
        return (self.coord[0]/lx+0.5, self.coord[1]/ly+0.5, self.coord[2]/lz+0.5)
        
    
def main():
    atom = Atom(0,13,0.0,0.0,0.0)
    print(atom)

if __name__ == "__main__":
    main()
