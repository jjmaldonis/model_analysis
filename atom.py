import znum2sym
import copy

class VoronoiPoly(object):
    def __init__(self):
        super(VoronoiPoly,self).__init__()
        self.index = None # Will be a list
        self.type = None # User defined string, e.g. 'Crystal-like'
        self.nnabsp = None # Will be a dictionary
        self.neighs = None # Will be a list of this atoms neighbors
        self.vol = None # This will be a float
        # Note: the center atom is this atom!

    def copy(self):
        new = VoronoiPoly()
        new.index = copy.copy(self.index)
        new.type = self.type
        new.nnabsp = copy.copy(self.nnabsp)
        new.neighs = copy.copy(self.neighs)
        new.vol = self.vol
        return new

    def __repr__(self):
        return '<{0}>'.format(','.join(str(x) for x in self.index))
    def __str__(self):
        return '<{0}>'.format(','.join(str(x) for x in self.index))

class Atom(object):
    def __init__(self, id, znum, x, y, z):
        super(Atom, self).__init__()
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

    def __eq__(self,other):
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

    #def __eq__(self,item):
    #    #if( self.id == item.id and self.z == item.z and self.coord == item.coord):
    #    if( self.z == item.z and self.coord == item.coord):
    #    #if( self.z == item.z and round(self.coord[0],10) == round(item.coord[0],10) and round(self.coord[1],10) == round(item.coord[1],10) and round(self.coord[2],10) == round(item.coord[2],10)):
    #        return True
    #    else:
    #        return False

    def __repr__(self):
        #return str(self.id)
        return 'atomID={0}'.format(self.id)
        #return str(self.id)+'\t'+str(self.z)+'\t('+str(round(self.coord[0],3))+','+str(round(self.coord[1],3))+','+str(round(self.coord[2],3))+')'
        #return str(self.id)
        #if(type(self.z) != type('hi')):
        #    return str(self.id) + '\t' + str(self.z) + '\t' + str(self.coord[0]) + '\t' + str(self.coord[1]) + '\t' + str(self.coord[2])
        #elif(type(self.z) == type('hi')):
        #    return str(self.id) + '\t' + self.z + '\t' + str(self.coord[0]) + '\t' + str(self.coord[1]) + '\t' + str(self.coord[2])

        #Format: (no znum): id x y z
        #return str(self.id) + '\t' + str(self.coord[0]) + '\t' + str(self.coord[1]) + '\t' + str(self.coord[2])

    @property
    def frac11(self,lx,ly,lz):
        """ Returns the coordinatese in fractional form between -1 and 1. Inputs are the world sizes """
        return str(self.coord[0]/lx*2)+'\t'+str(self.coord[1]/ly*2)+'\t'+str(self.coord[2]/lz*2)
    @property
    def frac01(self,lx,ly,lz):
        """ Returns the coordinatese in fractional form between 0 and 1. Inputs are the world sizes """
        return str(self.coord[0]/lx)+'\t'+str(self.coord[1]/ly)+'\t'+str(self.coord[2]/lz)
        
    def set_vp(self,index,center,type=None):
        self.vp.index = index
        self.vp.center = center
        self.vp.type = type

    def compute_vp_type(self,vp_dict):
        """ Computes, sets, and returns the vp type based on the input dictionary """
        try:
            self.vp
        except:
            raise Exception("The voronoi polyhedra for atom {0} has not been set!".format(self))
        for key in vp_dict:
            for vps in vp_dict[key]:
                found = True
                for i in range(0,4):
                    if(vps[i] != '*' and vps[i] != self.vp.index[i]):
                        found = False
                if(found):
                    self.vp.type = key[:-1]
                    return key[:-1]
        #raise Exception("The voronoi polyhedra for atom {0} could not be determined!".format(self))
        return "Undef"
    def realxyz(self):
        return str(znum2sym.z2sym(self.z))+'\t'+str(self.coord[0])+'\t'+str(self.coord[1])+'\t'+str(self.coord[2])
    
def main():
    atom = Atom(0,13,0.0,0.0,0.0)
    print(atom)
    print(atom.convert_to_sym())

if __name__ == "__main__":
    main()
