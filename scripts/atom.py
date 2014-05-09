import znum2sym

class VoronoiPoly(object):
    def __init__(self):
        super(VoronoiPoly,self).__init__()
        self.index = None # Will be a list
        self.type = None # User defined string, e.g. 'Crystal-like'
        self.nnabsp = None # Will be a dictionary
        self.neighs = None # Will be a list of this atoms neighbors
        self.vol = None # This will be a float
        # Note: the center atom is this atom!

class Atom(object):
    """ atom class """

    def __init__(self, id, znum, x, y, z):
        """ constructor 
        id = atom id (from 0 to natoms)
        znum = atomic number
        x = x coordinate
        y = y coordinate
        z = z coordinate """

        super(Atom, self).__init__()
        self.id = id
        if(type(znum) == type(0)):
            self.z = znum
        elif(type(znum) == type('string')):
            self.z = znum2sym.sym2z(znum)
        self.coord = (x,y,z)
        self.vp = VoronoiPoly()
        self.neighs = None # Neighbors. List when set
        self.cn = None # Coordination number. int when set
        self.sym = znum2sym.z2sym(self.z) #atomic symbol
    

    def __eq__(self,item):
        #if( self.id == item.id and self.z == item.z and self.coord == item.coord):
        if( self.z == item.z and self.coord == item.coord):
            return True
        else:
            return False
    def __ne__(self,item):
        return not self == item

    def __repr__(self):
        return str(self.id)
        #return str(self.id)+'\t'+str(self.z)+'\t('+str(round(self.coord[0],3))+','+str(round(self.coord[1],3))+','+str(round(self.coord[2],3))+')'
        #return str(self.id)
        #if(type(self.z) != type('hi')):
        #    return str(self.id) + '\t' + str(self.z) + '\t' + str(self.coord[0]) + '\t' + str(self.coord[1]) + '\t' + str(self.coord[2])
        #elif(type(self.z) == type('hi')):
        #    return str(self.id) + '\t' + self.z + '\t' + str(self.coord[0]) + '\t' + str(self.coord[1]) + '\t' + str(self.coord[2])

        #Format: (no znum): id x y z
        #return str(self.id) + '\t' + str(self.coord[0]) + '\t' + str(self.coord[1]) + '\t' + str(self.coord[2])
    def vesta(self):
        """ returns a string with the atoms information in the format for the vesta input file """
        if(type(self.z) != type('hi')):
            return znum2sym.z2sym(self.z) + '\t' + str(self.coord[0]) + '\t' + str(self.coord[1]) + '\t' + str(self.coord[2])
        elif(type(self.z) == type('hi')):
            return self.z + '\t' + str(self.coord[0]) + '\t' + str(self.coord[1]) + '\t' + str(self.coord[2])
    def ourxyz(self):
        return str(self.z)+'\t'+str(self.coord[0])+'\t'+str(self.coord[1])+'\t'+str(self.coord[2])
    def realxyz(self):
        return str(znum2sym.z2sym(self.z))+'\t'+str(self.coord[0])+'\t'+str(self.coord[1])+'\t'+str(self.coord[2])
    def frac11(self,lx,ly,lz):
        """ returns the coordinatese in fractional form between -1 and 1
            inputs are the world sizes """
        return str(self.coord[0]/lx*2)+'\t'+str(self.coord[1]/ly*2)+'\t'+str(self.coord[2]/lz*2)
    def frac01(self,lx,ly,lz):
        """ returns the coordinatese in fractional form between 0 and 1
            inputs are the world sizes """
        return str(self.coord[0]/lx)+'\t'+str(self.coord[1]/ly)+'\t'+str(self.coord[2]/lz)
        
    def set_vp(self,index,center,type=None):
        self.vp.index = index
        self.vp.center = center
        self.vp.type = type
    def set_vp_type(self,type):
        self.vp.type = type

    def compute_vp_type(self,vp_dict):
        """ Computes, sets, and returns the vp type based on the input dictionary, which is constructed from a file """
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
    
    def convert_to_sym(self):
        self.z = znum2sym.z2sym(self.z)
        return self

    def set_id(self,id):
        self.id = id
        return self
    def set_coord(self,x,y,z):
        self.coord = (x,y,z)
        return self

    def set_znum(self,z):
        self.z = z
        return self
    def get_znum(self):
        return self.z

def main():
    atom = Atom(0,13,0.0,0.0,0.0)
    print(atom)
    print(atom.convert_to_sym())

if __name__ == "__main__":
    main()
