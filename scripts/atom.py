import znum2sym

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
        self.z = znum
        self.coord = (x,y,z)
    

    def __eq__(self,item):
        #if( self.id == item.id and self.z == item.z and self.coord == item.coord):
        if( self.z == item.z and self.coord == item.coord):
            return True
        else:
            return False
    def __ne__(self,item):
        return not self == item

    def __repr__(self):
        return str(self.id)+'\t'+str(self.z)+'\t('+str(round(self.coord[0],3))+','+str(round(self.coord[1],3))+','+str(round(self.coord[2],3))+')'
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
    def frac11(self,lx,ly,lz):
        """ returns the coordinatese in fractional form between -1 and 1
            inputs are the world sizes """
        return str(self.coord[0]/lx*2)+'\t'+str(self.coord[1]/ly*2)+'\t'+str(self.coord[2]/lz*2)
    def frac01(self,lx,ly,lz):
        """ returns the coordinatese in fractional form between 0 and 1
            inputs are the world sizes """
        return str(self.coord[0]/lx)+'\t'+str(self.coord[1]/ly)+'\t'+str(self.coord[2]/lz)
        
    def set_vp(self,vp):
        self.vp = vp

    def get_vp_type(self,vp_dict):
        try:
            self.vp
        except:
            raise Exception("The voronoi polyhedra for atom {0} has not been set!".format(self))
        for key in vp_dict:
            for vps in vp_dict[key]:
                found = True
                for i in range(0,4):
                    if(vps[i] != '*' and vps[i] != self.vp[i]):
                        found = False
                if(found):
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
