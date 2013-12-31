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
    
    def __repr__(self):
        #return str(self.id)+','+str(self.z)+',('+str(round(self.coord[0],3))+','+str(round(self.coord[1],3))+','+str(round(self.coord[2],3))+')'
        #return str(self.id)
        #if(type(self.z) != type('hi')):
        #    return str(self.id) + '\t' + str(self.z) + '\t' + str(self.coord[0]) + '\t' + str(self.coord[1]) + '\t' + str(self.coord[2])
        #elif(type(self.z) == type('hi')):
        #    return str(self.id) + '\t' + self.z + '\t' + str(self.coord[0]) + '\t' + str(self.coord[1]) + '\t' + str(self.coord[2])

        #VESTA Format: (no id)
        if(type(self.z) != type('hi')):
            return str(self.z) + '\t' + str(self.coord[0]) + '\t' + str(self.coord[1]) + '\t' + str(self.coord[2])
        elif(type(self.z) == type('hi')):
            return self.z + '\t' + str(self.coord[0]) + '\t' + str(self.coord[1]) + '\t' + str(self.coord[2])

    
    def convert_to_sym(self):
        self.z = znum2sym.z2sym(self.z)
        return self


def main():
    atom = Atom(0,13,0.0,0.0,0.0)
    print(atom)
    print(atom.convert_to_sym())

if __name__ == "__main__":
    main()
