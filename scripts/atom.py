

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
        return str(self.id)
