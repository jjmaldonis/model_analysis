from model import Model
import sys

def recenter_model(m):
    #xmax = 0.0
    #xmin = 0.0
    #ymax = 0.0
    #ymin = 0.0
    #zmax = 0.0
    #zmin = 0.0
    xmin = min(atom.coord[0] for atom in m.atoms)
    xmax = max(atom.coord[0] for atom in m.atoms)
    print("Original x-min/max = ({0},{1}".format(xmin,xmax))
    ymin = min(atom.coord[1] for atom in m.atoms)
    ymax = max(atom.coord[1] for atom in m.atoms)
    print("Original y-min/max = ({0},{1}".format(ymin,ymax))
    zmin = min(atom.coord[2] for atom in m.atoms)
    zmax = max(atom.coord[2] for atom in m.atoms)
    print("Original z-min/max = ({0},{1}".format(zmin,zmax))
    #for atom in m.atoms:
    #    if(atom.coord[0] > xmax): xmax = atom.coord[0]
    #    if(atom.coord[1] > ymax): ymax = atom.coord[1]
    #    if(atom.coord[2] > zmax): zmax = atom.coord[2]
    #    if(atom.coord[0] < xmin): xmin = atom.coord[0]
    #    if(atom.coord[1] < ymin): ymin = atom.coord[1]
    #    if(atom.coord[2] < zmin): zmin = atom.coord[2]
    xcenter = (xmax - xmin)/2.0
    ycenter = (ymax - ymin)/2.0
    zcenter = (zmax - zmin)/2.0
    print("Center is at ({0},{1},{2})".format(xcenter,ycenter,zcenter))

    for i in range(0,m.natoms):
        m.atoms[i].coord = ( m.atoms[i].coord[0] - xcenter, m.atoms[i].coord[1] - ycenter, m.atoms[i].coord[2] - zcenter )

def main():
    modelfile = sys.argv[1]
    m = Model(modelfile)

    recenter_model(m)
    m.write_our_xyz('temp.xyz')

if __name__ == "__main__":
    main()
