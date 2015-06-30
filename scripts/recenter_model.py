from model import Model
import sys

def recenter_model(m,debug=False):
    xmin = min(atom.coord[0] for atom in m.atoms)
    xmax = max(atom.coord[0] for atom in m.atoms)
    ymin = min(atom.coord[1] for atom in m.atoms)
    ymax = max(atom.coord[1] for atom in m.atoms)
    zmin = min(atom.coord[2] for atom in m.atoms)
    zmax = max(atom.coord[2] for atom in m.atoms)
    if(debug):
        print("Original x-min/max = ({0},{1}".format(xmin,xmax))
        print("Original y-min/max = ({0},{1}".format(ymin,ymax))
        print("Original z-min/max = ({0},{1}".format(zmin,zmax))
    xcenter = xmin + (xmax - xmin)/2.0
    ycenter = ymin + (ymax - ymin)/2.0
    zcenter = zmin + (zmax - zmin)/2.0
    if(debug):
        print("Center was at ({0},{1},{2})".format(xcenter,ycenter,zcenter))

    for i in range(0,m.natoms):
        m.atoms[i].coord = ( m.atoms[i].coord[0] - xcenter, m.atoms[i].coord[1] - ycenter, m.atoms[i].coord[2] - zcenter )
    if(debug):
        xmin = min(atom.coord[0] for atom in m.atoms)
        xmax = max(atom.coord[0] for atom in m.atoms)
        print("Original x-min/max = ({0},{1}".format(xmin,xmax))
        ymin = min(atom.coord[1] for atom in m.atoms)
        ymax = max(atom.coord[1] for atom in m.atoms)
        print("Original y-min/max = ({0},{1}".format(ymin,ymax))
        zmin = min(atom.coord[2] for atom in m.atoms)
        zmax = max(atom.coord[2] for atom in m.atoms)
        print("Original z-min/max = ({0},{1}".format(zmin,zmax))
        xcenter = round(xmin + (xmax - xmin)/2.0,15)
        ycenter = round(ymin + (ymax - ymin)/2.0,15)
        zcenter = round(zmin + (zmax - zmin)/2.0,15)
        print("Center is at ({0},{1},{2})".format(xcenter,ycenter,zcenter))

def main():
    modelfile = sys.argv[1]
    m = Model(modelfile)

    recenter_model(m)
    m.write_our_xyz('temp.xyz')

if __name__ == "__main__":
    main()
