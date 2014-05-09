from model import Model
import sys

def frac2full(m):
    for i in range(0,m.natoms):
        m.atoms[i].coord = ( m.atoms[i].coord[0]*m.lx/2.0, m.atoms[i].coord[1]*m.ly/2.0, m.atoms[i].coord[2]*m.lz/2.0 )

def main():
    modelfile = sys.argv[1]
    m = Model(modelfile)

    frac2full(m)
    m.write_real_xyz()

if __name__ == "__main__":
    main()
