import sys
import znum2sym
from  model import Model


def main():
    m = Model(sys.argv[1])

    for i,atom in enumerate(m.atoms):
        m.atoms[i] = atom.convert_to_sym()

    content = []
    content.append(str(m.natoms))
    content.append(str(m.lx)+' '+str(m.ly)+' '+str(m.lz)+' '+m.comment)
    for i,atom in enumerate(m.atoms):
        content.append( atom.z + ' ' + str(atom.coord[0]) + ' ' + str(atom.coord[1]) + ' ' + str(atom.coord[2]) )

    for line in content:
        print(line.strip())



if __name__ == "__main__":
    main()
