import sys, os
from model import Model

def main():
    modelfiles = sys.argv[1:]
    head, tail = os.path.split(modelfiles[0])
    print(modelfiles)
    m0 = Model(modelfiles[0])
    for f in modelfiles[1:]:
        m = Model(f)
        for atom in m.atoms:
            if atom not in m0.atoms:
                m0.add(atom)
    print(m0.natoms)
    print(os.path.join(head, 'xtal.xyz'))
    m0.write(os.path.join(head, 'xtal.xyz'))

if __name__ == '__main__':
    main()
