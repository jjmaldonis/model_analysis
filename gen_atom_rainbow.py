

from atom import Atom


def main():
    atoms = []
    for i in range(1,101):
        atoms.append(Atom(i,i,0,0,i))
        atoms[i-1] = atoms[i-1].convert_to_sym()
        print atoms[i-1]


if __name__ == '__main__':
    main()
