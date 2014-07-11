import sys
from model import Model

def main():
    mf1 = sys.argv[1]
    mf2 = sys.argv[2]

    m1 = Model(mf1)
    m2 = Model(mf2)

    for atom1 in m1.atoms:
        found = False
        for atom2 in m2.atoms:
            if(atom1 == atom2):
                found = True
                break
        else:
            raise Exception("Atom not found! {0}".format(atom1))
    else:
        print("All atoms were found!")

if __name__ == "__main__":
    main()

