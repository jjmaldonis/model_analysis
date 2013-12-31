
import os
import sys


def gen(modelfile,cutoff):
    with open(modelfile) as f:
        content = f.readlines()

    content = [line.strip().split() for line in content]
    content.pop(0) # remove comment header
    content.pop(-1) # remove last '-1' line from model file

    # convert to ints and floats
    for j in range(0,len(content)):
        for i in range(0,len(line)):
            try:
                content[j][i] = int(content[j][i])
            except:
                try:
                    content[j][i] = float(content[j][i])
                except:
                    pass
    
    box = content[0][0]
    content.pop(0) # remove world size line
    natoms = len(content)
    
    atom_types = []
    natom_types = {}
    for atom in content:
        if atom[0] not in atom_types:
            atom_types.append(atom[0])
            natom_types[atom[0]] = 1
        else:
            natom_types[atom[0]] += 1

    natom_types = [natom_types[i] for i in atom_types]
    for i in range(0,len(natom_types)):
        natom_types[i] = str(natom_types[i])
        atom_types[i] = str(atom_types[i])

    out = []
    out.append("# atom types")
    out.append(str(len(atom_types)) + ' ' + ' '.join(atom_types))
    out.append("# steps, # total atoms, # atoms of type1, # atoms of type2")
    out.append('1 '+str(natoms)+' '+" ".join(natom_types))
    out.append("# box size, # cut-off of neighbor")
    out.append(str(box)+' '+str(cutoff))

    return "\n".join(out)


def main():
    print(gen(sys.argv[1],3.5))

if __name__ == "__main__":
    main()

## atom types
#3 Zr Cu Al
## steps, # total atoms, # atoms of type1, # atoms of type2
#1 1425 713 642 70
## box size, # cut-off of neighbor
#   28.28450  3.5

