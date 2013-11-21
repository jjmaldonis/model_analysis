import sys

def main():
    with open( sys.argv[2]) as f:
        lines = [line.split() for line in f]
    for line in lines:
        for i in range(0,len(line)):
            try:
                line[i] = int(line[i])
            except:
                try:
                    line[i] = float(line[i])
                except:
                    pass
    # Now all the lines are read in and in integer or float form
    # We care about columns 7-10, which are i3, i4, i5, i6.

    # Now open the input parameter file.
    # Should be of the form:
    # Crystal:
    #     0,2,8,*
    with open( sys.argv[1]) as f:
        vp_params = [line.split(',') for line in f]
    # For each voronoi polyhedra 'structure', change it to an int if it's not a *
    vp_dict = {}
    for vps in vp_params:
        for i in range(0,len(vps)):
            vps[i] = vps[i].strip()
        if(len(vps) == 1):
            current_entry = vps[0]
            vp_dict[current_entry] = []
        elif(len(vps) == 4):
            for i in range(0,4):
                if(vps[i] != '*'):
                    vps[i] = int(vps[i])
            vp_dict[current_entry].append(vps)
        else:
            sys.exit("ERROR IN INPUT VP PARAMETERS. STOPPING.")

    # Create a new dictionary to start all the atoms that are crystalline, icosahedra, etc.
    atom_dict = {}
    for key in vp_dict:
        atom_dict[key] = []

    # Sort atoms
    for line in lines:
        for key in vp_dict:
            for vps in vp_dict[key]:
                found = True
                for i in range(0,4):
                    if(vps[i] != '*' and vps[i] != line[i+6]):
                        found = False
                if(found):
                    atom_dict[key].append(line)

    print("Inputs:")
    for key in vp_dict:
        print(key)
        for vps in vp_dict[key]:
            print vps

    print("\nOutputs:")

    sum = 0
    for key in atom_dict:
        sum = sum + len(atom_dict[key])
    print("Total number of atoms categorized: " + str(sum))
    for key in atom_dict:
        print(key + ' ' + str(len(atom_dict[key])))
        for line in atom_dict[key]:
            for i in range(0,len(line)):
                if(type(line[i]) != type('hi')):
                    line[i] = str(line[i])
            line = "\t".join(line)
            print line


if __name__ == "__main__":
    main()
