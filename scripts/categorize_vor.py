
# Categorizes the voronoi polyhedra according to the input parameter file, 
# which must be the first arg passed to this script.
# The second arg must be a _index.out file.
# Output is printed to screen.

import sys
import vor
from model import Model
from voronoi_3d import voronoi_3d

def load_index_file(indexfile):
    # Open _index.out file.
    with open(indexfile) as f:
        index = f.readlines()
    return index

def load_param_file(paramfile):
    """ Returns a dictionary in the form:
        {'Crystal-like:': [[0, 4, 4, '*'], [0, 5, 2, '*']], 'Icosahedra-like:': [[0, 2, 8, '*'], [0, 1, 10, '*'], [0, 0, 12, 0], [0, 0, 12, 2]], 'Mixed:': [[0, 3, 6, '*'], [0, 3, 7, 2], [1, 2, 5, 4]]} """
    # Now open the input parameter file.
    # Should be of the form:
    # Crystal:
    #     0,2,8,*
    with open(paramfile) as f:
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
    return vp_dict

def generate_atom_dict(indexes,vp_dict):
    """ Generate a new dictionary to store all the atoms that are crystalline, icosahedra, etc.
        Returns a dictionary in the form:
        { 'Mixed:': [list of mixed atoms], 'Crystal-like:', [list of crystal-like atoms], etc}
        where each atom has the form of the _index.out lines:
        [id, .., .., .., .., .., i3, i4, i5, i6, .....] """
    indexes = [line.strip().split() for line in indexes]
    for line in indexes:
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

    # Generate a new dictionary to store all the atoms that are crystalline, icosahedra, etc.
    atom_dict = {}
    for key in vp_dict:
        atom_dict[key] = []
    atom_dict["Undef:"] = []

    # Sort atoms, putting them in their correct dictionary place
    for line in indexes:
        if('id,znum,nneighs,nneighst1,nneighst2,nneighs3,n3,n4,n5,n6,n7,n8,n9,n10,vol' not in line):
            for key in vp_dict:
                for vps in vp_dict[key]:
                    found = True
                    for i in range(0,4):
                        if(vps[i] != '*' and vps[i] != line[i+6]):
                            found = False
                    if(found):
                        atom_dict[key].append(line)
            found = False
            for key in vp_dict:
                if line in atom_dict[key]:
                    found = True
            if not found:
                atom_dict["Undef:"].append(line)
    return atom_dict


def set_atom_vp_types(model,vp_dict):
    """ saves the voronoi polyhedra type for each atom to the atom in the model """
    for atom in model.atoms:
        for key in vp_dict:
            for vps in vp_dict[key]:
                found = True
                for i in range(0,4):
                    if(vps[i] != '*' and vps[i] != atom.vp.index[i]):
                        found = False
                if(found):
                    atom.vp.type = key[:-1]
        if atom.vp.type == None:
            atom.vp.type = 'Undef'
    #for i,atomi in enumerate(model.atoms):
    #    found = False
    #    for key in atom_dict:
    #        for line in atom_dict[key]:
    #            if line[0] == i:
    #                vp = tuple(line[6:14]) # This saves n3 - n10
    #                center = line[1]
    #                found = True
    #    if not found:
    #        raise Exception("Error! Couldn't find an atom when trying to save it's VP.")
    #    model.atoms[i].set_vp(vp)


def printVorCats(atom_dict,vp_dict):
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
    atom_dict = atom_dict
    for key in atom_dict:
        print(key + ' ' + str(len(atom_dict[key])) + ' ' + str(round(len(atom_dict[key])/float(sum)*100,1))+'%')
        for line in atom_dict[key]:
            #print(line)
            for i in range(0,len(line)):
                if(type(line[i]) != type('hi')):
                    line[i] = str(line[i])
            line = "\t".join(line)
            print(line)
        for line in atom_dict[key]:
            for i in range(0,len(line)):
                try:
                    line[i] = int(line[i])
                    #print("Converted "+str(line[i]))
                except:
                    try:
                        line[i] = float(line[i])
                        #print("Converted "+str(line[i]))
                    except:
                        pass
                        #print("Could not convert: "+str(line[i]))


def vor_stats(m):
    #counter = {}
    #species_counter = {}
    #for atom in m.atoms:
    #    counter[atom.vp.type] = counter.get(atom.vp.type, 0) + 1
    #    if atom.sym not in species_counter:
    #        species_counter[atom.sym] = {}
    #    species_counter[atom.sym][atom.vp.type] = species_counter[atom.sym].get(atom.vp.type, 0) + 1
    #print(counter)
    #print(species_counter)
    cats = {}
    for atom in m.atoms:
        if atom.vp.type not in cats:
            cats[atom.vp.type] = {}
        cats[atom.vp.type][atom.sym] = cats[atom.vp.type].get(atom.sym,0) + 1
        cats[atom.vp.type]["Total"] = cats[atom.vp.type].get("Total",0) + 1
    for key in cats:
        print("{0}: {1}".format(key,cats[key]["Total"]))
        for elem in cats[key]:
            if(elem != "Total"): print("   {0}: {1}".format(elem,cats[key][elem]))

def print_all(m):
    for atom in m.atoms:
        print("{0} {1} {2}".format(atom,atom.vp.index,atom.vp.type))

        


def main():
    ##### sys.argv == [categorize_parameters.txt, *_index.out]
    # sys.argv == [categorize_parameters.txt, modelfile]
    #if(len(sys.argv) <= 2): sys.exit("\nERROR! Fix your inputs!\n\nArg 1:  input param file detailing each voronoi 'structure'.\nShould be of the form:\nCrystal:\n    0,2,8,*\n\nArg2: an _index.out file.\n\nOutput is printed to screen.")
    if(len(sys.argv) <= 2): sys.exit("\nERROR! Fix your inputs!\n\nArg 1:  input param file detailing each voronoi 'structure'.\nShould be of the form:\nCrystal:\n    0,2,8,*\n\nArg2: a model file.\n\nOutput is printed to screen.")

    paramfile = sys.argv[1]
    modelfile = sys.argv[2]

    m = Model(modelfile)

    #voronoi_3d(m,3.5)

    vorrun = vor.Vor()
    vorrun.runall(modelfile,3.5)
    vorrun.set_atom_vp_indexes(m)

    vp_dict = load_param_file(paramfile)
    #atom_dict = generate_atom_dict(vorrun.index,vp_dict)
    #vor_cats.load_index_file(sys.argv[2])
    #printVorCats(atom_dict,vp_dict)

    set_atom_vp_types(m,vp_dict)

    index = 0
    m.generate_neighbors(3.5)
    for atom in m.atoms:
        if(atom.vp.index[0:4] == [0, 0, 12, 0]):
            temp_model = Model("comment", m.lx, m.ly, m.lz, atom.neighs+[atom])
            temp_model.write_cif("temp{0}.cif".format(index))
            index += 1
            print("temp{0}.cif".format(index))

    #print("Undefined indexes:")
    #for atom in m.atoms:
    #    if atom.vp.type == "Undef":
    #        print("{0} {1} {2}".format(atom.id,atom.z,atom.vp.index))

    #print_all(m)
    #vor_stats(m) # Prints what you probably want


if __name__ == "__main__":
    main()
