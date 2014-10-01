
# Categorizes the voronoi polyhedra according to the input parameter file, 
# which must be the first arg passed to this script.
# The second arg must be a _index.out file.
# Output is printed to screen.

import sys
import vor
from model import Model
from voronoi_3d import voronoi_3d
from vor import Vor, fortran_voronoi_3d

""" Functions:
load_index_file(indexfile)
load_param_file(paramfile)
generate_atom_dict(indexes,vp_dict)
categorize_atoms(m,paramfile)
categorize_index(ind, vp_dict)
set_atom_vp_types(model,vp_dict)
vor_stats(m)
print_all(m)
save_vp_cluster_with_index(m,index) """


def load_index_file(indexfile):
    """ Simply reads in an index file into a list and returns it """
    with open(indexfile) as f:
        index = f.readlines()
    return index

def load_param_file(paramfile):
    """ Returns a dictionary in the form:
        {'Crystal-like': [[0, 4, 4, '*'], [0, 5, 2, '*']], 
        'Icosahedra-like': [[0, 2, 8, '*'], [0, 1, 10, '*'], [0, 0, 12, 0], [0, 0, 12, 2]],
        'Mixed': [[0, 3, 6, '*'], [0, 3, 7, 2], [1, 2, 5, 4]]} """
    # Open the input parameter file. Should be of the form:
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
            # If there is a : at the end of the vps, get rid of it
            if(':' == vps[0][-1]):
                current_entry = vps[0][:-1]
            else:
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

def generate_atom_dict(model):
    """ Generate a new dictionary to store all the atoms that are crystalline, icosahedra, etc.
    Returns a dictionary in the form:
    { 'Mixed:': [list of mixed atoms], 'Crystal-like:', [list of crystal-like atoms], etc}.
    All atoms must be assigned a VP type prior to this function. """
    atom_dict = {}
    for atom in model.atoms:
        atom_dict[atom.vp.type] = atom_dict.get(atom.vp.type,[]) + [atom]
    return atom_dict


def categorize_atoms(m,paramfile):
    """ Shortcut to run load_param_file
    and set_atom_vp_types in one function.
    Also stores vp_dict in the model. """
    vp_dict = load_param_file(paramfile)
    set_atom_vp_types(m,vp_dict)
    for key in vp_dict.keys():
        vp_dict[key[:key.index(':')]] = vp_dict.pop(key)
    m.vp_dict = vp_dict

def categorize_index(ind, vp_dict):
    """ ind should be an integer list
    vp_dict is returned by load_param_file.
    This function returns the type of index
    that ind is given vp_dict. """
    for key in vp_dict:
        for vps in vp_dict[key]:
            found = True
            for i in range(0,4):
                if(vps[i] != '*' and vps[i] != ind[i]):
                    found = False
            if(found):
                return key
    return 'Undef'

def set_atom_vp_types(model,vp_dict):
    """ saves the voronoi polyhedra type for each atom to the atom in the model """
    for atom in model.atoms:
        atom.vp.type = categorize_index(atom.vp.index,vp_dict)


def vor_stats(m):
    """ Prints the number of atoms in each VP category """
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
    return cats

def print_all(m):
    """ Prints the index and type of each atom in m """
    for atom in m.atoms:
        print("{0} {1} {2}".format(atom,atom.vp.index,atom.vp.type))

def save_vp_cluster_with_index(m,index):
    """ Index should be a 4-list, e.g. [0,0,12,0].
    This function goes thru the model and finds all
    atoms with index "index" and saves that atom's
    VP as a new model, with name "temp{atom.id}.cif """
    for atom in m.atoms:
        if(atom.vp.index[0:4] == index):
            temp_model = Model("VP with index {0}".format(index), m.lx, m.ly, m.lz, atom.neighs+[atom])
            temp_model.write_cif("temp{0}.cif".format(atom.id))
            print("Saved VP cluster to modelfile temp{0}.cif".format(atom.id))
        


def main():
    # sys.argv == [categorize_parameters.txt, modelfile]
    if(len(sys.argv) <= 2): sys.exit("\nERROR! Fix your inputs!\n\nArg 1:  input param file detailing each voronoi 'structure'.\nShould be of the form:\nCrystal:\n    0,2,8,*\n\nArg2: a model file.\n\nOutput is printed to screen.")

    paramfile = sys.argv[1]
    modelfile = sys.argv[2]

    m = Model(modelfile)

    cutoff = {}
    #cutoff[(40,40)] = 3.5
    #cutoff[(13,29)] = 3.5
    #cutoff[(29,13)] = 3.5
    #cutoff[(40,13)] = 3.5
    #cutoff[(13,40)] = 3.5
    #cutoff[(29,40)] = 3.5
    #cutoff[(40,29)] = 3.5
    #cutoff[(13,13)] = 3.5
    #cutoff[(29,29)] = 3.5
    cutoff[(41,41)] = 3.7
    cutoff[(28,28)] = 3.7
    cutoff[(41,28)] = 3.7
    cutoff[(28,41)] = 3.7

    voronoi_3d(m,cutoff)
    #m = fortran_voronoi_3d(modelfile,3.5)

    #vorrun = vor.Vor()
    #m = vorrun.runall(modelfile,3.5)
    #vorrun.set_atom_vp_indexes(m)

    vp_dict = load_param_file(paramfile)
    set_atom_vp_types(m,vp_dict)

    #atom_dict = generate_atom_dict(vorrun.index,vp_dict)
    #vor_cats.load_index_file(sys.argv[2])
    #printVorCats(atom_dict,vp_dict)

    #m.generate_neighbors(3.5)
    #save_vp_cluster_with_index(m,[0,0,12,0])

    #cutoff = 3.5
    #m2 = Model(sys.argv[3])
    #vor_instance = Vor()
    #vor_instance.runall(modelfile,cutoff)
    #vor_instance.set_atom_vp_indexes(m)
    #nbins = 6
    #del_bin = 100.0/nbins
    #fracs = []
    #for atom in m2.atoms:
    #    fracs.append((float(m.atoms[m.atoms.index(atom)].vp.index[2])/float(sum(m.atoms[m.atoms.index(atom)].vp.index))*100.0))
    #    bin =    int((float(m.atoms[m.atoms.index(atom)].vp.index[2])/float(sum(m.atoms[m.atoms.index(atom)].vp.index))*100.0) /(100.0/(nbins-1)))
    #    atom.z = bin+1
    #fracs.sort()
    #print('Min %: {0}. Max %: {1}'.format(min(fracs),max(fracs)))
    #for atom in m2.atoms:
    #    atom.vp.type = m.atoms[m.atoms.index(atom)].vp.type
    #atoms = []
    #for atom in m2.atoms:
    #    if(atom.vp.type == "Icosahedra-like"):
    #        atoms.append(atom)
    #m3 = Model(str(len(atoms)),m.lx,m.ly,m.lz,atoms)
    #m3.write_real_xyz()


    #print_all(m)
    vor_stats(m) # Prints what you probably want
    for atom in m.atoms:
        if( atom.vp.type == "Undef"):
            print( atom.vp.index )
    print( "")
    for atom in m.atoms:
         if( atom.vp.type == "Crystal-like"):
             print(atom.vp.index )
    print (m.atoms[163].vp.index)
    print (m.atoms[163].vp.neighs)
    print (m.atoms[163].coord)
    print (m.atoms[163].id)

if __name__ == "__main__":
    main()
