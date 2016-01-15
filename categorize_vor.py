
# Categorizes the voronoi polyhedra according to the input parameter file, 
# which must be the first arg passed to this script.
# The second arg must be a _index.out file.
# Output is printed to screen.

import sys, math, copy
from collections import OrderedDict
from znum2sym import z2sym
from model import Model
from voronoi_3d import voronoi_3d

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
    #with open(paramfile) as f:
    #    vp_params = [line.split(',') for line in f]
    vp_params = open(paramfile).readlines()
    # For each voronoi polyhedra 'structure', change it to an int if it's not a *
    vp_dict = {}
    for line in vp_params:
        if(':' in line): # Title line
            current_entry = line.strip()[:-1]
            vp_dict[current_entry] = []
        else: # Index line
            vps = line.strip().split(',')
            vps = [x.strip() for x in vps] # remove spaces if ", " rather than just ","
            for i in range(len(vps)):
                if(vps[i] != '*'):
                    vps[i] = int(vps[i])
            vp_dict[current_entry].append(vps)
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


class VPStatistics(object):
    def __init__(self, model):
        self.m = model
        self.natoms = self.m.natoms
        self.indices = self.index_stats()
        self.indexes = self.indices
        self.categories = self.vor_stats()

    def vor_stats(self):
        cats = {}
        for atom in self.m.atoms:
            if atom.vp.type not in cats:
                cats[atom.vp.type] = {}
            cats[atom.vp.type][atom.sym] = cats[atom.vp.type].get(atom.sym,0) + 1
            cats[atom.vp.type]["Total"] = cats[atom.vp.type].get("Total",0) + 1
        return cats

    def print_categories(self):
        """ Prints the number of atoms in each VP category """
        for key in sorted(self.categories):
            print("{0}:\nTotal:\t{1}\t{2}%".format(key, self.categories[key]["Total"], round(100.0*self.categories[key]["Total"]/self.natoms,2)))
            for elem in sorted(self.categories[key]):
                if(elem != "Total"):
                    print("{0}:\t{1}\t{2}%".format(elem, self.categories[key][elem],
                        round(100.0*self.categories[key][elem]/self.categories[key]['Total'],2)))

    def index_stats(self):
        cats = {}
        for atom in self.m.atoms:
            if atom.vp.index not in cats:
                cats[atom.vp.index] = OrderedDict()
                for typ in sorted([z2sym(x) for x in self.m.atomtypes], reverse=True):
                    cats[atom.vp.index][typ] = 0
            cats[atom.vp.index][atom.sym] = cats[atom.vp.index].get(atom.sym,0) + 1
            cats[atom.vp.index]["Total"]  = cats[atom.vp.index].get("Total",0) + 1
        return cats

    def print_indices(self, cutoff=0.005):
        """ Prints the number of atoms in each VP index"""
        for val,key in sorted( ((v,k) for k,v in self.indices.iteritems()), key=lambda t: t[0]['Total']): 
            if val['Total'] < self.natoms*cutoff: continue
            output = copy.copy(val)
            for k,v in output.items():
                if k != 'Total':
                    output[k] = '{0}  ({1}%) '.format(v, int(round(v*100./output['Total'])))
            output = [(k,outputue) for k,outputue in output.items()]
            output = [str(x) for row in output for x in row]
            output = '\t'.join(output)
            print("{0}:\tCN: {1}  \t{2}".format(key, sum(key), output))

    def print_indexes(self):
        self.print_indices()

    def __add__(self, other):
        if isinstance(other, Model):
            self.natoms += other.natoms
            cats = self.indices
            for atom in other.atoms:
                if atom.vp.index not in cats:
                    cats[atom.vp.index] = OrderedDict()
                    for typ in sorted([z2sym(x) for x in other.atomtypes], reverse=True):
                        cats[atom.vp.index][typ] = 0
                cats[atom.vp.index][atom.sym] = cats[atom.vp.index].get(atom.sym,0) + 1
                cats[atom.vp.index]["Total"]  = cats[atom.vp.index].get("Total",0) + 1

            cats = self.categories
            for atom in other.atoms:
                if atom.vp.type not in cats:
                    cats[atom.vp.type] = {}
                cats[atom.vp.type][atom.sym] = cats[atom.vp.type].get(atom.sym,0) + 1
                cats[atom.vp.type]["Total"] = cats[atom.vp.type].get("Total",0) + 1
        #elif isinstance(other, VPStatistics):
        else:
            raise Exception("Cannot add type {0} to VPStatistics!".format(type(other)))
        return self


def print_all(m):
    """ Prints the index and type of each atom in m """
    for atom in m.atoms:
        print("{0} {1} {2}".format(atom, atom.vp.index, atom.vp.type))

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


def fix_cluster_pbcs(m):
    # First recenter to first octant
    m.recenter()
    meanx, meany, meanz = m.lx/4.0, m.ly/4.0, m.lz/4.0
    for atom in m.atoms[1:]:
        atom.coord = (atom.coord[0]+meanx-m.atoms[0].coord[0], atom.coord[1]+meany-m.atoms[0].coord[1], atom.coord[2]+meanz-m.atoms[0].coord[2])
    m.atoms[0].coord = (m.atoms[0].coord[0]+meanx-m.atoms[0].coord[0], m.atoms[0].coord[1]+meany-m.atoms[0].coord[1], m.atoms[0].coord[2]+meanz-m.atoms[0].coord[2])

    # See if we need to fix
    fix = False
    for atom in m.atoms[1:]:
        if round(m.dist(m.atoms[0], atom)) != round(m.dist(m.atoms[0], atom, pbc=False)):
            fix = True
            break
    else:
        m.recenter()
        return m
    # If so, fix
    for atom in m.atoms:
        new = []
        if round(m.dist(m.atoms[0], atom)) != round(m.dist(m.atoms[0], atom, pbc=False)):
            for c in atom.coord:
                if c < 0:
                    new.append(c+m.lx)
                elif c > m.lx:
                    new.append(c-m.lx)
                else:
                    new.append(c)
            atom.coord = tuple(new)
    m.recenter()
    return m


def main():
    # sys.argv == [categorize_parameters.txt, modelfile]
    if(len(sys.argv) <= 2): sys.exit("\nERROR! Fix your inputs!\n\nArg 1:  input param file detailing each voronoi 'structure'.\nShould be of the form:\nCrystal:\n    0,2,8,*\n\nArg2: a model file.\n\nOutput is printed to screen.")

    paramfile = sys.argv[1]
    modelfiles = sys.argv[2:]

    from cutoff import cutoff

    vp_dict = load_param_file(paramfile)

    m0 = Model(modelfiles[0])
    m0.generate_neighbors(cutoff)
    voronoi_3d(m0, cutoff)
    set_atom_vp_types(m0, vp_dict)
    stats0 = VPStatistics(m0)
    print(modelfiles[0])
    #stats0.print_indexes()
    stats0.print_categories()
    return

    if len(modelfiles) > 1:
        for modelfile in modelfiles[1:]:
            print(modelfile)
            m = Model(modelfile)
            voronoi_3d(m, cutoff)
            set_atom_vp_types(m, vp_dict)
            stats = VPStatistics(m)
            stats.print_categories()

            #stats0 = stats0 + m
            #stats0.print_indexes()
            #stats0.print_categories()




if __name__ == "__main__":
    main()
