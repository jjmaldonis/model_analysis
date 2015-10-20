
# Categorizes the voronoi polyhedra according to the input parameter file, 
# which must be the first arg passed to this script.
# The second arg must be a _index.out file.
# Output is printed to screen.

import sys, math, copy, os
from collections import OrderedDict, defaultdict
from znum2sym import z2sym
import vor
from model import Model
from voronoi_3d import voronoi_3d
from vor import Vor, fortran_voronoi_3d
from recenter_model import recenter_model
from nearest_atoms import find_center_atom

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

def set_atom_vp_types(model, vp_dict, submodel=None):
    """ saves the voronoi polyhedra type for each atom to the atom in the model """
    if submodel is None:
        for atom in model.atoms:
            atom.vp.type = categorize_index(atom.vp.index,vp_dict)
    else:
        for count, subatom in enumerate(submodel.atoms):
            i = model.atoms.index(subatom)
            atom = model.atoms[i]
            subatom.vp.type = atom.vp.type


def vor_stats(m, verbose=True):
    """ Prints the number of atoms in each VP category """
    cats = {}
    for atom in m.atoms:
        if atom.vp.type not in cats:
            cats[atom.vp.type] = {}
        cats[atom.vp.type][atom.sym] = cats[atom.vp.type].get(atom.sym,0) + 1
        cats[atom.vp.type]["Total"] = cats[atom.vp.type].get("Total",0) + 1
    if not verbose: return cats
    # Print
    for key in sorted(cats):
        print("{0}:\nTotal:\t{1}\t{2}%".format(key,cats[key]["Total"], round(100.0*cats[key]["Total"]/m.natoms,2)))
        for elem in sorted(cats[key]):
            #if(elem != "Total"): print("   {0}: {1}".format(elem,cats[key][elem]))
            if(elem != "Total"):
                print("{0}:\t{1}\t{2}%".format(elem, cats[key][elem],
                    round(100.0*cats[key][elem]/cats[key]['Total'],2)))
    return cats

def index_stats(m):
    """ Prints the number of atoms in each VP index"""
    cats = {}
    for atom in m.atoms:
        #if atom.sym == 'Zr' or atom.sym == 'Al': continue
        #if atom.sym == 'Zr' or atom.sym == 'Cu': continue
        if atom.vp.index not in cats:
            cats[atom.vp.index] = OrderedDict()
            for typ in sorted([z2sym(x) for x in m.atomtypes], reverse=True):
                cats[atom.vp.index][typ] = 0
        cats[atom.vp.index][atom.sym] = cats[atom.vp.index].get(atom.sym,0) + 1
        cats[atom.vp.index]["Total"]  = cats[atom.vp.index].get("Total",0) + 1
    # Print
    for val,key in sorted( ((v,k) for k,v in cats.iteritems()), key=lambda t: t[0]['Total']): 
        #for atom in m.atoms:
        #    if(atom.vp.index == key):
        #        typ = atom.vp.type
        #        break
        #print("{0}: \t{1}\t{2}".format(key,val,typ))
        #if list(key)[0:4] != [0,2,8,0]:
        #    continue
        if val['Total'] < m.natoms*0.005: continue
        val = [(k,value) for k,value in val.items()]
        val = [str(x) for row in val for x in row]
        val = '\t'.join(val)
        print("{0}: \t{1}".format(key,val))
        #print("{0}: \t{1}\t{2}".format(key,val,round(val['Total']/float(m.natoms)*100,3)))
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


def dist(atom1,atom2):
    x = (atom1.coord[0] - atom2.coord[0])
    y = (atom1.coord[1] - atom2.coord[1])
    z = (atom1.coord[2] - atom2.coord[2])
    return math.sqrt(x**2+y**2+z**2)

def fix_cluster_pbcs(m):
    # First recenter to first octant
    recenter_model(m)
    meanx, meany, meanz = m.lx/4.0, m.ly/4.0, m.lz/4.0
    for atom in m.atoms[1:]:
        atom.coord = (atom.coord[0]+meanx-m.atoms[0].coord[0], atom.coord[1]+meany-m.atoms[0].coord[1], atom.coord[2]+meanz-m.atoms[0].coord[2])
    m.atoms[0].coord = (m.atoms[0].coord[0]+meanx-m.atoms[0].coord[0], m.atoms[0].coord[1]+meany-m.atoms[0].coord[1], m.atoms[0].coord[2]+meanz-m.atoms[0].coord[2])

    # See if we need to fix
    fix = False
    for atom in m.atoms[1:]:
        if round(m.dist(m.atoms[0], atom)) != round(dist(m.atoms[0], atom)):
            fix = True
            break
    else:
        recenter_model(m)
        return m
    # If so, fix
    for atom in m.atoms:
        new = []
        if round(m.dist(m.atoms[0], atom)) != round(dist(m.atoms[0], atom)):
            for c in atom.coord:
                if c < 0:
                    new.append(c+m.lx)
                elif c > m.lx:
                    new.append(c-m.lx)
                else:
                    new.append(c)
            atom.coord = copy.copy(new)
    recenter_model(m)
    return m

def normalize_bond_distances(m):
    """ Rescales a cluster so that the average bond length is 1.0 """
    center = find_center_atom(m)
    for atom in m.atoms:
        atom.coord = (atom.coord[0]-center.coord[0], atom.coord[1]-center.coord[1], atom.coord[2]-center.coord[2])

    avg = 0.
    for atom in m.atoms:
        if atom.id != center.id:
            avg += m.dist(center, atom)
    avg /= (m.natoms-1)
    for atom in m.atoms:
        if atom.id != center.id:
            atom.coord = (atom.coord[0]/avg, atom.coord[1]/avg, atom.coord[2]/avg)
    recenter_model(m)
    return avg


def main():
    cutoff = {}
    cutoff[(40,40)] = 3.6
    cutoff[(13,29)] = 3.6
    cutoff[(29,13)] = 3.6
    cutoff[(40,13)] = 3.6
    cutoff[(13,40)] = 3.6
    cutoff[(29,40)] = 3.6
    cutoff[(40,29)] = 3.6
    cutoff[(13,13)] = 3.6
    cutoff[(29,29)] = 3.6

    cutoff[(41,41)] = 3.7
    cutoff[(28,28)] = 3.7
    cutoff[(41,28)] = 3.7
    cutoff[(28,41)] = 3.7

    cutoff[(46,46)] = 3.45
    cutoff[(14,14)] = 3.45
    cutoff[(46,14)] = 3.45
    cutoff[(14,46)] = 3.45

    paramfile = sys.argv[1]
    vp_dict = load_param_file(paramfile)

    modelfiles = sys.argv[2:]
    count = defaultdict(int) # Stores how many VPs have been found of each index type
    count = 0
    direc = 'ZrCuAl/md_80k/'
    for modelfile in modelfiles:
        print(modelfile)
        m = Model(modelfile)
        m.generate_neighbors(cutoff)
        #voronoi_3d(m,cutoff)
        #set_atom_vp_types(m,vp_dict)
        #vor_stats(m)
        #cats = index_stats(m)
        #for atom in m.atoms:
        #    #new = Model('0,0,12,0; number of atoms is {0};'.format(count), m.lx, m.ly, m.lz, atom.neighs + [atom])
        #    new = Model('{0}; number of atoms is {1};'.format(atom.vp.index, count[atom.vp.index]), m.lx, m.ly, m.lz, atom.neighs + [atom])
        #    fix_cluster_pbcs(new)
        #    val = normalize_bond_distances(new)
        #    new.comment = '{0}; number of atoms is {1}; bond length scaling factor is {2}'.format(atom.vp.index, count,val)
        #    center = find_center_atom(new)
        #    new.remove(center)
        #    new.add(center)
        #    vp_str = ''.join([str(x) for x in atom.vp.index])
        #    if not os.path.exists(direc+vp_str):
        #        os.makedirs(direc+vp_str)
        #    new.write(direc+'{0}/{0}.{1}.xyz'.format(vp_str, count[atom.vp.index]))
        #    count[atom.vp.index] += 1
        #print(count)
        cn = 0.0
        for atom in m.atoms:
            cn += atom.cn
        cn = float(cn)/m.natoms
        print(cn)

        for atom in m.atoms:
            new_cut = copy.copy(cutoff)
            old_cn = atom.cn
            inc = 0.0
            while atom.cn < 12:
                for key,val in new_cut.items(): new_cut[key] = val + 0.1
                inc += 0.1
                atom.neighs = m.get_atoms_in_cutoff(atom, new_cut)
                if(atom in atom.neighs): atom.neighs.remove(atom)
                atom.cn = len(atom.neighs)
            new = Model('CN changed from {0} to {1};'.format(old_cn, atom.cn), m.lx, m.ly, m.lz, atom.neighs + [atom])
            new.write('temp/{0}.xyz'.format(count))
            if inc > 0.0: print("Increased shell by {0} Ang. for atom {1}".format(inc, count))
            count += 1
        cn = 0.0
        for atom in m.atoms:
            cn += atom.cn
        cn = float(cn)/m.natoms
        print(cn)

    return 0

    modelfile = sys.argv[2]
    m = Model(modelfile)
    xtal_atoms = sys.argv[3]
    xtal_atoms = Model(xtal_atoms).atoms

    #x,y,z = (round(x,6) for x in self.coord)
    #a,b,c = (round(x,6) for x in other.coord)
    #print("HERE")
    #print([round(x,7) for x in xtal_atoms[13].coord])
    #print([round(x,7) for x in m.atoms[153].coord])
    #print(type(xtal_atoms[13]))
    #print(type(m.atoms[153]))
    #print(xtal_atoms[13] == m.atoms[153])
    #return 0

    glassy_atoms = []
    for atom in m.atoms:
        if atom not in xtal_atoms:
            glassy_atoms.append(atom)
    print(len(glassy_atoms))
    print(len(xtal_atoms))
    assert len(glassy_atoms) + len(xtal_atoms) == m.natoms
    voronoi_3d(m, cutoff)
    set_atom_vp_types(m, vp_dict)
    m.generate_neighbors(3.45)
    head, tail = os.path.split(modelfile)
    head = head + '/'
    if not os.path.exists(head+'glassy/'):
        os.makedirs(head+'glassy/')
    if not os.path.exists(head+'xtal/'):
        os.makedirs(head+'xtal/')
    for count, atom in enumerate(xtal_atoms):
        i = m.atoms.index(atom)
        atom = m.atoms[i]
        new = Model('{0}'.format(atom.vp.index), m.lx, m.ly, m.lz, atom.neighs + [atom])
        fix_cluster_pbcs(new)
        val = normalize_bond_distances(new)
        new.comment = '{0}; bond length scaling factor is {1}'.format(atom.vp.index, val)
        center = find_center_atom(new)
        new.remove(center)
        new.add(center)
        vp_str = ''.join([str(x) for x in atom.vp.index])
        new.write(head+'xtal/{0}.xyz'.format(count))
    for count, atom in enumerate(glassy_atoms):
        i = m.atoms.index(atom)
        atom = m.atoms[i]
        new = Model('{0}'.format(atom.vp.index), m.lx, m.ly, m.lz, atom.neighs + [atom])
        fix_cluster_pbcs(new)
        val = normalize_bond_distances(new)
        new.comment = '{0}; bond length scaling factor is {1}'.format(atom.vp.index, val)
        center = find_center_atom(new)
        new.remove(center)
        new.add(center)
        vp_str = ''.join([str(x) for x in atom.vp.index])
        new.write(head+'glassy/{0}.xyz'.format(count))
    return 0
        




    for atom in volume_atoms.atoms:
        for i,atom2 in enumerate(m.atoms):
            if(atom.z == atom2.z and [round(x, 5) for x in atom.coord] == [round(x, 5) for x in atom2.coord]): good[i] = True
    count = defaultdict(int) # Stores how many VPs have been found of each index type
    for modelfile in modelfiles:
        print(modelfile)
        m = Model(modelfile)
        voronoi_3d(m,cutoff)
        set_atom_vp_types(m,vp_dict)
        #vor_stats(m)
        #cats = index_stats(m)
        for i,atom in enumerate(m.atoms):
            if not good[i]: continue
            #new = Model('0,0,12,0; number of atoms is {0};'.format(count), m.lx, m.ly, m.lz, atom.neighs + [atom])
            new = Model('{0}; number of atoms is {1};'.format(atom.vp.index, count), m.lx, m.ly, m.lz, atom.neighs + [atom])
            fix_cluster_pbcs(new)
            val = normalize_bond_distances(new)
            new.comment = '{0}; number of atoms is {1}; bond length scaling factor is {2}'.format(atom.vp.index, count,val)
            center = find_center_atom(new)
            new.remove(center)
            new.add(center)
            vp_str = ''.join([str(x) for x in atom.vp.index])
            if not os.path.exists(vp_str):
                os.makedirs(vp_str)
            new.write('{0}/{0}.{1}.xyz'.format(vp_str, count[atom.vp.index]))
            count[atom.vp.index] += 1
        print(count)
        for c,v in count.items():
            print("{0}: {1}".format(c,v))
        print(sum(count.values()))
        cn = 0.0
        for atom in m.atoms:
            cn += atom.cn
        cn = float(cn)/m.natoms
        print(cn)

if __name__ == "__main__":
    main()
