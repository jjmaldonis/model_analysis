#import matplotlib
#matplotlib.use('PDF')
import copy
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Voronoi, voronoi_plot_2d
from pprint import pprint
from atom import Atom
from hutch import Hutch
import znum2sym
import sys
import math
import masses

class Model(object):
    """ xyz model file class 
    functions:
        write_our_xyz(outfile)
        write_cif(outfile)
        generate_neighbors(cutoff)
        get_atoms_in_cutoff(atom,cutoff)
        nearest_neigh(atom)
        save_vp_dict(vp_dict) """
    
    def __init__(self, *args, **kwargs):
        """ sets:
                self.comment
                self.lx
                self.ly
                self.lz
                self.atoms
                self.natoms
                self.atomtypes
                self.xx, self.yy, self.zz
            Extra routines can set:
                self.coord_numbers """

        super(Model,self).__init__()
        if(len(args) == 1):
            modelfile = args[0]
            if modelfile[-4:] == '.xyz':
                self.read_xyz(modelfile)
            elif modelfile[-4:] == '.dat':
                self.read_dat(modelfile)
            else:
                raise Exception("Unknown input model type! File={0}".format(modelfile))
        elif(len(args) == 5):
            self.comment = args[0]
            self.lx = args[1]
            self.ly = args[2]
            self.lz = args[3]
            self.atoms = [atom.copy() for atom in args[4]]
            self.natoms = len(args[4])
        else:
            raise Exception("Unknown input parameters to Model()!")
        self.hutch = Hutch(self)

        self.atomtypes = {}
        for atom in self.atoms:
            self.atomtypes[atom.z] = self.atomtypes.get(atom.z, 0) + 1

        self.generate_position_arrays()


    def read_xyz(self,modelfile):
        with open(modelfile) as f:
            content = f.readlines()

        self.comment = content.pop(0) # Comment line
        if('-1' in content[-1] and '.' not in content[-1]): # '-1' line
            content.pop(-1)
        self.lx,self.ly,self.lz = tuple([float(x) for x in content.pop(0).strip().split()])

        self.natoms = len(content)
        content = [x.strip().split() for x in content]
        for i in range(0,len(content)):
            for j in range(0,len(content[i])):
                try:
                    content[i][j] = int(content[i][j])
                except:
                    try:
                        content[i][j] = float(content[i][j])
                    except:
                        pass
        self.atoms = []
        for i,atom in enumerate(content):
            self.atoms.append(Atom(i,atom[0],atom[1],atom[2],atom[3]))
            #if type(self.atoms[i].z) == type('hi'):
            #    self.atoms[i].z = znum2sym.sym2z(self.atoms[i].z)

    def read_dat(self,modelfile):
        with open(modelfile) as f:
            content = f.readlines()
        self.comment = content.pop(0) # Comment line
        content = [x for x in content if not x.startswith('#')]

        for line in content:
            if('atoms' in line): self.natoms = int(line.split()[0])
            if('xlo' in line and 'xhi' in line):
                self.lx = abs(float(line.split()[0])) + abs(float(line.split()[1]))
            if('ylo' in line and 'yhi' in line):
                self.ly = abs(float(line.split()[0])) + abs(float(line.split()[1]))
            if('zlo' in line and 'zhi' in line):
                self.lz = abs(float(line.split()[0])) + abs(float(line.split()[1]))
            if('atom types' in line): nelems = int(line.split()[0])
            if('Masses' in line): mflag = content.index(line) + 1
            if('Atoms' in line): aflag = content.index(line) + 1
        try:
            mflag
        except NameError:
            raise Exception("ERROR! You need to define the masses in the .dat file.")
        atomtypes = {}
        while(nelems > 0):
            if(len(content[mflag].split()) == 2):
                atomtypes[int(content[mflag].split()[0])] = masses.get_znum(float(content[mflag].split()[1]))
                nelems -= 1
            mflag += 1
        self.atoms = []
        natoms = self.natoms
        while(natoms > 0):
            sline = content[aflag].split()
            if(len(sline) >= 5):
                # We found an atom
                id = int(sline[0])
                type = int(sline[1])
                x = float(sline[2])
                y = float(sline[3])
                z = float(sline[4])
                znum = atomtypes[type]
                # Add it to the model
                self.atoms.append(Atom(id,znum,x,y,z))
                natoms -= 1
            aflag += 1

    def write_dat(self,outfile=None):
        if outfile is None:
            of = sys.stdout
        else:
            of = open(outfile,'w')
        of.write(self.comment+'\n')
        of.write('{0} atoms\n\n'.format(self.natoms))
        of.write('{0} atom types\n\n'.format(len(self.atomtypes)))
        of.write('{0} {1} xlo xhi\n'.format(-self.lx/2,self.lx/2))
        of.write('{0} {1} ylo yhi\n'.format(-self.ly/2,self.ly/2))
        of.write('{0} {1} zlo zhi\n\n'.format(-self.lz/2,self.lz/2))
        of.write('Masses\n\n')
        atomtypes = list(self.atomtypes)
        atomtypes.sort()
        atomtypes.reverse()
        for i,z in enumerate(atomtypes):
            of.write('{0} {1}\n'.format(i+1,round(masses.get_mass(z),2)))
        of.write('\n')
        of.write('Atoms\n\n')
        for i,atom in enumerate(self.atoms):
            of.write('{0} {1} {2} {3} {4}\n'.format(atom.id+1, atomtypes.index(atom.z)+1, atom.coord[0], atom.coord[1], atom.coord[2]))


    def write_our_xyz(self,outfile=None):
        if outfile is None:
            of = sys.stdout
        else:
            of = open(outfile,'w')
        of.write(self.comment.strip()+'\n')
        of.write("{0} {1} {2}\n".format(self.lx, self.lx, self.lz))
        for atom in self.atoms:
            of.write(atom.ourxyz()+'\n')
        of.write('-1')
        of.close()

    def write_real_xyz(self,outfile=None):
        if outfile is None:
            of = sys.stdout
        else:
            of = open(outfile,'w')
        of.write(str(self.natoms)+'\n')
        of.write("{1} {2} {3} {0}\n".format(self.comment.strip(),self.lx, self.lx, self.lz))
        for i,atom in enumerate(self.atoms):
            of.write(atom.realxyz()+'\n')

    def write_cif(self,outfile=None):
        if outfile is None:
            of = sys.stdout
        else:
            of = open(outfile,'w')
        of.write('_pd_phase_name\t'+self.comment+'\n')
        of.write('_cell_length_a '+str(self.lx)+'\n')
        of.write('_cell_length_b '+str(self.ly)+'\n')
        of.write('_cell_length_c '+str(self.lz)+'\n')
        of.write('_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n')
        of.write('_symmetry_space_group_name_H-M         \'P 1\'\n_symmetry_Int_Tables_number            1\n\n')
        of.write('loop_\n_symmetry_equiv_pos_as_xyz\n   \'x, y, z\'\n\n')
        of.write('loop_\n   _atom_site_label\n   _atom_site_occupancy\n   _atom_site_fract_x\n   _atom_site_fract_y\n   _atom_site_fract_z\n   _atom_site_adp_type\n   _atom_site_B_iso_or_equiv\n   _atom_site_type_symbol\n')

        atomtypes = {}
        for atom in self.atoms:
            atomtypes[atom.z] = atomtypes.get(atom.z, 0) + 1
            of.write('   '+znum2sym.z2sym(atom.z)+str(atomtypes[atom.z])+'\t1.0\t'+str(atom.coord[0]/self.lx+0.5)+'\t'+str(atom.coord[1]/self.ly+0.5)+'\t'+str(atom.coord[2]/self.lz+0.5)+'\tBiso\t1.000000\t'+znum2sym.z2sym(atom.z)+'\n')
        of.write('\n')
        of.close()

    
        
    def add_atom(self,atom):
        # This doesn't update anything besides the hutch and self.atoms
        self.atoms.append(atom)
        self.hutch.add_atom(atom)
        self.natoms += 1

    def remove_atom(self,atom):
        # This doesn't update anything besides the hutch and self.atoms
        self.hutch.remove_atom(atom)
        self.atoms.remove(atom)
        self.natoms -= 1

    def generate_neighbors(self,cutoff):
        for atom in self.atoms:
            #if(atom.z == 40):
            #    atom.neighs = self.get_atoms_in_cutoff(atom,3.82)
            #else:
            #    atom.neighs = self.get_atoms_in_cutoff(atom,cutoff)
            atom.neighs = []
            atom.neighs = self.get_atoms_in_cutoff(atom,cutoff)
            if(atom in atom.neighs): atom.neighs.remove(atom)
            atom.cn = len(atom.neighs)

    def check_neighbors(self):
        for atom in self.atoms:
            for n in atom.neighs:
                if atom not in n.neighs:
                    print("You're neighbors are screwed up! Atom IDs are {0}, {1}.".format(atom.id,n.id))
                    print("Neighbors of: {0}".format(atom))
                    print(atom.neighs)
                    print("Neighbors of: {0}".format(n))
                    print(n.neighs)
                    print("Dist = {0}".format(self.dist(atom,n)))
                    print("Dist = {0}".format(self.dist(n,atom)))


    def generate_coord_numbers(self):
        """ atom.neighs must be defined first for all atoms """
        self.coord_numbers = {}
        # Form will be:
        #   {'Cu-Al': 4.5}, etc.
        for typea in self.atomtypes:
            self.coord_numbers[znum2sym.z2sym(typea)] = 0
            for typeb in self.atomtypes:
                self.coord_numbers[znum2sym.z2sym(typea)+'-'+znum2sym.z2sym(typeb)] = 0
        for atom in self.atoms:
            for n in atom.neighs:
                self.coord_numbers[znum2sym.z2sym(atom.z)] += 1
                self.coord_numbers[znum2sym.z2sym(atom.z)+'-'+znum2sym.z2sym(n.z)] += 1
        self.bonds = self.coord_numbers.copy()
        for key in self.coord_numbers:
            elem = znum2sym.sym2z(key.split('-')[0])
            self.coord_numbers[key] /= float(self.atomtypes[elem])


    def get_atoms(self):
        return self.atoms
    def get_natoms(self):
        return self.natoms
    def get_box_size(self):
        return (self.lx,self.ly,self.lz)

    def get_atoms_in_cutoff(self,atom,cutoff):
        x_start = atom.coord[0] - cutoff
        y_start = atom.coord[1] - cutoff
        z_start = atom.coord[2] - cutoff
        x_end   = atom.coord[0] + cutoff
        y_end   = atom.coord[1] + cutoff
        z_end   = atom.coord[2] + cutoff
        if(x_start < -self.lx/2.0): x_start = x_start + self.lx #PBC
        if(y_start < -self.ly/2.0): y_start = y_start + self.ly #PBC
        if(z_start < -self.lz/2.0): z_start = z_start + self.lz #PBC
        if(x_end > self.lx/2.0): x_end = x_end - self.lx #PBC
        if(y_end > self.ly/2.0): y_end = y_end - self.ly #PBC
        if(z_end > self.lz/2.0): z_end = z_end - self.lz #PBC
        atom_s = Atom(0,0,x_start,y_start,z_start)
        hutch_s = self.hutch.hutch_position(atom_s)
        i_start = hutch_s[0]
        j_start = hutch_s[1]
        k_start = hutch_s[2]
        atom_e = Atom(0,0,x_end,y_end,z_end)
        hutch_e = self.hutch.hutch_position(atom_e)
        i_end = hutch_e[0]
        j_end = hutch_e[1]
        k_end = hutch_e[2]
        list = []

        for i in range(0,self.hutch.nhutchs):
            if(i_start <= i_end):
                if(i < i_start or i > i_end): continue
            else:
                if(i < i_start and i > i_end): continue
            for j in range(0,self.hutch.nhutchs):
                if(j_start <= j_end):
                    if(j < j_start or j > j_end): continue
                else:
                    if(j < j_start and j > j_end): continue
                for k in range(0,self.hutch.nhutchs):
                    if(k_start <= k_end):
                        if(k < k_start or k > k_end): continue
                    else:
                        if(k < k_start and k > k_end): continue
                    list = list + self.hutch.get_atoms_in_hutch((i,j,k))
        #if(atom in list): list.remove(atom)
        list = [atomi for atomi in list if self.dist(atom,atomi) <= cutoff]
        return list

    def dist(self,atom1,atom2):
        x = (atom1.coord[0] - atom2.coord[0])
        y = (atom1.coord[1] - atom2.coord[1])
        z = (atom1.coord[2] - atom2.coord[2])
        x = x - self.lx*round(x/self.lx)
        y = y - self.ly*round(y/self.ly)
        z = z - self.lz*round(z/self.lz)
        #if(math.sqrt(x**2+y**2+z**2) > self.lx/2.0*math.sqrt(3)):
        #    print("Error in dist")
        #    print(atom1,atom2)
        #    print(x,y,z,math.sqrt(x**2+y**2+z**2))
        #while(x > self.lx/2): x = self.lx - x
        #while(y > self.ly/2): y = self.ly - y
        #while(z > self.lz/2): z = self.lz - z
        #while(x < -self.lx/2): x = self.lx + x
        #while(y < -self.ly/2): y = self.ly + y
        #while(z < -self.lz/2): z = self.lz + z
        return math.sqrt(x**2+y**2+z**2)
    def get_all_dists(self):
        dists = []
        for atomi in self.atoms:
            for atomj in self.atoms[self.atoms.index(atomi)+1:]:
                dists.append([atomi,atomj,self.dist(atomi,atomj)])
        return dists
    def nearest_neigh_of_same_type(self,atom):
        """ returns the nearest at of the same species as 'atom' """
        cutoff = 3.5
        atoms = []
        while(len(atoms) == 0):
            atoms = self.get_atoms_in_cutoff(atom,cutoff)
            atoms = [x for x in atoms if x.z == atom.z]
            #if atom in atoms: atoms.remove(atom)
            cutoff *= 2
        cutoff /= 2 # set back to the value used in case I want it later
        d = float("inf")
        for atomi in atoms:
            dt = self.dist(atom,atomi)
            if dt < d:
                d = dt
                a = atomi
        if(a.z != atom.z): raise Exception("Error! Function 'nearest_neigh_of_same_type' didn't work!")
        return a
    def nearest_neigh(self,atom):
        """ returns an atoms nearest neighbor """
        hutch = self.hutch.hutch_position(atom)
        atoms = self.hutch.get_atoms_in_hutch(hutch)[:]
        if atom in atoms: atoms.remove(atom)

        # This generation of nearby hutches isn't perfect but it will work
        rots = [(1,0,0),(0,1,0),(0,0,1)]
        i = 0
        while len(atoms) == 0:
            hutch = ((hutch[0]+rots[i][0])%self.hutch.nhutchs,(hutch[1]+rots[i][1])%self.hutch.nhutchs,(hutch[2]+rots[i][2])%self.hutch.nhutchs)
            i = (i+1) % 3
            atoms = self.hutch.get_atoms_in_hutch(hutch)
            if atom in atoms: atoms.remove(atom)
        start = atoms[0]

        atoms = self.get_atoms_in_cutoff(atom,self.dist(atom,start))
        #if atom in atoms: atoms.remove(atom)
        d = float("inf")
        for atomi in atoms:
            dt = self.dist(atom,atomi)
            if dt < d:
                d = dt
                a = atomi
        return a

    def save_vp_dict(self,vp_dict):
        self.vp_dict = vp_dict

    def composition(self):
        d = {}
        for key in self.atomtypes:
            d[znum2sym.z2sym(key)] = self.atomtypes[key]/float(self.natoms)
        return d

    def print_bond_stats(self):
        comp = self.composition()
        nbonds = 0.0
        for key in self.bonds:
            if '-' not in key:
                nbonds += self.bonds[key]
        bond_stats = self.bonds.copy()
        for key in bond_stats.keys():
            if '-' in key:
                del bond_stats[key]
            else:
                #bond_stats[key] /= nbonds/comp[znum2sym.sym2z(key)]
                bond_stats[key] *= 1.0/(nbonds*comp[key])
        print('Bond statistics:')
        pprint(self.bonds)
        print('Ratio of actual/expected bonds:')
        pprint(bond_stats)

    def atom_coords(self):
        coords = []
        for atom in self.atoms:
            coords.append(atom.coord)
        return coords

    def nn_vor(self):
        coords = self.atom_coords()
        vor = Voronoi(coords)
        print(vor.vertices)
        print(vor.regions)
        voronoi_plot_2d(vor)
        plt.show()

    def generate_position_arrays(self):
        self.xx = []
        self.yy = []
        self.zz = []
        for i,atomi in enumerate(self.atoms):
            self.xx.append(atomi.coord[0])
            self.yy.append(atomi.coord[1])
            self.zz.append(atomi.coord[2])
        self.xx = np.array(self.xx)
        self.yy = np.array(self.yy)
        self.zz = np.array(self.zz)





def main():
    m = Model(sys.argv[1])
    outflag = False
    if(len(sys.argv) > 3):
        if(sys.argv[2] == '-o'):
            outtype = sys.argv[3]
            outflag = True
    if(outflag):
        if(outtype == 'dat' or outtype == '.dat'):
            m.write_dat()
        elif(outtype == 'cif' or outtype == '.cif'):
            m.write_cif()
        elif(outtype == 'xyz' or outtype == '.xyz'):
            m.write_our_xyz()
        elif(outtype == 'realxyz' or outtype == '.realxyz'):
            m.write_real_xyz()

    ## create histogram of distances
    #dists = np.array(m.get_all_dists())
    #dists = [dist[2] for dist in dists]
    #hist,bin = np.histogram(dists,200)
    #hist = list(hist)
    #bin = list(bin)
    #bin.pop()
    #print(hist)
    #print(bin)
    #print(hist.index(2172))
    #print(bin[hist.index(2172)])
    #print(hist.index(2053))
    #print(bin[hist.index(2053)])
    #plt.plot(bin,hist)
    #plt.savefig('temp.png',format='png')

    #m.nn_vor()

    #m.write_dat()
    #m.write_cif(sys.argv[1][:-3]+'cif')
    #m.write_our_xyz(sys.argv[1][:sys.argv[1].rindex('.')]+'.our.xyz')

    #dists = []
    #for atom in m.atoms:
    #    dists.append(m.dist(atom,m.nearest_neigh(atom)))
    #print(sum(dists)/len(dists))
    #print(max(dists))
    #print(min(dists))

if __name__ == "__main__":
    main()


