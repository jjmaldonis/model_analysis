#import matplotlib
#matplotlib.use('PDF')
import copy
#import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Voronoi, voronoi_plot_2d
from pprint import pprint
from atom import Atom
from hutch import Hutch
import znum2sym
import sys
import math
import masses
from collections import defaultdict

def drange(start, stop, step):
    r = start
    while r < stop:
        yield r
        r += step

def get_file_extension(f):
    return f.split('.')[-1]

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
                self.filename
            Extra routines can set:
                self.coord_numbers """

        if(len(args) == 1):
            modelfile = args[0]
            self.filename = modelfile
            self._load()
        elif(len(args) == 5):
            self.comment = args[0]
            self.lx = args[1]
            self.ly = args[2]
            self.lz = args[3]
            self.atoms = [atom.copy() for atom in args[4]]
            self.natoms = len(args[4])
            self.filename = 'None'
        else:
            raise Exception("Unknown input parameters to Model()!")
        if(self.lx and self.ly and self.lz):
            self.hutch = Hutch(self)

        if(not hasattr(self,'atomtypes')):
            self.atomtypes = {} # TODO make into a Counter
            for atom in self.atoms:
                self.atomtypes[atom.z] = self.atomtypes.get(atom.z, 0) + 1

    def __contains__(self,key):
        return key in self.atoms

    def __getitem__(self,atomId):
        try:
            if(self.atoms[atomId].id == atomId):
                return self.atoms
        except:
            pass
        for atom in self.atoms:
            if(atom.id == atomId):
                return atom
        return None

    def __len__(self):
        return self.natoms
        
    def add(self,atom):
        """ Adds atom 'atom' to the model. Note that there is no error checking to prevent
        the user from adding an atom twice. Be careful not to do that. """
        self.atoms.append(atom)
        self.natoms += 1
        self.atomtypes[atom.z] += 1
        try:
            self.hutch.add_atom(atom)
        except AttributeError:
            pass

    def remove(self,atom):
        """ Removes atom 'atom' from the model """
        try:
            self.hutch.remove_atom(atom)
        except AttributeError:
            pass
        self.atoms.remove(atom)
        self.natoms -= 1
        self.atomtypes[atom.z] -= 1

    def _load(self):
        ext = get_file_extension(self.filename)
        if(ext == 'xyz'):
            self._load_xyz()
        elif(ext == 'dat'):
            self._load_dat()
        else:
            raise Exception("Unknown input model type! File {0} has extension {1}".format(modelfile,ext))

    def _load_xyz(self):
        modelfile = self.filename
        self.atoms = []
        self.atomtypes = defaultdict(int)
        self.natoms = 0
        with open(modelfile) as f:
            natoms = int(f.readline().strip())
            comment = f.readline().strip()
            try:
                self.lx,self.ly,self.lz = tuple([float(x) for x in comment.split()[:3]])
            except:
                self.lx,self.ly,self.lz = (None,None,None)
            for i in range(natoms):
                znum,x,y,z = tuple(f.readline().strip().split())
                x = float(x); y = float(y); z = float(z)
                atom = Atom(i,znum,x,y,z)
                self.add(atom)

    def _load_dat(self):
        """ Still need to refactor """
        modelfile = self.filename
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

    def write(self,outfile=None,ext=None):
        if(outfile is not None and ext is None): ext = get_file_extension(outfile)
        elif(ext is None): ext = 'xyz'
        if(ext == 'xyz'):   self._write_xyz(outfile)
        elif(ext == 'dat'): self._write_dat(outfile)
        elif(ext == 'cif'): self._write_cif(outfile)

    def _write_dat(self,outfile=None):
        if outfile is None: of = sys.stdout
        else:               of = open(outfile,'w')
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

    def _write_xyz(self,outfile=None):
        if outfile is None: of = sys.stdout
        else:               of = open(outfile,'w')
        of.write(str(self.natoms)+'\n')
        of.write("{1} {2} {3} {0}\n".format(self.comment.strip(),self.lx, self.lx, self.lz))
        for i,atom in enumerate(self.atoms):
            of.write(atom.realxyz()+'\n')

    def _write_cif(self,outfile=None):
        if outfile is None: of = sys.stdout
        else:               of = open(outfile,'w')
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

    def set_vp_dict(self,vp_dict):
        self.vp_dict = vp_dict

    def generate_neighbors(self,cutoff):
        for atom in self.atoms:
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

    def generate_average_coord_numbers(self):
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

    def get_atoms_in_cutoff(self,atom,cutoff):
        return self.hutch.get_atoms_in_radius(atom,cutoff)

    def dist(self,atom1,atom2):
        x = (atom1.coord[0] - atom2.coord[0])
        y = (atom1.coord[1] - atom2.coord[1])
        z = (atom1.coord[2] - atom2.coord[2])
        x = x - self.lx*round(x/self.lx)
        y = y - self.ly*round(y/self.ly)
        z = z - self.lz*round(z/self.lz)
        return math.sqrt(x**2+y**2+z**2)

    def dist2(self,atom1,atom2):
        x = (atom1.coord[0] - atom2.coord[0])
        y = (atom1.coord[1] - atom2.coord[1])
        z = (atom1.coord[2] - atom2.coord[2])
        x = x - self.lx*round(x/self.lx)
        y = y - self.ly*round(y/self.ly)
        z = z - self.lz*round(z/self.lz)
        return x**2+y**2+z**2

    def get_all_dists(self):
        dists = [0. for i in range(self.natoms*self.natoms)]
        count = 0
        for atomi in self.atoms:
            for atomj in self.atoms[self.atoms.index(atomi)+1:]:
                dists[count] = (atomi, atomj, self.dist(atomi,atomj))
                count += 1
        return dists

    def nearest_neigh_of_same_type(self,atom):
        """ returns the nearest atom of the same species as 'atom' """
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
        atoms = self.hutch.get_atoms_in_same_hutch(atom)[:]
        if atom in atoms: atoms.remove(atom)

        # This generation of nearby hutches isn't perfect but it will work
        rots = [(1,0,0),(0,1,0),(0,0,1)]
        i = 0
        while len(atoms) == 0:
            hutch = ((hutch[0]+rots[i][0])%self.hutch.nhutchs,(hutch[1]+rots[i][1])%self.hutch.nhutchs,(hutch[2]+rots[i][2])%self.hutch.nhutchs)
            i = (i+1) % 3
            atoms = self.hutch.hutchs[hutch]
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

    @property
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


    def radial_composition(self, outfile):
        """ Creates 1D waves stored in outfile for each element in the model. Each 
        wave is a histogram of the number of atoms between two radial positions
        starting at the center of model and radiating outward. """
        npix = 16
        keys = self.atomtypes.keys()
        #histo = [[0.0 for x in range(npix)] for key in keys]
        histo = [{} for x in range(npix)]
        dx = (self.lx/2.0)/npix # Cube assumed
        for i,r in enumerate(drange(dx, npix*dx-dx, dx)):
            atoms = self.get_atoms_in_cutoff( (0.0,0.0,0.0), r)
            print(r, len(atoms))
            comp = {}
            for atom in atoms:
                comp[str(atom.z)] = comp.get(str(atom.z),0) + 1.0
            for type in self.atomtypes:
                if( str(type) not in comp.keys()):
                    comp[str(type)] = 0.0
            comp['Total'] = len(atoms)
            histo[i] = comp
        of = open(outfile,'w')
        of.write('IGOR\n')
        for atomtype in keys:
            of.write('\nWAVES/N=({0})\t {1}\nBEGIN\n'.format(npix,'partial_radial_comp_'+znum2sym.z2sym(atomtype)))
            for i in range(npix-2):
                if(i != 0):
                    of.write("{0} ".format((histo[i][str(atomtype)] - histo[i-1][str(atomtype)])/( 4.0/3.0*np.pi*( (i*dx)**3 - ((i-1)*dx)**3 ))))
                    #print("{1}  {0} ".format(histo[i][str(atomtype)],i*dx))
                    #print("  {1}  {0} ".format(histo[i][str(atomtype)] - histo[i-1][str(atomtype)],i*dx))
                else:
                    of.write("{0} ".format(histo[i][str(atomtype)]))
                    #print("{1}  {0} ".format(histo[i][str(atomtype)],i*dx))
            of.write("\n")
            of.write('END\n')
            of.write('X SetScale x 0,{1}, {0};\n'.format('partial_comp_'+znum2sym.z2sym(atomtype),npix*dx-dx))
        of.close()

            
    def local_composition(self, outfile):
        """ Variable radius sliding average Goes pixel by pixel and calculates the
        composition around that pixel (within some radius) and assigns the center
        pixel that composition. Use a 256 x 256 x 256 matrix. """
        radius = 3.6 * 2
        npix = 64
        #mat = np.zeros((npix,npix,npix),dtype=np.float)
        #mat = np.zeros((npix,npix,npix),dtype={'names':['col1', 'col2', 'col3'], 'formats':['f4','f4','f4']})
        #mat = np.zeros((npix,npix,npix),dtype={'names':['40', '13', '29'], 'formats':['f4','f4','f4']})
        #mat = np.zeros((npix,npix,npix),dtype={'names':['id','data'], 'formats':['f4','f4']})
        #names = ['id','data']
        #formats = ['i4',('f4','f4','f4')]
        #mat = np.zeros((npix,npix,npix),dtype=dict(names = names, formats=formats))
        #mat = np.zeros((npix,npix,npix),dtype={'40':('i4',0), '29':('f4',0), '13':('f4',0)})
        print("Creating matrix...")
        mat = [[[{} for i in range(npix)] for j in range(npix)] for k in range(npix)]
        print("Finished creating matrix.")
        #print(repr(mat))
        dx = self.lx/npix
        dy = self.ly/npix
        dz = self.lz/npix
        for ii,i in enumerate(drange(-npix/2*dx,npix/2*dx-dx,dx)):
            print("On ii = {0}".format(ii))
            for jj,j in enumerate(drange(-npix/2*dy,npix/2*dy-dy,dy)):
                for kk,k in enumerate(drange(-npix/2*dz,npix/2*dz-dz,dz)):
                    atoms = self.get_atoms_in_cutoff( (i,j,k), radius )
                    comp = {}
                    for atom in atoms:
                        comp[str(atom.z)] = comp.get(str(atom.z),0) + 1.0
                    for key in comp:
                        comp[key] /= len(atoms)
                    #print(comp)
                    #mat[ii][jj][kk] = copy.copy(comp)
                    mat[ii][jj][kk] = comp
        of = open(outfile,'w')
        of.write('IGOR\n')
        for atomtype in self.atomtypes:
            of.write('\nWAVES/N=({0},{1},{2})\t {3}\nBEGIN\n'.format(npix,npix,npix,'partial_comp_'+znum2sym.z2sym(atomtype)))
            for layer in mat:
                for column in layer:
                    for value in column:
                        try:
                            of.write("{0} ".format(value[str(atomtype)]))
                        except KeyError:
                            of.write("{0} ".format(0.0))
                of.write("\n")
            of.write('END\n')
            of.write('X SetScale/P x 0,1,"", {0}; SetScale/P y 0,1,"", {0}; SetScale/P z 0,1,"", {0}; SetScale d 0,0,"", {0}\n'.format('partial_comp_'+znum2sym.z2sym(atomtype)))
        of.close()
        return mat



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

if __name__ == "__main__":
    main()
