#import matplotlib
#matplotlib.use('PDF')
import sys, os
import copy
#import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Voronoi, voronoi_plot_2d
from pprint import pprint
from atom import Atom
from hutch import Hutch
import znum2sym
import math
from collections import defaultdict, Counter
from tools import drange


""" Add these functions into Model:
"""

class Masses(object):
    def __init__(self):
        self.masses = {1:1.007947, 2:4.0026022, 3:6.9412, 4:9.0121823, 5:10.8117, 6:12.01078, 7:14.00672, 8:15.99943, 9:18.99840325, 10:20.17976, 11:22.989769282, 12:24.30506, 13:26.98153868, 14:28.08553, 15:30.9737622, 16:32.0655, 17:35.4532, 18:39.9481, 19:39.09831, 20:40.0784, 21:44.9559126, 22:47.8671, 23:50.94151, 24:51.99616, 25:54.9380455, 26:55.8452, 27:58.9331955, 28:58.69342, 29:63.5463, 30:65.4094, 31:69.7231, 32:72.641, 33:74.921602, 34:78.963, 35:79.9041, 36:83.7982, 37:85.46783, 38:87.621, 39:88.905852, 40:91.2242, 41:92.906382, 42:95.942, 43:98, 44:101.072, 45:102.905502, 46:106.421, 47:107.86822, 48:112.4118, 49:114.8183, 50:118.7107, 51:121.7601, 52:127.603, 53:126.904473, 54:131.2936, 55:132.90545192, 56:137.3277, 57:138.905477, 58:140.1161, 59:140.907652, 60:144.2423, 61:145, 62:150.362, 63:151.9641, 64:157.253, 65:158.925352, 66:162.5001, 67:164.930322, 68:167.2593, 69:168.934212, 70:173.043, 71:174.9671, 72:178.492, 73:180.947882, 74:183.841, 75:186.2071, 76:190.233, 77:192.2173, 78:195.0849, 79:196.9665694, 80:200.592, 81:204.38332, 82:207.21, 83:208.980401, 84:210, 85:210, 86:220, 87:223, 88:226, 89:227, 91:231.035882, 90:232.038062, 93:237, 92:238.028913, 95:243, 94:244, 96:247, 97:247, 98:251, 99:252, 100:257, 101:258, 102:259, 103:262, 104:261, 105:262, 106:266, 107:264, 108:277, 109:268, 110:271, 111:272, 112:285, 113:284, 114:289, 115:288, 116:292, 118:293}

    def get_mass(self, znum):
        return self.masses[znum]

    def get_znum(self, mass):
        prec1 = len(str(mass)[int(math.log10(mass))+2:])
        for z, m in self.masses.items():
            prec2 = len(str(m)[int(math.log10(m))+2:])
            p = min(prec1, prec2)
            if round(mass, p) == round(m, p): return z
        raise Exception("Mass not found!")


masses = Masses()


class Model(object):
    """
        Holds an atomic model and defines a set of helper functions for it.
        functions:
            write(outfile=stdout):  Supported extensions are xyz, dat, and a simple version of cif
            generate_neighbors(cutoff):  Generates neighbors for every atom using the supplied cutoff (which can be a float or a dictionary)
            get_atoms_in_cutoff(atom,cutoff)
            nearest_neigh(atom)
    """
    
    def __init__(self, modelfilename=None, comment=None, xsize=None, ysize=None, zsize=None, atoms=None):
        """ sets:
                self.comment
                self.xsize
                self.ysize
                self.zsize
                self.atoms
                self.natoms
                self.atomtypes
                self.xx, self.yy, self.zz
                self.filename
            Extra routines can set:
                self.coord_numbers """

        self.filename = modelfilename
        if self.filename is not None:
            self._load()
        else:
            self.comment = comment
            self.xsize = xsize
            self.ysize = ysize
            self.zsize = zsize
            self.atoms = [atom.copy() for atom in atoms]
            self.natoms = len(self.atoms)

        if(self.xsize and self.ysize and self.zsize):
            self.hutch = Hutch(self)

        if(not hasattr(self,'atomtypes')):
            self.atomtypes = Counter(atom.z for atom in self.atoms)

    def __contains__(self, key):
        return key in self.atoms

    def __getitem__(self, atomid):
        try:
            if(self.atoms[atomid].id == atomid):
                return self.atoms[atomid]
            else:
                raise Exception("Atom positions have been shuffled or you gave an out of bounds index.")
        except:
            for atom in self.atoms:
                if(atom.id == atomid):
                    return atom
        return None

    def __len__(self):
        assert self.natoms == len(self.atoms)
        return self.natoms
        
    def add(self, atom, reset_id=False):
        """ Adds atom 'atom' to the model. Note that there is no error checking to prevent
        the user from adding an atom twice. Be careful not to do that. """
        if reset_id:
            atom.id = self.natoms
        self.atoms.append(atom)
        self.natoms += 1
        self.atomtypes[atom.z] += 1
        try:
            self.hutch.add_atom(atom)
        except AttributeError:
            pass

    def remove(self, atom):
        """ Removes atom 'atom' from the model """
        try:
            self.hutch.remove_atom(atom)
        except:# AttributeError or ValueError:
            pass
        self.atoms.remove(atom)
        self.natoms -= 1
        self.atomtypes[atom.z] -= 1

    def _load(self):
        filename, ext = os.path.splitext(self.filename)
        if(ext == 'xyz' or ext == '.xyz'):
            self._load_xyz()
        elif(ext == 'dat' or ext == '.dat'):
            self._load_dat()
        else:
            raise Exception("Unknown input model type! File {0} has extension {1}".format(self.filename,ext))

    def _load_xyz(self):
        modelfile = self.filename
        self.atoms = []
        self.atomtypes = defaultdict(int)
        self.natoms = 0
        with open(modelfile) as f:
            natoms = int(f.readline().strip())
            self.comment = f.readline().strip()
            try:
                self.xsize,self.ysize,self.zsize = tuple([float(x) for x in self.comment.split()[:3]])
            except:
                self.xsize,self.ysize,self.zsize = (None,None,None)
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
                self.xsize = abs(float(line.split()[0])) + abs(float(line.split()[1]))
            if('ylo' in line and 'yhi' in line):
                self.ysize = abs(float(line.split()[0])) + abs(float(line.split()[1]))
            if('zlo' in line and 'zhi' in line):
                self.zsize = abs(float(line.split()[0])) + abs(float(line.split()[1]))
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

    def write(self, outfile=None, ext=None):
        if(outfile is not None and ext is None): _,ext = os.path.splitext(outfile)
        elif(ext is None): ext = '.xyz'
        if(ext == '.xyz' or ext == 'xyz'):
            self._write_xyz(outfile)
        elif(ext == '.dat' or ext == 'dat'):
            self._write_dat(outfile)
        elif(ext == '.cif' or ext == 'cif'):
            self._write_cif(outfile)
        return ''

    def _write_dat(self, outfile=None):
        if outfile is None: of = sys.stdout
        else:               of = open(outfile,'w')
        of.write(self.comment+'\n')
        of.write('{0} atoms\n\n'.format(self.natoms))
        of.write('{0} atom types\n\n'.format(len(self.atomtypes)))
        of.write('{0} {1} xlo xhi\n'.format(-self.xsize/2,self.xsize/2))
        of.write('{0} {1} ylo yhi\n'.format(-self.ysize/2,self.ysize/2))
        of.write('{0} {1} zlo zhi\n\n'.format(-self.zsize/2,self.zsize/2))
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

    def _write_xyz(self, outfile=None):
        if outfile is None: of = sys.stdout
        else:               of = open(outfile,'w')
        of.write(str(self.natoms)+'\n')
        of.write("{1} {2} {3} {0}\n".format(self.comment.strip(),self.xsize, self.xsize, self.zsize))
        for i,atom in enumerate(self.atoms):
            of.write(atom.realxyz()+'\n')

    def _write_cif(self, outfile=None):
        if outfile is None: of = sys.stdout
        else:               of = open(outfile,'w')
        of.write('_pd_phase_name\t'+self.comment+'\n')
        of.write('_cell_length_a '+str(self.xsize)+'\n')
        of.write('_cell_length_b '+str(self.ysize)+'\n')
        of.write('_cell_length_c '+str(self.zsize)+'\n')
        of.write('_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n')
        of.write('_symmetry_space_group_name_H-M         \'P 1\'\n_symmetry_Int_Tables_number            1\n\n')
        of.write('loop_\n_symmetry_equiv_pos_as_xyz\n   \'x, y, z\'\n\n')
        of.write('loop_\n   _atom_site_label\n   _atom_site_occupancy\n   _atom_site_fract_x\n   _atom_site_fract_y\n   _atom_site_fract_z\n   _atom_site_adp_type\n   _atom_site_B_iso_or_equiv\n   _atom_site_type_symbol\n')
        atomtypes = {}
        for atom in self.atoms:
            atomtypes[atom.z] = atomtypes.get(atom.z, 0) + 1
            of.write('   '+znum2sym.z2sym(atom.z)+str(atomtypes[atom.z])+'\t1.0\t'+str(atom.coord[0]/self.xsize+0.5)+'\t'+str(atom.coord[1]/self.ysize+0.5)+'\t'+str(atom.coord[2]/self.zsize+0.5)+'\tBiso\t1.000000\t'+znum2sym.z2sym(atom.z)+'\n')
        of.write('\n')
        of.close()

    def __str__(self):
        return self.write()

    @property
    def composition(self):
        d = {}
        for key in self.atomtypes:
            d[znum2sym.z2sym(key)] = self.atomtypes[key]/float(self.natoms)
        return d

    @property
    def coordinates(self):
        xx = np.fromfunction(lambda i: self.atoms[i].coord[0], self.natoms, dtype=np.float)
        yy = np.fromfunction(lambda i: self.atoms[i].coord[1], self.natoms, dtype=np.float)
        zz = np.fromfunction(lambda i: self.atoms[i].coord[2], self.natoms, dtype=np.float)
        return xx,yy,zz

    def generate_neighbors(self, cutoff):
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
        """ atom.neighs must be defined first for all atoms
            Form will be:
                {'Cu-Al': 4.5, ... }
        """
        coord_numbers = {}
        for typea in self.atomtypes:
            coord_numbers[znum2sym.z2sym(typea)] = 0
            for typeb in self.atomtypes:
                coord_numbers[znum2sym.z2sym(typea)+'-'+znum2sym.z2sym(typeb)] = 0
        for atom in self.atoms:
            for n in atom.neighs:
                coord_numbers[znum2sym.z2sym(atom.z)] += 1
                coord_numbers[znum2sym.z2sym(atom.z)+'-'+znum2sym.z2sym(n.z)] += 1
        for key in coord_numbers:
            elem = znum2sym.sym2z(key.split('-')[0])
            coord_numbers[key] /= float(self.atomtypes[elem])
        return coord_numbers

    def get_atoms_in_cutoff(self,atom,cutoff):
        return self.hutch.get_atoms_in_radius(atom,cutoff)

    def dist(self, atom1, atom2, pbc=True):
        x = (atom1.coord[0] - atom2.coord[0])
        y = (atom1.coord[1] - atom2.coord[1])
        z = (atom1.coord[2] - atom2.coord[2])
        if pbc:
            x = x - self.xsize*round(x/self.xsize)
            y = y - self.ysize*round(y/self.ysize)
            z = z - self.zsize*round(z/self.zsize)
        return math.sqrt(x**2+y**2+z**2)

    def dist2(self, atom1, atom2, pbc=True):
        x = (atom1.coord[0] - atom2.coord[0])
        y = (atom1.coord[1] - atom2.coord[1])
        z = (atom1.coord[2] - atom2.coord[2])
        if pbc:
            x = x - self.xsize*round(x/self.xsize)
            y = y - self.ysize*round(y/self.ysize)
            z = z - self.zsize*round(z/self.zsize)
        return x**2+y**2+z**2

    def get_all_dists(self):
        dists = [0. for i in range(self.natoms*self.natoms)]
        count = 0
        for atomi in self.atoms:
            for atomj in self.atoms[self.atoms.index(atomi)+1:]:
                dists[count] = (atomi, atomj, self.dist(atomi,atomj))
                count += 1
        return dists

    def nearest_neigh_of_same_type(self, atom, cutoff=3.5):
        """ returns the nearest atom of the same species as 'atom' """
        atoms = []
        while(len(atoms) == 0):
            atoms = self.get_atoms_in_cutoff(atom, cutoff)
            atoms = [x for x in atoms if x.z == atom.z]
            #if atom in atoms: atoms.remove(atom)
            cutoff *= 2
        cutoff /= 2 # set back to the value used in case I want it later
        d = float("inf")
        for atomi in atoms:
            dt = self.dist(atom, atomi)
            if dt < d:
                d = dt
                a = atomi
        if(a.z != atom.z): raise Exception("Error! Function 'nearest_neigh_of_same_type' didn't work!")
        return a

    def nearest_neigh(self, atom):
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

    def print_bond_stats(self):
        # TODO Rewrite this function if I ever need it again. There was a "self.bonds = composition" line before the last for loop in "generate_average_coord_numbers" to define self.bonds, but I deleted that line.
        nbonds = 0.0
        for key in self.bonds:
            if '-' not in key:
                nbonds += self.bonds[key]
        bond_stats = self.bonds.copy()
        for key in bond_stats.keys():
            if '-' in key:
                del bond_stats[key]
            else:
                bond_stats[key] *= 1.0/(nbonds*self.composition[key])
        print('Bond statistics:')
        pprint(self.bonds)
        print('Ratio of actual/expected bonds:')
        pprint(bond_stats)

    def radial_composition(self, outfile):
        """ Creates 1D waves stored in outfile for each element in the model. Each 
        wave is a histogram of the number of atoms between two radial positions
        starting at the center of model and radiating outward. """
        # TODO Rewrite if I ever need this again
        npix = 16
        keys = self.atomtypes.keys()
        #histo = [[0.0 for x in range(npix)] for key in keys]
        histo = [{} for x in range(npix)]
        dx = (self.xsize/2.0)/npix # Cube assumed
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
        # TODO Rewrite if I ever need this again
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
        dx = self.xsize/npix
        dy = self.ysize/npix
        dz = self.zsize/npix
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

    def compare(self, m2):
        """ Compares self and m2. An exception is raised if every atom in self is not found in m2 (but not vice versa). """
        assert self.natoms == m2.natoms
        for atom1 in self.atoms:
            found = False
            for atom2 in m2.atoms:
                if(atom1 == atom2):
                    found = True
                    break
            else:
                raise Exception("Atom not found! {0}".format(atom1))
        else:
            pass

    def recenter(self, debug=False):
        xmin = min(atom.coord[0] for atom in self.atoms)
        xmax = max(atom.coord[0] for atom in self.atoms)
        ymin = min(atom.coord[1] for atom in self.atoms)
        ymax = max(atom.coord[1] for atom in self.atoms)
        zmin = min(atom.coord[2] for atom in self.atoms)
        zmax = max(atom.coord[2] for atom in self.atoms)
        if(debug):
            print("Original x-min/max = ({0},{1}".format(xmin,xmax))
            print("Original y-min/max = ({0},{1}".format(ymin,ymax))
            print("Original z-min/max = ({0},{1}".format(zmin,zmax))
        xcenter = xmin + (xmax - xmin)/2.0
        ycenter = ymin + (ymax - ymin)/2.0
        zcenter = zmin + (zmax - zmin)/2.0
        if(debug):
            print("Center was at ({0},{1},{2})".format(xcenter,ycenter,zcenter))
        for i in range(self.natoms):
            self.atoms[i].coord = ( self.atoms[i].coord[0] - xcenter, self.atoms[i].coord[1] - ycenter, self.atoms[i].coord[2] - zcenter )
        if(debug):
            xmin = min(atom.coord[0] for atom in self.atoms)
            xmax = max(atom.coord[0] for atom in self.atoms)
            print("Original x-min/max = ({0},{1}".format(xmin,xmax))
            ymin = min(atom.coord[1] for atom in self.atoms)
            ymax = max(atom.coord[1] for atom in self.atoms)
            print("Original y-min/max = ({0},{1}".format(ymin,ymax))
            zmin = min(atom.coord[2] for atom in self.atoms)
            zmax = max(atom.coord[2] for atom in self.atoms)
            print("Original z-min/max = ({0},{1}".format(zmin,zmax))
            xcenter = round(xmin + (xmax - xmin)/2.0,15)
            ycenter = round(ymin + (ymax - ymin)/2.0,15)
            zcenter = round(zmin + (zmax - zmin)/2.0,15)
            print("Center is at ({0},{1},{2})".format(xcenter,ycenter,zcenter))

    def crop(self, xstart, xend, ystart, yend, zstart, zend):
        atoms = []
        for atom in self.atoms:
            if( xstart < atom.coord[0] < xend and
                ystart < atom.coord[1] < yend and
                zstart < atom.coord[2] < zend):
                atoms.append(atom)
        newm = Model('cropped model',2*abs(xend)+2*abs(xstart),2*abs(yend)+2*abs(ystart),2*abs(zend)+2*abs(zstart),atoms)
        newm.recenter()
        return newm

    def min_dist(self):
        # TODO Can use the function "find_nearest_neighbor" and that should make this much faster
        m = float("inf")
        for i,atomi in enumerate(self.atoms):
            for atomj in self.atoms[self.atoms.index(atomi)+1:]:
                d = self.dist(atomi, atomj)
                if(d < m):
                    m = d
                    t = (atomi,atomj)
        print("\nMinimimum atomic spacing: {0} from atoms {1}\n".format(m,t))
        return m, t

    def atoms_in_box(self):
        for atom in self.atoms:
            if( atom.coord[0] < -self.xsize/2 or
                atom.coord[0] >  self.xsize/2 or
                atom.coord[1] < -self.ysize/2 or
                atom.coord[1] >  self.ysize/2 or
                atom.coord[2] < -self.zsize/2 or
                atom.coord[2] >  self.zsize/2):
                print('ERROR! Atom {0} is out of your box! Coord: {1} Box: {2}'.format(atom, atom.coord, (self.xsize/2, self.ysize/2, self.zsize/2)))
        print("\nYour atom coords are all inside the box.\n")

    def atom_density(self):
        ad = self.natoms / (self.xsize * self.ysize * self.zsize)
        print("\nAtom density: {0}".format(ad))
        return ad

    def min_max_positions(self):
        print('\nBox sizes: {0} x {1} x {2}'.format(self.xsize, self.ysize, self.zsize))
        print('Half box sizes: {0}, {1}, {2}'.format(self.xsize/2, self.ysize/2, self.zsize/2))
        xx = [atom.coord[0] for atom in self.atoms]
        yy = [atom.coord[1] for atom in self.atoms]
        zz = [atom.coord[2] for atom in self.atoms]
        print('x-min: {0}\t x-max: {1}'.format(min(xx),max(xx)))

    def combine(self, *models):
        m0 = self.copy()
        for m in models:
            for atom in m.atoms:
                if atom not in m0.atoms:
                    m0.add(atom)
        return m0

    def generate_larger_model(self, mult):
        # This is slow because I am adding atoms one at a time
        model = Model(xsize=self.xsize*mult, ysize=self.ysize*mult, zsize=self.zsize*mult, atoms=[])
        natoms = 0
        for atom in self.atoms:
            for i in range(mult):
                for j in range(mult):
                    for k in range(mult):
                        if(i == j == k == 0): continue
                        natoms += 1
                        x = atom.coord[0] + i*model.xsize
                        y = atom.coord[1] + j*model.ysize
                        z = atom.coord[2] + k*model.zsize
                        model.add_atom(Atom(natoms, atom.z, x, y, z))

        # Shift right 1/2 world and left mult/2 worlds
        # which is equivalent to left (mult-1)/2 worlds
        for i,atom in enumerate(model.atoms):
            model.atoms[i].set_coord(model.atoms[i].coord[0] - (mult-1)/2.0*model.xsize, model.atoms[i].coord[1] - (mult-1)/2.0*model.ysize, model.atoms[i].coord[2] - (mult-1)/2.0*model.zsize)
        return model

    def icofrac(self):
        """ Sets the atom.z for all atoms in self to be a number between 0 and nbins based on the fraction of pentagonal VP faces. """
        nbins = 6
        del_bin = 100.0/nbins
        fracs = []
        for atom in m.atoms:
            fracs.append((float(atom.vp.index[2])/float(sum(atom.vp.index))*100.0))
            bin = int( (float(atom.vp.index[2])/float(sum(atom.vp.index))*100.0) /(100.0/(nbins-1)))
            atom.z = bin+1
        fracs.sort()
        print('Min %: {0}. Max %: {1}'.format(min(fracs),max(fracs)))

    def bond_angle_distribution(m, nbins=None, dtheta=None):
        if nbins is None:
            nbins = np.pi/dtheta
        elif dtheta is None:
            dtheta = np.pi/nbins
        else:
            raise Exception("Either nbins or dtheta must be specified.")
        divisor = 0
        hist = np.zeros(nbins+1, np.int)
        for atomi in m.atoms:
            for j in range(0,atomi.cn):
                rij = m.dist(atomi,atomi.neighs[j])
                for k in range(j+1,atomi.cn):
                    rik = m.dist(atomi,atomi.neighs[k])
                    rjk = m.dist(atomi.neighs[j],atomi.neighs[k])
                    temp = round( (rij**2 + rik**2 - rjk**2)/(2.0*rij*rik) , 10 )
                    tijk = acos(temp)
                    bin = int(tijk/dtheta)
                    hist[bin] += 1
            divisor += atomi.cn*(atomi.cn-1)
        #if(divisor > 0):
        #    hist = [float(x)/float(divisor)*200 for x in hist]

        dtheta = np.arange(0.0,180.0+dtheta*180.0/np.pi,dtheta*180.0/np.pi)
        return dtheta,hist



def main():
    m = Model(sys.argv[1])
    outflag = False
    if(len(sys.argv) > 3):
        if(sys.argv[2] == '-o'):
            outtype = sys.argv[3]
            m.write(ext=outtype)

if __name__ == "__main__":
    main()
