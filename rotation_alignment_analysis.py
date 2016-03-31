import sys, os, math
import warnings
import numpy as np
import pickle
from model import Model
from atom import Atom
#from voronoi_3d import calculate_atom
import voronoi_3d as VT
import matplotlib
if sys.platform == 'darwin': matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from tabulate import tabulate
from collections import defaultdict, Counter
HOME = os.environ['HOME']


def cart2sph(x,y,z):
    x2y2 = x**2 + y**2
    r = math.sqrt(x2y2 + z**2)               # r
    elev = math.atan2(z, math.sqrt(x2y2))    # theta
    az = math.atan2(y,x)                     # phi
    if az < -math.pi/2.: az += math.pi
    elif az > math.pi/2.: az -= math.pi
    if elev < -math.pi/2.: elev += math.pi
    elif elev > math.pi/2.: elev -= math.pi
    return r, elev, az



class Group(object):
    """ This object is designed to hold a group of AlignedData's and includes functions to perform operations on that group as a whole. """
    def __init__(self, clusters, comment=None):
        if not isinstance(clusters, list) or isinstance(clusters, np.ndarray):
            raise Exception("The clusters need to be contained in an ordered list so that the 'successful' mask can be applied correctly.")
        self.clusters = clusters
        self.comment = comment
        # Check to see if every atom has the same number of atoms
        counter = Counter(cluster.natoms for cluster in self.clusters if cluster.successful)
        self._natoms = counter.most_common(1)[0][0]
        if len(counter) != 1:
            warnings.warn("Not all clusters in this group have the same number of atoms. {0}\n  Ignoring those that comprise < 1% of the total.".format(counter), RuntimeWarning)
            total = sum(val for k,val in counter.items())
            for k,v in counter.items():
                if v/float(total) < 0.01: counter[k] = False
            for cluster in self.clusters:
                if cluster.successful and not counter[cluster.natoms]:
                    cluster.successful = False
        self._combined = None
        self._average = None

    def __getitem__(self, key):
        return self.clusters[key]

    def __len__(self):
        return self.nclusters

    @property
    def nclusters(self):
        """ The number of clusters in this Group. """
        return len(self.clusters)


    @property
    def successful(self):
        """ A mask to operate only on successfully aligned clusters. """
        return np.array([cluster.successful for cluster in self.clusters])
    

    @property
    def nsuccessful(self):
        """ Returns the number of successfully aligned clusters in this Group. """
        return sum(self.successful)


    @property
    def coordinates(self):
        """ A numpy array of shape (3, nclusters) with the coordination positions of each atom for each cluster. """
        return np.array( [ [a.coord[0], a.coord[1], a.coord[2]] for a in cluster.atoms for cluster in self.clusters] ).T


    @property
    def natoms(self):
        """ Returns the number of atoms in the clusters. If there are different numbers of atoms in the clusters, an exception is raised. """
        return self._natoms


    def average_structure(self, force_update=False):
        """ Calculates the atom positions of the average structure of the aligned clusters in the group. """
        if not force_update and self._average is not None:
            return self._average

        nsuccessful = sum(self.successful)
        avg_coords = [[0.0, 0.0, 0.0] for i in range(self.natoms)]
        for cluster in self.clusters:
            if not cluster.successful: continue
            # The order of atoms in aligned_targets is the same for every target so we don't need to reorder by result.ind. That has been done during the alignment.
            for i,atom in enumerate(cluster.aligned_target):
                if i == len(cluster.aligned_target): break # This is a weird bug...
                avg_coords[i][0] += atom.coord[0]/nsuccessful
                avg_coords[i][1] += atom.coord[1]/nsuccessful
                avg_coords[i][2] += atom.coord[2]/nsuccessful

        m = Cluster(center_included=False, comment='averaged structure', xsize=100.,ysize=100.,zsize=100., atoms=[Atom(i, 'Si', *coord) for i,coord in enumerate(avg_coords)])
        m.add(Atom(m.natoms, 'Si', *[0., 0., 0.])) # Add a center atom at (0,0,0)
        m.center = m.atoms[-1]
        m.center.neighs = m.atoms[:-1]
        #m.rescale_bond_distances(3.5)
        VT.calculate_atom(m, atom=m.center, cutoff=max(m.dist(m.atoms[-1],atom) for atom in m.atoms[:-1])+0.1) #VP calculation
        #print("  Center atom's VP index is {0}".format(m.center.vp))
        self._average = m
        return m


    def combine(self, colorize=True, force_update=False):
        """ Combines the group of aligned clusters into a single model. """
        if not force_update and self._combined is not None:
            return self._combined

        m = Model(comment='combined structures', xsize=100.,ysize=100.,zsize=100., atoms=[])
        count = 0
        atom_types = ['Si', 'Na', 'Mg', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'B', 'Al', 'Ga', 'C', 'Sn', 'Pb', 'O']
        for i,cluster in enumerate(self.clusters):
            if not cluster.successful: continue
            for j,atom in enumerate(cluster.aligned_target):
                if j == len(cluster.aligned_target.atoms): break
                if colorize:
                    new = Atom(count, atom_types[j], *atom.coord)
                else:
                    new = Atom(count, 'Si', *atom.coord)
                count += 1
                m.add(new)
        self._combined = m
        return m


    @property
    def site_occupancy(self):
        """ """
        c = Counter()
        for a in self.clusters:
            if not a.successful: continue
            c += Counter(a.order)
        return c


    def stdev(self, avg=None):
        if avg is None:
            avg = self.average_structure()
        for cluster in self.clusters:
            if not cluster.successful: continue
            for atom in cluster.aligned_target.atoms:
                closest = None
                mindist = float('inf')
                for center in avg.atoms:
                    dist = avg.dist(atom, center)
                    if dist < mindist:
                        mindist = dist
                        closest = center
                if mindist == float('inf'): raise Exception("Closest center atom not found!")
                if closest.neighs is None: closest.neighs = []
                closest.neighs.append(atom)
        stdevs = []
        for i,atom in enumerate(avg.atoms[:-1]):
            if atom.neighs is None: continue
            X = np.array([a.coord[0] for a in atom.neighs])
            Y = np.array([a.coord[1] for a in atom.neighs])
            Z = np.array([a.coord[2] for a in atom.neighs])
            #X = np.array([cart2sph(*a.coord)[0] for a in atom.neighs])
            #Y = np.array([cart2sph(*a.coord)[1] for a in atom.neighs])
            #Z = np.array([cart2sph(*a.coord)[2] for a in atom.neighs])
            std = (np.std(X) + np.std(Y) +np.std(Z))/3.
            stdevs.append(std)
        return np.mean(stdevs)#, np.std(stdevs)


    @property
    def mean_error(self):
        return sum(cluster.error for cluster in self.clusters if cluster.successful) / sum(self.successful)



class AlignedData(object):
    """ This object takes in the data from the alignment procedure and provides functions to operate on that data. """
    def __init__(self, data, reorder=False, fix_pbcs=True, recenter=True, comment=None):
        # Add a center atom to self.model, self.target, and self.aligned_target at (0,0,0)
        self.comment = comment
        self.scale_factor = None
        if data.ind is not None:
            self.successful = True
            self.model = Cluster(center_included=False, fix_pbcs=fix_pbcs, recenter=recenter, comment='data.model', xsize=100.,ysize=100.,zsize=100., atoms=[Atom(i, 'Si', *coord) for i,coord in enumerate(data.model)])
            self.target = Cluster(center_included=False, fix_pbcs=fix_pbcs, recenter=recenter, comment='data.target', xsize=100.,ysize=100.,zsize=100., atoms=[Atom(i, 'Si', *coord) for i,coord in enumerate(data.target)])
            self.order = data.ind[0]
            self.rotation = data.transformation_R
            self.translation = data.transformation_t.T[0]

            if not reorder:
                self.aligned_target = Cluster(center_included=False, fix_pbcs=fix_pbcs, recenter=recenter, comment='data.aligned_target', xsize=100.,ysize=100.,zsize=100., atoms=[Atom(i, 'Si', *coord) for i,coord in enumerate(data.aligned_target)])
            else:
                #print(self.order)
                while max(self.order) > len(self.order)-1:
                    for i, x in enumerate(self.order):
                        if i not in self.order:
                            break
                    for j, x in enumerate(self.order):
                        if x > i:
                            self.order[j] = x-1
                self.order = [i+1 for i in self.order]
                self.aligned_target = Cluster(center_included=False, fix_pbcs=fix_pbcs, recenter=recenter, comment='data.aligned_target', xsize=100.,ysize=100.,zsize=100., atoms=[Atom(i, 'Si', *coord) for i,coord in enumerate([data.aligned_target[j-1] for j in self.order])])

            self.error = data.error # L1-norm
        else:
            self.successful = False
            self.model = []
            self.target = []


    @property
    def natoms(self):
        return self.model.natoms # Or target? or ind? What happens when the target and model do not have the same number of atoms?



class Cluster(Model):
    """ An extension of Model that defines helper functions for small clusters. """
    def __init__(self, center_included=True, fix_pbcs=True, recenter=True, **kwargs):
        super(Cluster, self).__init__(**kwargs) # Send input args to Model constructor

        # If we don't run these two lines first, self.find_center_atom() may return the incorrect atom!
        if fix_pbcs:
            self.fix_cluster_pbcs()
        if recenter:
            self.recenter()

        if center_included:
            # Remove the center atom and put it on the end (if it isn't already there)
            center = self.find_center_atom()
            if center != self.atoms[-1]:
                self.remove(center)
                self.add(center)
            self.center_included = True
        else:
            self.center_included = False


    def fix_cluster_pbcs(self):
        self.recenter()

        # Move all the atoms on the - side to the + side if we need to fix the cluster
        for atom1 in self.atoms:
            for atom2 in self.atoms:
                x = (atom1.coord[0] - atom2.coord[0])
                y = (atom1.coord[1] - atom2.coord[1])
                z = (atom1.coord[2] - atom2.coord[2])
                if round(x/self.xsize) == -1:
                    atom1.coord = (atom1.coord[0] +self.xsize, atom1.coord[1], atom1.coord[2])
                if round(y/self.ysize) == -1:
                    atom1.coord = (atom1.coord[0], atom1.coord[1] +self.ysize, atom1.coord[2])
                if round(z/self.zsize) == -1:
                    atom1.coord = (atom1.coord[0], atom1.coord[1], atom1.coord[2] +self.zsize)
                dpbc = round(self.dist(atom1, atom2))
                npbc = round(self.dist(atom1, atom2, pbc=False))

        for atom1 in self.atoms:
            for atom2 in self.atoms:
                if dpbc != npbc:
                    raise Exception("PBC fix failed!", dpbc, npbc)

        self.recenter()
        return self


    def rescale_bond_distances(self, avg):
        """ Rescales a cluster so that the average bond length is 'avg' """
        center = self.find_center_atom()

        # Place the center atom at (0,0,0) and move every other atom relative to that translation
        for atom in self.atoms:
            atom.coord = (atom.coord[0]-center.coord[0], atom.coord[1]-center.coord[1], atom.coord[2]-center.coord[2])

        current_avg = 0.
        for atom in self.atoms:
            if atom != center:
                current_avg += self.dist(center, atom)
        current_avg /= (self.natoms-1)
        for atom in self.atoms:
            if atom.id != center.id:
                atom.coord = (atom.coord[0]/current_avg*avg, atom.coord[1]/current_avg*avg, atom.coord[2]/current_avg*avg)
        self.recenter()
        return avg/current_avg


    def normalize_bond_distances(self):
        return self.rescale_bond_distances(avg=1.0)


    def find_center_atom(self):
        """ Returns the center atom in self.atoms. """
        return sorted([(atom.coord[0]**2 + atom.coord[1]**2 + atom.coord[2]**2, atom) for atom in self.atoms])[0][1]


    def find_closest(self, atom, other_model=None):
        """ Returns the index of the closest atom to 'atom' in self.atoms or in other_model if specified. """
        if other_model is None:
            other_model = self
        return sorted([((atom.coord[0]-atom2.coord[0])**2 + (atom.coord[1]-atom2.coord[1])**2 + (atom.coord[2]-atom2.coord[2])**2, atom2) for atom2 in other_model.atoms])[0][1]
        

    def closest_list(self, atom):
        """ Returns a sorted list of (dist**2, atom) of the closest atoms to 'atom' in self.atoms. """
        return sorted([((atom.coord[0]-atom2.coord[0])**2 + (atom.coord[1]-atom2.coord[1])**2 + (atom.coord[2]-atom2.coord[2])**2, atom2) for atom2 in self.atoms])


    def L1Norm(self, other):
        natoms = min(self.natoms, other.natoms) - self.center_included
        L1 = 0.0
        for i in range(natoms):
            L1 += abs(self.atoms[i].coord[0] - self.find_closest(self.atoms[i], other_model=other).coord[0])
            L1 += abs(self.atoms[i].coord[1] - self.find_closest(self.atoms[i], other_model=other).coord[1])
            L1 += abs(self.atoms[i].coord[2] - self.find_closest(self.atoms[i], other_model=other).coord[2])
        return L1/natoms


    def L2Norm(self, other):
        natoms = min(self.natoms, other.natoms) - self.center_included
        L2 = 0.0
        for i in range(natoms):
            L2 += ((self.atoms[i].coord[0] - other.atoms[i].coord[0])**2 + 
                   (self.atoms[i].coord[1] - other.atoms[i].coord[1])**2 + 
                   (self.atoms[i].coord[2] - other.atoms[i].coord[2])**2)
        return math.sqrt(L2)/natoms


    def L2Norm2(self, other):
        natoms = min(self.natoms, other.natoms) - self.center_included
        L2 = 0.0
        for i in range(natoms):
            L2 += ((self.atoms[i].coord[0] - other.atoms[i].coord[0])**2 + 
                  (self.atoms[i].coord[1] - other.atoms[i].coord[1])**2 + 
                  (self.atoms[i].coord[2] - other.atoms[i].coord[2])**2)
        return L2/natoms


    def LinfNorm(self, other):
        natoms = min(self.natoms, other.natoms) - self.center_included
        Linf = [0.0 for i in range(natoms)]
        for i in range(natoms):
            Linf[i] = (abs(self.atoms[i].coord[0] - other.atoms[i].coord[0]) + 
                       abs(self.atoms[i].coord[1] - other.atoms[i].coord[1]) + 
                       abs(self.atoms[i].coord[2] - other.atoms[i].coord[2]))
        return max(Linf)/natoms



class Aligned_data:
    ind = None
    model = None
    target = None
    transformation_R = None
    transformation_t = None
    mapped_target = None
    error = None
    def any(self):
        if self.ind == None:
            return False
        else:
            return True



def log_normal(x, y0, x0, amp, width):
    return y0 + amp * np.exp( -(np.log(x/x0)/width)**2 )


def log_normal_large(x):
    return log_normal(x, y0=0, amp=3183.3, x0=0.061456, width=0.26866)
    #return log_normal(x, y0=-4.8386, amp=3183.3, x0=0.061456, width=0.26866)


def log_normal_small(x):
    return log_normal(x, y0=0, amp=373.16, x0=0.037255, width=0.17297)
    #return log_normal(x, y0=-4.8386, amp=373.16, x0=0.037255, width=0.17297)


def load_pkl(filename):
    with warnings.catch_warnings(record=True) as w: # Suppress the warnings...
        o = open(filename, 'rb')
        data = pickle.load(o)
        o.close()
    return data


def load_vp(filename, start=0):
    content = open(filename, 'r').readlines()
    for i,line in enumerate(reversed(content)):
        if 'CN:' in line:
            end = len(content)-i
            break
    vp_lines = content[start:end]
    vp_lines = [line.split(':')[0] for line in vp_lines]
    #vp = reversed([eval(tup) for tup in vp_lines])
    vp = [eval(tup) for tup in vp_lines]
    return vp


def get_bond_scaling(dir, i):
    f = '/home/jjmaldonis/working/ZrCuAl/md_80k/{dir}/{dir}.{i}.xyz'.format(dir=dir, i=i)
    m = Model(f)
    index = m.comment.index('bond length scaling factor is')
    index += len('bond length scaling factor is ')
    x = float(m.comment[index:])
    return x




def main():
    #pkl_files = ['00120_aligned_to_average_polys/00120_aligned_to_avg0012_8_as_model.pkl', '00120_aligned_to_average_polys/00120_aligned_to_avg0012_as_model.pkl', '00120_aligned_to_average_polys/00120_aligned_to_avg0282_9_as_model.pkl', '00120_aligned_to_average_polys/00120_aligned_to_avg0282_as_model.pkl',
    #             '0282_aligned_to_average_polys/0282_aligned_to_avg0012_8_as_model.pkl', '0282_aligned_to_average_polys/0282_aligned_to_avg0012_as_model.pkl', '0282_aligned_to_average_polys/0282_aligned_to_avg0282_9_as_model.pkl', '0282_aligned_to_average_polys/0282_aligned_to_avg0282_as_model.pkl']
    #pkl_files = [os.path.join(root, name) if '.pkl' in name for name in files for root, dirs, files in os.walk('.')]

    #pkl_file = './pkls/04410000/04410000-icos.pkl'
    #pkl_file = './pkls/04531000/04531000-icos.pkl'
    #pkl_file = './pkls/02820000/02820000-icos.pkl'
    #pkl_file = './pkls/04440000/04440000-icos.pkl'
    #pkl_file = './pkls/05260000/05260000-icos.pkl'
    #aligned = load_pkl(pkl_file)
    #pkl_files = [os.path.join('./pkls/', f, os.listdir('./pkls/'+f)[0]) for f in os.listdir('./pkls/')]
    #pkl_files = ['./pkls/001200000/001200000-icos.pkl']
    #pkl_files = ['./pkls/04440000/04440000-icos.pkl']

    VP = load_vp(HOME+'/working/Arash_uploads/vp.txt', start=2)
    #pkl_files = ['./pkls/' + ''.join([str(x) for x in vp]) + '/' + ''.join([str(x) for x in vp]) + '-icos.pkl' for vp in VP]
    #pkl_files = ['./pkls/' + ''.join([str(x) for x in vp]) + '/' + ''.join([str(x) for x in vp]) + '-less_than_12.pkl' for vp in VP]
    pkl_files =[]
    for vp in VP:
        f = HOME+'/working/Arash_uploads/pkls/' + ''.join([str(x) for x in vp]) + '/' + ''.join([str(x) for x in vp]) + '-icos.pkl'
        pkl_files.append(f)
        f = HOME+'/working/Arash_uploads/pkls/' + ''.join([str(x) for x in vp]) + '/' + ''.join([str(x) for x in vp]) + '-less_than_12.pkl'
        pkl_files.append(f)
    temp = []
    for f in pkl_files:
        try:
            open(f)
            temp.append(f)
        except:
            pass
    pkl_files = temp

    #pkl_files = [os.path.join(root, name) for root, dirs, files in os.walk('./pkls/') for name in files if '.pkl' in name]

    #total = 0
    #d = defaultdict(int)
    #for f,pkl_file in enumerate(pkl_files):
    #    #if sum(VP[f]) >= 12: continue
    #    aligned = load_pkl(pkl_file)
    #    print("Finished loading data pkl file {0}.".format(pkl_file))
    #    #g = Group([AlignedData(a) for a in aligned], comment = ''.join([str(x) for x in VP[f]]))
    #    g = Group([AlignedData(a) for a in aligned])
    #    for cluster in g.clusters:
    #        if cluster.successful and cluster.error > 0.15:
    #            #d[cluster.natoms] += 1
    #            d[pkl_file] += 1
    #    #g.combine().write(os.path.basename(pkl_file)+'.xyz')
    #    print(d)
    #    total += g.nclusters
    #    print(total)
    #return 0

    #pkl_files = []
    #pkl_files += ['/home/jjmaldonis/working/Arash_uploads/pkls/02800000/02800000-less_than_12.pkl']
    #pkl_files += ['/home/jjmaldonis/working/Arash_uploads/pkls/001200000/001200000-icos.pkl']
    #pkl_files += ['/home/jjmaldonis/working/Arash_uploads/pkls/04440000/04440000-icos.pkl']

    print(pkl_files)

    table = {"Name":[], "Mean Error":[], "VP":[], "Number of clusters":[], "Stdev":[], "Large":[], "Small":[]}
    groups = []
    for f,pkl_file in enumerate(pkl_files):
        head, tail = os.path.split(pkl_file)

        aligned = load_pkl(pkl_file)
        print("Finished loading data pkl file {0}.".format(pkl_file))
        #g = Group([AlignedData(a, reorder=False) for a in aligned], comment = ''.join([str(x) for x in VP[f]]))
        g = Group([AlignedData(a, reorder=False) for a in aligned])
        groups.append(g)
        print("Finished rendering into Group!")

        print("Rescaling every cluster...")
        vp_string = os.path.split(head)[-1]
        i = 0
        for cluster in g.clusters:
            if not cluster.successful: continue
            scaling = get_bond_scaling(vp_string, i)
            cluster.model.rescale_bond_distances(scaling)
            cluster.target.rescale_bond_distances(scaling)
            cluster.aligned_target.rescale_bond_distances(scaling)
            l2 =  cluster.aligned_target.L2Norm(cluster.model)
            #print(l2 / cluster.error)
            cluster.error = l2
            i += 1
        print("Finished rescaling every cluster!")


        print("  {0}".format(VP[f]))
        table["Name"].append('<{0}>'.format(','.join(str(x) for x in VP[f])))
        print("  Number of clusters: {0}".format(g.nclusters))
        table["Number of clusters"].append(g.nclusters)
        print("  Number of atoms in those clusters: {0}".format(g.natoms))

        avg = g.average_structure()
        avg.write(os.path.join(head, tail[:-4] + '-average_structure.xyz'))
        print("  Saved average structure to {0}".format(os.path.join(head, tail[:-4] + '.xyz')))
        combined = g.combine()
        combined.write(os.path.join(head, tail[:-4] + '-combined.xyz'))
        try:
            print("  VP of average structure: {0}".format(avg.center.vp))
            table["VP"].append(avg.center.vp)
        except Exception as error:
            print("  VP calculation of avg structure failed.")
            table["VP"].append(None)

        print("  Mean error: {0}".format(g.mean_error))
        table["Mean Error"].append(g.mean_error)
        std = g.stdev(avg=avg)
        print("  Stdev of cluster positions: {0}".format(std))
        table["Stdev"].append(std)

        amp_large = 3183.3
        amp_small = 373.16
        #large = np.zeros(np.sum(g.successful), dtype=np.float)
        #small = np.zeros(np.sum(g.successful), dtype=np.float)
        i = 0
        of = open(os.path.join(head, 'output.txt'), 'w')
        #of.write('cluster_number, error, P_small, Linf, L1\n')
        of.write('cluster_number, error, Linf, L1\n')
        for cluster in g.clusters:
            if not cluster.successful: continue
            #g_large = log_normal_large(cluster.error)#/amp_large
            #g_small = log_normal_small(cluster.error)#/amp_small
            #large[i] = g_large / (g_large+g_small)
            #small[i] = g_small / (g_large+g_small)
            #print(i, cluster.error, small[i], large[i], cluster.aligned_target.LinfNorm(cluster.model), cluster.aligned_target.L1Norm(cluster.model))
            #of.write('{0}, {1}, {2}, {3}, {4}\n'.format(i, cluster.error, small[i], cluster.aligned_target.LinfNorm(cluster.model), cluster.aligned_target.L1Norm(cluster.model)))
            of.write('{0}, {1}, {2}, {3}\n'.format(i, cluster.error, cluster.aligned_target.LinfNorm(cluster.model), cluster.aligned_target.L1Norm(cluster.model)))
            i += 1
        of.close()
        #print("  Probability that one of these clusters might be in the large peak: {0}".format(np.mean(large)))
        #table["Large"].append(np.mean(large))
        #print("  Probability that one of these clusters might be in the small peak: {0}".format(np.mean(small)))
        #table["Small"].append(np.mean(small))

        print('')
        print(tabulate(table, headers="keys"))
        print('')
        open('table.txt', 'w').write(tabulate(table, headers="keys"))
    return 0




    # THIS IS OLD CODE but I haven't decided whether I want to make it or part of it into a helper function or not
    pkl_files = [os.path.join(root, name) for root, dirs, files in os.walk('.') for name in files if '.pkl' in name]
    aligned_files = sorted([(pkl, load_pkl(pkl)) for pkl in pkl_files], key=lambda tup: len(tup[1]))
    print("Finished unpacking pickles.")

    avg_errors = []
    other = []
    filecount = 0
    for filename, aligned in aligned_files:
        output_filename = filename[:-4] + '.xyz'
        print("Analyzing {1}: {0}...".format(filename, filecount))
        fix_results(aligned)

        # Basic error analysis
        nbad = sum(1 for a in aligned if a.error is None)
        errors = [a.error/len(a.ind) for a in aligned if a.error is not None]
        avg_error = np.mean(errors, axis=0)
        stdev_error = np.std(errors, axis=0)
        print("  Error: {0} +- {1}  ({2} %)".format(avg_error, stdev_error, 100.0*stdev_error/avg_error))
        print("  Number of unalignable models was {0}  ({1} %)".format(nbad, float(nbad)/len(aligned)*100))
        print("  Total number of models was {0}  ({1} %)".format(len(aligned), 100.0*len(aligned)/91200))
        print('')
        if aligned[0].ind is not None:
            if len(aligned[0].ind) > 11:
                avg_errors.append((avg_error, stdev_error))
            else:
                other.append((avg_error, stdev_error))
        filecount += 1

    out = open("output.txt", 'w')
    out.write("error, stdev\n")
    for err, std in avg_errors:
        out.write("{0}, {1}\n".format(err, std))
    for err, std in other:
        out.write("{0}, {1}\n".format(err, std))
    out.close()
    return 0

    output_filename = 'test'

    m = combine_aligned_into_single_model(aligned)
    path, output_filename = os.path.split(output_filename)
    m.write(os.path.join(path, 'all_'+output_filename))

    m = average_structure(aligned)
    m.write(os.path.join(path, 'averaged_'+output_filename))

    return 0



if __name__ == '__main__':
    main()


