import os
import math
import numpy as np
import scipy.spatial.distance
from collections import Counter
import warnings

from lazy_property import lazyproperty


def load_alignment_data(all_data, prepend_path=''):
    new_data = [None for _ in range(len(all_data))]
    for i, data in enumerate(all_data):
        if isinstance(data['model'], str) or isinstance(data['model'], unicode):
            model_file = os.path.join(prepend_path, data['model'])
            model_coords = None
        else:
            model_file = None
            model_coords = data['model']
        if isinstance(data['target'], str) or isinstance(data['target'], unicode):
            target_file = os.path.join(prepend_path, data['target'])
            target_coords = None
        else:
            target_file = None
            target_coords = data['target']

        new_data[i] = AlignedData(
            R=data['R'], T=data['T'], mapping=data['mapping'], inverted=data['inverted'], error=data['error_lsq'], swapped=data['swapped'],
            model_file=model_file, model_coords=model_coords, model_symbols=None, model_scale=data['model_rescale'],
            target_file=target_file, target_coords=target_coords, target_symbols=None, target_scale=data['target_rescale'],
            aligned_target_coords=data['aligned_target'], aligned_target_symbols=None
        )
    new_data = AlignedGroup(new_data)
    return new_data


class AlignedGroup(object):
    def __init__(self, data):
        self.data = data
        if not isinstance(data, list) or isinstance(data, np.ndarray):
            raise Exception("The input aligned data need to be contained in an ordered list.")
        self._combined = None
        self._average = None


    def __getitem__(self, key):
        return self.data[key]

    def __len__(self):
        return self.nclusters

    @property
    def nclusters(self):
        """ The number of clusters in this Group. """
        return len(self.data)


    @property
    def successful(self):
        """ A mask to operate only on successfully aligned clusters. """
        return np.array([cluster.successful for cluster in self.data])
    

    @property
    def coordinates(self):
        """ A numpy array of shape (3, nclusters) with the coordination positions of each atom for each cluster. """
        return np.array( [ [a.coord[0], a.coord[1], a.coord[2]] for a in cluster.atoms for cluster in self.clusters] ).T


    @property
    def natoms(self):
        """ Returns the (most common) number of atoms in the aligned clusters."""
        if self._natoms is None:
            # Check to see if every atom has the same number of atoms
            counter = Counter(datum.aligned_target.natoms for datum in self.data if datum.successful)
            self._natoms = counter.most_common(1)[0][0]
            if len(counter) != 1:
                warnings.warn("Not all clusters in this group have the same number of atoms. {0}\n  Ignoring those that comprise < 1% of the total.".format(counter), RuntimeWarning)
                total = sum(val for k,val in counter.items())
                for k,v in counter.items():
                    if v/float(total) < 0.01: counter[k] = False
        return self._natoms


    def average_structure(self, force_update=False):
        """ Calculates the atom positions of the average structure of the aligned clusters in the group. """
        if not force_update and hasattr(self, '_average') and self._average is not None:
            return self._average

        nsuccessful = sum(self.successful)
        avg_coords = np.zeros((self.natoms, 3), dtype=float)
        for aligned in self.data:
            if not aligned.successful: continue
            positions = aligned.aligned_target.positions#[aligned.mapping]
            if aligned.inverted:
                positions = positions * -1
            #print(positions[0])
            for i, position in enumerate(positions):
                avg_coords[i,0] += position[0,0]/nsuccessful
                avg_coords[i,1] += position[0,1]/nsuccessful
                avg_coords[i,2] += position[0,2]/nsuccessful

        self._average = Cluster(symbols=['Si' for i in range(self.natoms)], positions=avg_coords)
        return self._average


    def combine(self, colorize=True, force_update=False):
        """ Combines the group of aligned clusters into a single model. """
        pass
        """
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
        """
        pass



class AlignedData(object):
    def __init__(self, R, T, mapping, error, inverted, swapped=None,
                 model_file=None, model_coords=None, model_symbols=None, model_scale=1.0,
                 target_file=None, target_coords=None, target_symbols=None, target_scale=1.0,
                 aligned_target_coords=None, aligned_target_symbols=None):
        """
             model_file OR (model_coords AND model_symbols)
             target_file OR (target_coords AND target_symbols)
             aligned_target_coords, OPTIONAL aligned_target_symbols
        """
        self.R = np.matrix(R)
        self.T = np.matrix(T)

        self.mapping = np.array(mapping)
        self.error = error

        self._inverted = inverted
        self._swapped = swapped
        self.model_scale = model_scale
        self.target_scale = target_scale

        self._model_symbols = model_symbols
        self._model = np.matrix(model_coords) if model_coords else None
        if model_file is None and self._model is not None and len(self._model) == 3:
            self._model = self._model.T
        else:
            self._model_file = model_file

        self._target_symbols = target_symbols
        self._target = np.matrix(target_coords) if target_coords else None
        if target_file is None and self._target is not None and len(self._target) == 3:
            self._target = self._target.T
        else:
            self._target_file = target_file

        self._aligned_target_symbols = aligned_target_symbols
        self._aligned_target = np.matrix(aligned_target_coords)
        if len(self._aligned_target) == 3:
            self._aligned_target = self._aligned_target.T


    @property
    def successful(self):
        return max(Counter(self.mapping).values()) == 1


    @property
    def swapped(self):
        if self._swapped is None:
            self._swapped = len(self.target.positions) < len(self.model.positions)
        return self._swapped


    @lazyproperty
    def inverted(self):
        if self._inverted is None:
            included = np.array([i for i in range(np.amax(self.mapping)+1) if i in self.mapping])
            if self.swapped:
                model = self.model.positions[included]
            else:
                model = self.target.positions[included]
            l2 = Cluster._L2Norm(-model, self.aligned_target.positions)
            self._inverted = round(l2, 12) == round(self.error, 12)
        return self._inverted


    @lazyproperty
    def model(self):
        if self._model_file is None:
            print(self._model)
            c = Cluster(positions=self._model, symbols=self._model_symbols)
            del self._model
            del self._model_symbols
        else:
            #print("Loading model from file {}".format(self._model_file))
            c = Cluster(filename=self._model_file)
            del self._model_file
        #self._len_model_positions = len(c.positions)
        #if not hasattr(self, '_len_target_positions'):
        #    self.target
        #swapped = len(c.positions) < self._len_target_positions
        #if swapped and self.inverted:
        #    c._positions = -c._positions
        return c


    @lazyproperty
    def target(self):
        swapped = False
        if self._target_file is None:
            c = Cluster(positions=_target, symbols=self._target_symbols)
            del self._target
            del self._target_symbols
        else:
            #print("Loading target from file {}".format(self._target_file))
            c = Cluster(filename=self._target_file)
            del self._target_file
        #self._len_target_positions = len(c.positions)
        #if not hasattr(self, '_len_model_positions'):
        #    self.model
        #swapped = self._len_model_positions < len(c.positions)
        #if swapped and self.inverted:
        #    c._positions = -c._positions
        return c


    @property
    def aligned_model(self):
        """Returns the model coordinates that were used during alignment to compare to the rotated/translated target."""
        included = np.array([i for i in range(np.amax(self.mapping)+1) if i in self.mapping])
        indices = np.argsort(self.mapping)
        if self.swapped:
            #model = self.target.positions.copy()
            #model = self.target.positions[self.mapping].copy()
            model = self.target.positions[included].copy()
            #model = self.target.positions[indices].copy()
            scale = self.target_scale
        else:
            #model = self.model.positions.copy()
            #model = self.model.positions[self.mapping].copy()
            model = self.model.positions[included].copy()
            #model = self.model.positions[indices].copy()
            scale = self.model_scale
        Cluster._rescale_coordinates(model, 1.0/scale)
        if self.inverted:
            model = -model
        return model


    @property
    def aligned_target_symbols(self):
        indices = np.argsort(self.mapping)
        if self._aligned_target_symbols is not None:
            return self._aligned_target_symbols
        elif hasattr(self, '_target_symbols') and self._target_symbols is not None:
            return list(np.array(self._target_symbols)[indices])
        else:
            return list(np.array(self.target.symbols)[indices])

    @lazyproperty
    def aligned_target(self):
        # The coordinates from the aligned_target have been reorderd already, but the target has not been.
        # I want the target to stay looking exactly like the original target from the input file, so I do not want to reorder that.
        # Therefore I will un-reorder the aligned_target's data instead.
        indices = np.argsort(self.mapping)
        c = Cluster(positions=self._aligned_target[indices], symbols=self.aligned_target_symbols)  # Un-reorder the input data
        #c = Cluster(positions=self._aligned_target, symbols=self.aligned_target_symbols)  # Leave the input data in the order it was given
        del self._aligned_target
        return c


    def align_target(self, swapped=None, apply_mapping=True):
        # I don't fully understand why I need to use [indices] rather than [included] here.
        # I expected that the input data for the aligned_target coordinates have already been reorderd.
        # I unorder them, however, when I load initialize the datastructure.
        # There I expected that the target should not need to be reorderd, I should only need to extract out the specific atoms that were used, in order of their occurance in the original file.
        # This doesn't seem to be the case; the target evidently needs to be reorderd...?
        # However, the L2Norm function and the above do not seem to be consistent. Still, the code here gives me the correct results, so I am leaving it.

        #included = np.array([i for i in range(np.amax(self.mapping)+1) if i in self.mapping])
        indices = np.argsort(self.mapping)
        if swapped is None:
            swapped = self.swapped
        if not swapped:
            #target = self.target.positions
            #target = self.target.positions[included]  # I expected this one to be the one we needed to use, but evidently not...?
            #target = self.target.positions[indices]  # Use this one
            #target = self.target.positions[self.mapping]
            if apply_mapping:
                target = self.target.positions[indices]
            else:
                target = self.target.positions
            scale = self.target_scale
        else:
            #target = self.model.positions
            #target = self.model.positions[included]  # I expected this one to be the one we needed to use, but evidently not...?
            #target = self.model.positions[indices]  # Use this one
            #target = self.model.positions[self.mapping]
            if apply_mapping:
                target = self.model.positions[indices]
            else:
                target = self.model.positions
            scale = self.model_scale
        target /= scale
        return (self.R * target.T).T + np.tile((self.T), [target.shape[0], 1])


    def to_dict(self):
        d = {}
        d['R'] = self.R
        d['T'] = self.T
        d['mapping'] = self.mapping
        d['inverted'] = self.inverted
        d['error'] = self.error

        if hasattr(self, '_model_file'):
            d['model'] = self._model_file
        elif hasattr(self, '_model'):
            d['model'] = self._model
        else:
            d['model'] = self.model.to_dict()

        if hasattr(self, '_target_file'):
            d['target'] = self._target_file
        elif hasattr(self, '_target'):
            d['target'] = self._target
        else:
            d['target'] = self.target.to_dict()

        if hasattr(self, '_aligned_target'):
            d['aligned_target'] = self._aligned_target
        else:
            d['aligned_target'] = self.aligned_target.to_dict()
            del d['aligned_target']['symbols']
        return d


    def to_json(self, **kwargs):
        import json
        d = self.to_dict()
        for key, value in d.items():
            if isinstance(value, np.array(0).__class__):
                d[key] = value.tolist()
        return json.dumps(d, **kwargs)


    def L2Norm(self):
        return Cluster._L2Norm(self.aligned_model, self.aligned_target.positions)

    def L2Norm2(self):
        return Cluster._L2Norm2(self.aligned_model, self.aligned_target.positions)

    def L1Norm(self):
        return Cluster._L1Norm(self.aligned_model, self.aligned_target.positions)

    def LinfNorm(self):
        return Cluster._LinfNorm(self.aligned_model, self.aligned_target.positions)

    def angular_variation(self, neighbor_cutoff):
        return Cluster._angular_variation(self.aligned_model, self.aligned_target.positions, neighbor_cutoff)



class Cluster(object):
    def __init__(self, symbols=None, positions=None, rescaling_constant=None, center_atom=None, filename=None):
        self.filename = filename
        if filename is None:
            assert len(symbols) == len(positions)
            self._positions = np.matrix(positions).copy()
            self.natoms = len(symbols)
            self.atoms = [Atom(symbol=symbol, position=position, id=i) for i,(symbol,position) in enumerate(zip(symbols, self._positions))]
            self.rescaling_constant = rescaling_constant
            self.center_atom = center_atom
        else:
            with open(filename) as f:
                self.natoms = int(f.readline().strip())
                self.comment = f.readline().strip()
                self._positions = np.matrix(np.zeros((self.natoms, 3), dtype=float))
                self.atoms = []
                self.rescaling_constant = rescaling_constant
                self.center_atom = center_atom
                for i in range(self.natoms):
                    znum, x,y,z = tuple(f.readline().strip().split())
                    x = float(x); y = float(y); z = float(z)
                    self._positions[i,0] = x
                    self._positions[i,1] = y
                    self._positions[i,2] = z
                    atom = Atom(symbol=znum, position=self._positions[i], id=i)
                    self.atoms.append(atom)

    def __getitem__(self, index):
        return self.atoms[i]

    def __len__(self):
        return self.natoms

    @property
    def symbols(self):
        return [atom.symbol for atom in self]


    @property
    def positions(self):
        return self._positions


    def __getitem__(self, index):
        return self.atoms[index]


    def to_dict(self):
        d = {}
        if getattr(self, 'filename', None) is not None:
            d['filename'] = self.filename
        else:
            d['positions'] = self._positions.tolist()
            d['symbols'] = self.symbols
            d['rescaling_constant'] = self.rescaling_constant
        return d


    def to_json(self, **kwargs):
        import json
        return json.dumps(self.to_dict(), **kwargs)


    def to_xyz(self, filename=None, comment='', xsize=None, ysize=None, zsize=None):
        lines = []
        lines.append('{}\n'.format(self.natoms))
        comment = '{} {} {} {}'.format(xsize or '', ysize or '', zsize or '', comment).strip()
        lines.append('{}\n'.format(comment))
        for atom in self.atoms:
            lines.append('{sym} {x} {y} {z}\n'.format(sym=atom.symbol, x=atom.x, y=atom.y, z=atom.z))
        if filename is not None:
            of = open(filename, 'w')
            for line in lines:
                of.write(line)
            of.close()
        return ''.join(lines)


    @staticmethod
    def _rescale_coordinates(coordinates, scale):
        coordinates *= scale


    def rescale_edge_lengths(self, avg_edge_len, verbose=False):
        pdists = scipy.spatial.distance.pdist(coordinates)
        mean = np.mean(pdists)
        Cluster._rescale_coordinates(self._positions, 1.0/mean*avg_edge_len)
        if verbose:
            print("Rescaled by {}".format(scale))
        return scale


    def normalize_edge_lengths(self, verbose=False):
        pdists = scipy.spatial.distance.pdist(coordinates)
        mean = np.mean(pdists)
        return self.rescale_edge_lengths(1.0/mean, verbose)


    @staticmethod
    def _L2Norm(c1, c2):
        assert len(c1) == len(c2)
        L2 = 0.0
        for i in range(len(c1)):
            L2 += ((c1[i,0] - c2[i,0])**2 + 
                   (c1[i,1] - c2[i,1])**2 + 
                   (c1[i,2] - c2[i,2])**2)
        return math.sqrt(L2)/len(c1)


    @staticmethod
    def _L2Norm2(c1, c2):
        assert len(c1) == len(c2)
        L2 = 0.0
        for i in range(len(c1)):
            L2 += ((c1[i,0] - c2[i,0])**2 + 
                   (c1[i,1] - c2[i,1])**2 + 
                   (c1[i,2] - c2[i,2])**2)
        return L2/len(c1)


    @staticmethod
    def _L1Norm(c1, c2):
        assert len(c1) == len(c2)
        L1 = 0.0
        for i in range(len(c1)):
            L1 += abs(c1[i,0] - c2[i,0]) + abs(c1[i,1] - c2[i,1]) + abs(c1[i,2] - c2[i,2])
        return L1/len(c1)


    @staticmethod
    def _LinfNorm(c1, c2):
        assert len(c1) == len(c2)
        Linf = [0.0 for i in range(len(c1))]
        for i in range(len(c1)):
            Linf[i] += (abs(c1[i,0] - c2[i,0]) +
                        abs(c1[i,1] - c2[i,1]) +
                        abs(c1[i,2] - c2[i,2]))
        return max(Linf)/len(c1)


    @staticmethod
    def _angle(a, b, c):
        a = np.array(a)
        b = np.array(b)
        c = np.array(c)
        ba = a - b
        bc = c - b
        if (ba == bc).all():
            return 0.0
        cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
        angle = np.arccos(cosine_angle)
        return np.degrees(angle)


    @staticmethod
    def _angular_variation(c1, c2, neighbor_cutoff):
        """This is a measure of the angular variation between the two clusters.
        Calculates Mean( abs(Delta(angle btwn all atoms and their neighbors, going through the center)) ).
        Both clusters must have had their center removed.
        Both clusters are assumed to have had their center at (0,0,0) before removal.
        Both cluster must have been reordered so that self[i] corresponds directly to other[i]."""

        #if c1.center_atom is not None or c2.center_atom is not None:
        #    raise NotImplementedError("The centers must be removed to use this method. It's also necessary that the cluster is recentered so that the center atom was at (0,0,0) before running this method.")
        natoms = min(len(c1), len(c2))

        # Generate neighbors using pdist and the array for both c1 and c2
        print((scipy.spatial.distance.pdist(c1) + scipy.spatial.distance.pdist(c2)) /2.0)
        c1_dist_matrix = scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(c1))
        c2_dist_matrix = scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(c2))
        dist_matrix = (c1_dist_matrix + c2_dist_matrix) / 2.0  # Calculate the average atom-to-atom distances
        neighbors = []
        for index, x in np.ndenumerate(dist_matrix):
            if x < neighbor_cutoff:
                neighbors.append(index)

        # Use those neighbors to calculate angles for c1 and c2 simultaneously, and calculate abs(Delta(angles))
        # Also calculate and return the mean
        mean = 0.
        for neighbor_pair in neighbors:
            i,j = neighbor_pair
            xi1,yi1,zi1 = c1[i,0], c1[i,1], c1[i,2]
            xi2,yi2,zi2 = c2[i,0], c2[i,1], c2[i,2]
            xj1,yj1,zj1 = c1[j,0], c1[j,1], c1[j,2]
            xj2,yj2,zj2 = c2[j,0], c2[j,1], c2[j,2]
            a1 = Cluster._angle((xi1,yi1,zi1), (0,0,0), (xj1,yj1,zj1))
            a2 = Cluster._angle((xi2,yi2,zi2), (0,0,0), (xj2,yj2,zj2))
            mean += abs(a1-a2)
        mean /= len(neighbors)
        return mean



class Atom(object):
    symbols = {1:"H", 2:"He", 3:"Li", 4:"Be", 5:"B", 6:"C", 7:"N", 8:"O", 9:"F", 10:"Ne", 11:"Na", 12:"Mg", 13:"Al", 14:"Si", 15:"P", 16:"S", 17:"Cl", 18:"Ar", 19:"K", 20:"Ca", 21:"Sc", 22:"Ti", 23:"V", 24:"Cr", 25:"Mn", 26:"Fe", 27:"Co", 28:"Ni", 29:"Cu", 30:"Zn", 31:"Ga", 32:"Ge", 33:"As", 34:"Se", 35:"Br", 36:"Kr", 37:"Rb", 38:"Sr", 39:"Y", 40:"Zr", 41:"Nb", 42:"Mo", 43:"Tc", 44:"Ru", 45:"Rh", 46:"Pd", 47:"Ag", 48:"Cd", 49:"In", 50:"Sn", 51:"Sb", 52:"Te", 53:"I", 54:"Xe", 55:"Cs", 56:"Ba", 57:"La", 58:"Ce", 59:"Pr", 60:"Nd", 61:"Pm", 62:"Sm", 63:"Eu", 64:"Gd", 65:"Tb", 66:"Dy", 67:"Ho", 68:"Er", 69:"Tm", 70:"Yb", 71:"Lu", 72:"Hf", 73:"Ta", 74:"W", 75:"Re", 76:"Os", 77:"Ir", 78:"Pt", 79:"Au", 80:"Hg", 81:"Tl", 82:"Pb", 83:"Bi", 84:"Po", 85:"At", 86:"Rn", 87:"Fr", 88:"Ra", 89:"Ac", 90:"Th", 91:"Pa", 92:"U", 93:"Np", 94:"Pu", 95:"Am", 96:"Cm", 97:"Bk", 98:"Cf", 99:"Es", 100:"Fm", 101:"Md", 102:"No", 103:"Lr", 104:"Rf", 105:"Db", 106:"Sg", 107:"Bh", 108:"Hs", 109:"Mt", 110:"Ds", 111:"Rg", 112:"Cn", 113:"Uut", 114:"Fl", 115:"Uup", 116:"Lv", 117:"Uus", 118:"Uu"}
    numbers = {v:k for k,v in symbols.items()}

    def __init__(self, symbol, position, id=None):
        self.id = id
        self._symbol = symbol
        self.position = position

    @property
    def symbol(self):
        if isinstance(self._symbol, str):
            return self._symbol
        elif isinstance(self._symbol, int):
            return Atom.symbols[self._symbol]

    @property
    def number(self):
        if isinstance(self._symbol, int):
            return self._symbol
        elif isinstance(self._symbol, str):
            return Atom.numbers[self._symbol]

    @property
    def x(self):
        return self.position[0,0]

    @property
    def y(self):
        return self.position[0,1]

    @property
    def z(self):
        return self.position[0,2]



