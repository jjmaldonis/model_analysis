import os
import math
import gzip, json
import numpy as np
import scipy.spatial.distance
from collections import Counter
import warnings

from lazy_property import lazyproperty



def load_gzipped_json(filename, verbose=True):
    with gzip.open(filename, "rb") as f:
        d = json.loads(f.read().decode("ascii"))
    if verbose:
        print("Loaded gzip successfully.")
    return d

# Use this one externally
def load_alignment_data(filename, prepend_path='', prepend_model_path='', prepend_target_path='', verbose=True):
    all_data = load_gzipped_json(filename, verbose)
    return render_alignment_data(all_data, prepend_path, prepend_model_path, prepend_target_path, verbose)

def render_alignment_data(all_data, prepend_path='', prepend_model_path='', prepend_target_path='', verbose=True):
    if prepend_path:
        prepend_model_path = prepend_path
        prepend_target_path = prepend_path
    new_data = [None for _ in range(len(all_data))]
    for i, data in enumerate(all_data):
        if isinstance(data['model'], str) or isinstance(data['model'], unicode):
            model_file = os.path.join(prepend_model_path, data['model'])
            model_coords = None
        else:
            model_file = None
            model_coords = data['model']
        if isinstance(data['target'], str) or isinstance(data['target'], unicode):
            target_file = os.path.join(prepend_target_path, data['target'])
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
    if verbose:
        print("Rendered data successfully")
    return new_data



class Positions(np.matrix):
    atom_types = ['Si', 'Na', 'Mg', 'Ti', 'V', 'Be', 'Mn', 'Fe', 'P', 'Ni', 'Cu', 'S', 'B', 'He', 'Ga', 'C', 'Sn', 'Pb', 'O']

    def __init__(self, *args, **kwargs):
        super(Positions, self).__init__(*args, **kwargs)
        self = self.view(Positions)
        assert isinstance(self, Positions)


    def apply_transformation(self, R, T, invert=False):
        if not invert:
            this = (R * self.transpose()).transpose() + np.tile(T, [self.shape[0], 1])
        else:
            this = (R * (self + np.tile(T, [self.shape[0], 1])).transpose()).transpose()
        return this.view(Positions)


    def to_xyz(self, symbols=None, comment=''):
        if symbols is None:
            symbols = ['Si' for _ in range(len(self))]
        elif isinstance(symbols, str):
                symbols = [symbols for _ in range(len(self))]
        lines = []
        natoms = len(self)
        lines.append(str(natoms))
        lines.append(comment)
        for i, row in enumerate(self):
            lines.append('{} {} {} {}'.format(symbols[i], row[0,0], row[0,1], row[0,2]))
        return '\n'.join(lines)


    def to_colorized_xyz(self, comment=''):
        lines = []
        natoms = len(self)
        lines.append(str(natoms))
        lines.append(comment)
        for i, row in enumerate(self):
            lines.append('{} {} {} {}'.format(Positions.atom_types[i], row[0,0], row[0,1], row[0,2]))
        return '\n'.join(lines)



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
        return Positions( [ [a.coord[0], a.coord[1], a.coord[2]] for a in cluster.atoms for cluster in self.clusters] ).T


    def average_structure(self, force_update=False):
        """ Calculates the atom positions of the average structure of the aligned clusters in the group. """
        if not force_update and hasattr(self, '_average') and self._average is not None:
            return self._average

        nsuccessful = sum(self.successful)
        natoms = np.max([a.target.natoms for a in self])
        avg_coords = np.zeros((natoms, 3), dtype=float)
        natoms_per_index = np.zeros((natoms,), dtype=int)
        for aligned in self:
            if not aligned.successful: continue
            target, model = aligned.rotate_target_onto_model(apply_mapping=True)
            natoms = min(len(target), len(model))
            model = model[:natoms]
            target = target[:natoms]
            for i, coord in enumerate(target):
                avg_coords[i,0] += coord[0,0]
                avg_coords[i,1] += coord[0,1]
                avg_coords[i,2] += coord[0,2]
                natoms_per_index[i] += 1
        for i, row in enumerate(avg_coords):
            if natoms_per_index[i] == 0:
                natoms_per_index[i] = 1
            row /= natoms_per_index[i]  # This modifies in-place

        self._average = Cluster(symbols=['Si' for i in range(len(avg_coords))], positions=avg_coords)
        return self._average


    def combine(self, colorize=True, force_update=False):
        """ Combines the group of aligned clusters into a single model. """
        if not force_update and hasattr(self, '_combined') and self._combined is not None:
            return self._combined

        #              0     1     2     3     4     5    6     7     8     9    10     11    12    13   14    15    16   17    18
        atom_types = ['Si', 'Na', 'Mg', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'B', 'Al', 'Ga', 'C', 'Sn', 'Pb', 'O']

        nsuccessful = sum(self.successful)
        coords = []
        for aligned in self.data:
            if not aligned.successful: continue
            target, model = aligned.rotate_target_onto_model(apply_mapping=True)
            natoms = min(len(target), len(model))
            model = model[:natoms]
            target = target[:natoms]
            for i, coord in enumerate(target):
                #coords.append((atom_types[len(target)], coord))
                coords.append((atom_types[i], coord))
        if colorize:
            symbols = [sym for sym, coord in coords]
        else:
            symbols = ['Si' for sym, coord in coords]
        #coords = Positions([coord for sym, coord in coords])
        coords = np.array([coord for sym, coord in coords])
        coords = Positions(coords)
        self._combined = Cluster(symbols=symbols, positions=coords)
        return self._combined



class AlignedData(object):
    def __init__(self, R, T, mapping, error, inverted, swapped,
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
        self._swapped = not swapped # TODO
        self.model_scale = model_scale
        self.target_scale = target_scale

        self._model_symbols = model_symbols
        self._model = Positions(model_coords) if model_coords else None
        if model_file is None and self._model is not None and len(self._model) == 3:
            self._model = self._model.T
        else:
            self._model_file = model_file

        self._target_symbols = target_symbols
        self._target = Positions(target_coords) if target_coords else None
        if target_file is None and self._target is not None and len(self._target) == 3:
            self._target = self._target.T
        else:
            self._target_file = target_file

        self._aligned_target_symbols = aligned_target_symbols
        self._aligned_target = Positions(aligned_target_coords)
        if len(self._aligned_target) == 3:
            self._aligned_target = self._aligned_target.T


    @property
    def successful(self):
        return max(Counter(self.mapping).values()) == 1


    @property
    def swapped(self):
        return self._swapped


    @property
    def inverted(self):
        return self._inverted


    @lazyproperty
    def model(self):
    #def __model(self):
        if self._model_file is None:
            c = Cluster(positions=self._model, symbols=self._model_symbols)
            del self._model
            del self._model_symbols
        else:
            c = Cluster(filename=self._model_file)
            del self._model_file
        return c


    @lazyproperty
    def target(self):
    #def __target(self):
        if self._target_file is None:
            c = Cluster(positions=_target, symbols=self._target_symbols)
            del self._target
            del self._target_symbols
        else:
            c = Cluster(filename=self._target_file)
            del self._target_file
        return c


    #@lazyproperty
    #def model(self):
    #    if not self.swapped:
    #        return self.__model()
    #    else:
    #        return self.__target()


    #@lazyproperty
    #def target(self):
    #    if not self.swapped:
    #        return self.__target()
    #    else:
    #        return self.__model()


    @property
    def aligned_model(self):
        """Returns the model coordinates that were used during alignment to compare to the rotated/translated target."""
        if self.swapped:
            model = self.target.positions
            scale = self.target_scale
        else:
            model = self.model.positions
            scale = self.model_scale
        #included = np.array([i for i in range(np.amax(self.mapping)+1) if i in self.mapping])
        indices = np.argsort(self.mapping)
        #model = model.copy()
        #model = model[self.mapping].copy()
        #model = model[included].copy()  # Use this one
        model = model[indices].copy()
        model = Cluster._rescale_coordinates(model, 1.0/scale)
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


    def align_target(self, rescale=True, apply_mapping=True):
        # I don't fully understand why I need to use [indices] rather than [included] here.
        # I expected that the input data for the aligned_target coordinates have already been reorderd.
        # I unorder them, however, when I load initialize the datastructure.
        # Therefore I expected that the target should not need to be reorderd, I should only need to extract out the specific atoms that were used, in order of their occurance in the original file.
        # This doesn't seem to be the case; the target evidently needs to be reorderd...?
        # However, the L2Norm function and the above do not seem to be consistent. Still, the code here gives me the correct results, so I am leaving it.

        if not self.swapped:
            target = self.target.positions.copy()
            scale = self.target_scale
        else:
            target = self.model.positions.copy()
            scale = self.model_scale
        if apply_mapping:
            included = np.array([i for i in range(np.amax(self.mapping)+1) if i in self.mapping])
            indices = np.argsort(self.mapping)
            #target = target
            target = target[included]  # I expected this one to be the one we needed to use, but evidently not...?
            #target = target[indices]  # Use this one
            #target = target[self.mapping]
        if rescale:
            target /= scale
        target = target.apply_transformation(self.R, self.T)
        return target


    # Either rotate_target_onto_model or rotate_model_onto_target will reproduce aligned_target.
    # It depends on whether swapped is true or false
    # if not swapped (which is equivalent to swapped is False) then use rotate_target_onto_model
    # if swapped, then use rotate_model_onto_target to reproduce aligned_target

    def rotate_target_onto_model(self, rescale=True, apply_mapping=True):
        included = np.array([i for i in range(np.amax(self.mapping)+1) if i in self.mapping])
        indices = np.argsort(self.mapping)

        # Rescale, apply mapping to, and tranform the target
        if rescale:
            target = self.target.rescale_edge_lengths(1.0/self.target_scale)
        else:
            target = self.target.positions.copy()
        if not self.swapped:
            if apply_mapping:
                target = target[indices]
            R = self.R
            T = self.T
            target = target.apply_transformation(R,T, invert=False)
        else:
            if self.inverted:
                target = -target
            if apply_mapping:
                target = target[self.mapping]
            R = np.linalg.inv(self.R)
            T = -self.T
            target = target.apply_transformation(R,T, invert=True)

        # Rescale and apply mapping to the model
        if rescale:
            model = self.model.rescale_edge_lengths(1.0/self.model_scale)
        else:
            model = self.model.positions.copy()
        if not self.swapped:
            if self.inverted:
                model = -model
            # This is weird; we want to include all the atoms even though in this case there are more
            # atoms in the model and the target. However, only some of the atoms in the model were
            # mapped to the target, so in order to color-compare the atoms between the targe and model,
            # we need to put the atoms in the model that were not included in the alignment on the end
            # of the positions array.
            not_included = np.array([i for i in range(len(model)) if i not in included])
            all = np.append(included, not_included)
            model = model[all]
        else:
            model = model

        if not self.swapped and self.inverted:
            target, model = -target, -model
        return target, model


    def rotate_model_onto_target(self, rescale=True, apply_mapping=True):
        included = np.array([i for i in range(np.amax(self.mapping)+1) if i in self.mapping])
        indices = np.argsort(self.mapping)

        # Rescale, apply mapping to, and tranform the model
        if rescale:
            model = self.model.rescale_edge_lengths(1.0/self.model_scale)
        else:
            model = self.model.positions.copy()
        if not self.swapped:
            if self.inverted:
                model = -model
            if apply_mapping:
                model = model[self.mapping]
            R = np.linalg.inv(self.R)
            T = -self.T
            model = model.apply_transformation(R,T, invert=True)
        else:
            if apply_mapping:
                model = model[indices]
            R = self.R
            T = self.T
            model = model.apply_transformation(R,T, invert=False)

        # Rescale and apply mapping to the target
        if rescale:
            target = self.target.rescale_edge_lengths(1.0/self.target_scale)
        else:
            target = self.target.positions.copy()
        if not self.swapped:
            target = target
        else:
            if self.inverted:
                target = -target
            # This is weird; we want to include all the atoms even though in this case there are more
            # atoms in the model and the target. However, only some of the atoms in the model were
            # mapped to the target, so in order to color-compare the atoms between the targe and model,
            # we need to put the atoms in the model that were not included in the alignment on the end
            # of the positions array.
            not_included = np.array([i for i in range(len(model)) if i not in included])
            all = np.append(included, not_included)
            target = target[all]

        if self.swapped and self.inverted:
            target, model = -target, -model
        return model, target


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
            self._positions = Positions(positions).copy()
            self._natoms = len(symbols)
            self.atoms = [Atom(symbol=symbol, position=position, id=i) for i,(symbol,position) in enumerate(zip(symbols, self._positions))]
            self.rescaling_constant = rescaling_constant
            self.center_atom = center_atom
        else:
            with open(filename) as f:
                self._natoms = int(f.readline().strip())
                self.comment = f.readline().strip()
                self._positions = Positions(np.zeros((self.natoms, 3), dtype=float))
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
    def natoms(self):
        return self._natoms


    @property
    def CN(self):
        center = self.find_center()
        if center is None:
            has_center = False
        else:
            has_center = True
        return self._natoms - has_center


    @property
    def positions(self):
        return self._positions


    def find_center(self):
        """This is beautiful."""
        if hasattr(self, '_center'):
            return self._center
        dist_matrix = scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(self.positions))
        dist_matrix /= np.mean(dist_matrix)
        normalized = np.divide(np.std(dist_matrix, axis=0), np.std(dist_matrix))
        inverse_normalized = np.abs(1.0 - 1.0/normalized)
        for i, value in enumerate(inverse_normalized):
            # i.e. if one atom's distances to all other atoms is at least 0.5 stdevs more than all the others, return it
            if value > 0.5:
                self._center = self[i]
                break
        else:
            self._center = None
        return self._center


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


    def vp_index(self):
        if not hasattr(self, '_vp_index'):
            from voropp import compute_index
            self._vp_index = compute_index(self.filename)
        return self._vp_index


    @staticmethod
    def _rescale_coordinates(coordinates, scale):
        return coordinates * scale


    def rescale_edge_lengths(self, scale, verbose=False):
        if verbose:
            print("Rescaled by {}".format(scale))
        return Cluster._rescale_coordinates(self.positions, scale)


    def normalize_edge_lengths(self, verbose=False):
        pdists = scipy.spatial.distance.pdist(self.positions)
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



