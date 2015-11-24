import sys
from model import Model
from atom import Atom


def create_model(translation_vectors, comment, lattice_parameter, units=1, atom_type='Si', filename=None):
    a = lattice_parameter
    eps = 1e-6
    points = []
    for x in range(units+1):
        for y in range(units+1):
            for z in range(units+1):
                for tv in translation_vectors:
                    point = ( (tv[0]+x)*a, (tv[1]+y)*a, (tv[2]+z)*a )
                    if point[0] <= units*a+eps and point[1] <= units*a+eps and point[2] <= units*a+eps:
                        points.append(point)

    atoms = [Atom(i, atom_type, x,y,z) for i,(x,y,z) in enumerate(points)]
    m = Model(comment=comment, xsize=units*a, ysize=units*a, zsize=units*a, atoms=atoms)
    m.recenter()
    return m

def simple_cubic(lattice_parameter, units=1, atom_type='Si', filename=None):
    comment = 'Simple cubic lattice'
    translation_vectors = [(0.0, 0.0, 0.0)]
    return create_model(translation_vectors=translation_vectors, comment=comment, lattice_parameter=lattice_parameter, units=units, atom_type=atom_type, filename=filename)

def fcc(lattice_parameter, units=1, atom_type='Si', filename=None):
    comment = 'FCC lattice'
    translation_vectors = [(0.0, 0.0, 0.0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5)]
    return create_model(translation_vectors=translation_vectors, comment=comment, lattice_parameter=lattice_parameter, units=units, atom_type=atom_type, filename=filename)

def bcc(lattice_parameter, units=1, atom_type='Si', filename=None):
    comment = 'BCC lattice'
    translation_vectors = [(0.0, 0.0, 0.0), (0.5, 0.5, 0.5)]
    return create_model(translation_vectors=translation_vectors, comment=comment, lattice_parameter=lattice_parameter, units=units, atom_type=atom_type, filename=filename)


def main():
    m = simple_cubic(4.050, units=2)
    m = fcc(4.050, units=2)
    m = bcc(4.050, units=2)

    m.save()

if __name__ == '__main__':
    main()
