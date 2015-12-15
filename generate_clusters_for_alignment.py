""" Pseudocode:
1) this script takes in a model (specified in-line not on the command line)
2) randomly picks an atom in that model
3) creates a nearest-neighbor cluster around that atom
4) fixes any periodic boundary problems where the cluster wraps around the edge of the model
5) normalizes the average bond distance of the cluster to be 1.0
6) saves the cluster to a directory of your choosing
7) steps 2-6 are repeated 'num_clusters' times, and no atom can be selected twice

These clusters are used as input for Arash's rotation alignment code.
"""

import sys, os, random
from model import Model
from rotation_alignment_analysis import Cluster

def main():
    # Load the MD model
    md_modelfile = 'model_0k.xyz'
    md = Model(md_modelfile)

    # Load the cutoff dictionary so that we can generate neighbors for every atom
    from cutoff import cutoff

    # Directory name to save the files to
    dir = 'clusters'

    # Set the number of clusters to randomly select
    num_clusters = 10000

    # Make a list of which atoms we have selected already
    picked = [False for i in range(md.natoms)]

    for n in range(num_clusters):
        while True:
            # Keep generating rand numbers between 0 and the total number of atoms until we get an atom that we haven't used yet
            rand = random.randrange(0, md.natoms)
            if not picked[rand]:
                # Set picked to True
                picked[rand] = True

                # Generate neighbors for this atom
                atom = md.atoms[rand]
                atom.neighs = md.get_atoms_in_cutoff(atom, cutoff)
                if(atom in atom.neighs): atom.neighs.remove(atom)
                atom.cn = len(atom.neighs)

                # Create the cluster, normalize the bond distances, and write to disk
                atoms = md.atoms[rand].neighs + [md.atoms[rand]]
                c = Cluster(
                    comment='cluster #{0} from atom {1}'.format(n, rand),
                    xsize=5.0,
                    ysize=5.0,
                    zsize=5.0,
                    atoms=atoms,
                    center_included = True
                )
                c.fix_cluster_pbcs()
                ratio = c.normalize_bond_distances()
                c.comment='cluster #{0} from atom {1}; normalized bond distances by {2}'.format(n, rand, ratio)
                c.write(os.path.join(dir, '{0}.xyz'.format(n)))
                break
        print(n)


if __name__ == '__main__':
    main()
