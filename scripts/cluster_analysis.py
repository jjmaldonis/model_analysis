
import sys
from model import Model
from hutch import Hutch
from atom import Atom
from atom_graph import AtomGraph
import numpy as np


def nnnic(atom,model,cluster):
    """ number of nearest neighbors in cluster 
    @param atom is the atom around which we want to find its neighbors
    @param model is the model that the cluster belongs to
    @param cluster is the cluster of atoms, in list form """
    n = 0
    for atomi in atom.neighs:
        if atomi in cluster:
            n += 1
    return n
    

def main():
    # sys.argv[1] should be the modelfile in the xyz format
    # sys.argv[2] should be the cutoff desired
    modelfile = sys.argv[1]
    cutoff = float(sys.argv[2])
    ag = AtomGraph(modelfile,cutoff)
    cluster_prefix = ''
    #clusters = ag.get_clusters('Icosahedra-like','Mixed')
    #clusters = ag.get_clusters('Icosahedra-like')
    #cluster_prefix = 'ico.golden.t3.'
    #clusters = ag.get_clusters('Full-icosahedra')
    #cluster_prefix = 'fullico.golden.t3.'
    clusters = ag.get_clusters('Crystal-like')
    cluster_prefix = 'xtal.t3.'


    ## Let's just get the biggest one for now.
    ## Comment this section to get them all.
    #l = 0
    #for cluster in clusters:
    #    if(len(cluster) > l):
    #        l = len(cluster)
    #        i = clusters.index(cluster)
    #clusters = [clusters[i]]

    orig_clusters = clusters[:]
    # Print orig clusters
    j = 0
    for i,cluster in enumerate(clusters):
        print("Orig cluster {0} contains {1} atoms.".format(i,len(cluster)))
        #for atom in cluster:
        #    if atom not in orig_clusters[i]:
        #        #print("{0} {1}".format(atom.vesta(),atom.vp))
        #        print(atom.vesta())
        #        j += 1
        #    else:
        #        print(atom.vesta())
        # Save cluster files
        cluster_model = Model("Orig cluster {0} contains {1} atoms.".format(i,len(cluster)),ag.model.lx, ag.model.ly, ag.model.lz, cluster)
        cluster_model.write_cif('{1}cluster{0}.cif'.format(i,cluster_prefix))
        cluster_model.write_our_xyz('{1}cluster{0}.xyz'.format(i,cluster_prefix))

    # This changes all the atoms in the original cluster
    # (ie of type set above) to have no atomic number, 
    # which displays as Te (gold) in vesta for me
    #for i,cluster in enumerate(clusters):
    #    for atom in cluster:
    #        atom.z = 0

    # Expand the clusters
    # This overwrites the old ones
    # Include all atoms already in the cluster and all atoms that have at least *5* neighbors in the cluster
    for i,cluster in enumerate(clusters):
        new_clust = [ atom for atom in ag.model.atoms if atom in clusters[i] or ( nnnic(atom,ag.model,clusters[i]) >= 5 ) ]
        new_clust = list(set(new_clust))
        clusters[i] = new_clust[:]

    j = 0
    for i,cluster in enumerate(clusters):
        print("New, bigger cluster {0} contains {1} atoms.".format(i,len(cluster)))
        #for atom in cluster:
        #    if atom not in orig_clusters[i]:
        #        #print("{0} {1}".format(atom.vesta(),atom.vp))
        #        print(atom.vesta())
        #        j += 1
        #    else:
        #        print(atom.vesta())
        # Save cluster files
        cluster_model = Model("New, bigger cluster {0} contains {1} atoms.".format(i,len(cluster)),ag.model.lx, ag.model.ly, ag.model.lz, cluster)
        cluster_model.write_cif('{1}cluster{0}.big.cif'.format(i,cluster_prefix))
        cluster_model.write_our_xyz('{1}cluster{0}.big.xyz'.format(i,cluster_prefix))

    #largest = max([x for x in clusters],key=len)
    #cluster_model = Model("New,bigger cluster contains {0} atoms.".format(len(largest)),ag.model.lx, ag.model.ly, ag.model.lz, largest)
    #cluster_model.write_cif('cluster.largest.cif')
    allclusters = np.array(clusters).ravel()
    allclusters = list(set([atom for cluster in clusters for atom in cluster]))
    savecluster = Model("All clusters.",ag.model.lx, ag.model.ly, ag.model.lz, allclusters)
    savecluster.write_cif('{0}allclusters.cif'.format(cluster_prefix))
    savecluster.write_our_xyz('{0}allclusters.xyz'.format(cluster_prefix))
    print("{0}allclusters.cif and {0}allclusters.xyz contain {1} atoms.".format(cluster_prefix, savecluster.natoms))

    x_cluster = ag.model.atoms[:]
    to_remove = []
    for i,atom in enumerate(x_cluster):
        if atom in allclusters:
            to_remove.append(i)
    to_remove.reverse()
    for i in to_remove:
        x_cluster.pop(i)
    x_cluster = Model("Opposite cluster of {0}".format(cluster_prefix),ag.model.lx, ag.model.ly, ag.model.lz, x_cluster)
    x_cluster.write_cif('{0}opposite.cif'.format(cluster_prefix))
    x_cluster.write_our_xyz('{0}opposite.xyz'.format(cluster_prefix))
    print('{0}opposite.cif and {0}opposite.xyz contain {1} atoms.'.format(cluster_prefix, x_cluster.natoms))


                


if __name__ == "__main__":
    main()
