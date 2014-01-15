
import sys
from model import Model
from hutch import Hutch
from atom import Atom
from atom_graph import AtomGraph


def nnnic(atom,model,cluster):
    """ number of nearest neighbors in cluster 
    @param atom is the atom around which we want to find its neighbors
    @param model is the model that the cluster belongs to
    @param cluster is the cluster of atoms, in list form """
    n = 0
    for atomi in model.neighs[atom]:
        if atomi in cluster:
            n += 1
    return n
    

def main():
    # sys.argv[1] should be the modelfile in the xyz format
    # sys.argv[2] should be the cutoff desired
    modelfile = sys.argv[1]
    cutoff = float(sys.argv[2])
    ag = AtomGraph(modelfile,cutoff)
    clusters = ag.get_clusters("Icosahedra-like")

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
        cluster_model.write_cif('cluster{0}.cif'.format(i))

    # Expand the clusters
    for i,cluster in enumerate(clusters):
        new_clust = [ atom for atom in ag.model.atoms if atom in clusters[i] or ( nnnic(atom,ag.model,clusters[i]) >= 4 ) ]
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
        cluster_model.write_cif('cluster{0}.big.cif'.format(i))

    #largest = max([x for x in clusters],key=len)
    #cluster_model = Model("New,bigger cluster contains {0} atoms.".format(len(largest)),ag.model.lx, ag.model.ly, ag.model.lz, largest)
    #cluster_model.write_cif('cluster.largest.cif')

                


if __name__ == "__main__":
    main()
