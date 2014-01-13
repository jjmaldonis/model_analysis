
import sys
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
    clusters = ag.get_clusters("Crystal-like")

    ## Let's just get the biggest one for now.
    ## Comment this section to get them all.
    #l = 0
    #for cluster in clusters:
    #    if(len(cluster) > l):
    #        l = len(cluster)
    #        i = clusters.index(cluster)
    #clusters = [clusters[i]]

    orig_clusters = clusters[:]
    for i,cluster in enumerate(clusters):
        new_clust = []
        for atom in ag.model.atoms:
            if( nnnic(atom,ag.model,clusters[i]) >= 4 ):
                new_clust.append(atom)
        for atom in clusters[i]:
            if atom not in new_clust:
                new_clust.append(atom)
        clusters[i] = new_clust[:]

    j = 0
    for i,cluster in enumerate(clusters):
        print("New, bigger cluster {0} contains {1} atoms.".format(i,len(cluster)))
        for atom in cluster:
            if atom not in orig_clusters[i]:
                print("{0} {1}".format(atom.vesta(),atom.vp))
                j += 1
            else:
                print(atom.vesta())
    print(j)

                


if __name__ == "__main__":
    main()
