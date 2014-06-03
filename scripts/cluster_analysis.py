import sys
from model import Model
from hutch import Hutch
from atom import Atom
from atom_graph import AtomGraph
from vor import Vor
import categorize_vor
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
    #for atom in ag.model.atoms:
    #    if(atom.z == 40):
    #        print(atom,atom.cn,atom.neighs,atom.vp.index,atom.vp.type)
    model = Model(modelfile)
    model.generate_neighbors(cutoff)

    #j = 0
    #atoms = []
    #for i,atom in enumerate(ag.model.atoms):
    #    if(atom.vp.index[0:3] == [0,3,6] and i > 200):
    #    #if(atom.vp.index[0:3] == [0,3,6]):
    #        atom.z = 0
    #        atoms.append(atom)
    #        for n in atom.neighs:
    #            atoms.append(n)
    #        j += 1
    #    if(j > 20): break
    #atoms = list(set(atoms))
    #m = Model('',ag.model.lx,ag.model.ly,ag.model.lz,atoms)
    #for i,atom in enumerate(m.atoms):
    #    atom.id = i
    #    #print(atom.coord)
    #m.generate_neighbors(3.5)
    ##for atom in m.atoms:
    ##    if(len(atom.neighs) == 10):
    ##        print(atom.vp.index)
    ##        #print(atom.neighs)
    ##print(len(m.atoms))
    #m.write_cif()
    #return

    golden = True
    cluster_prefix = ''
    #cluster_prefix = 'ico.golden.t3.'
    #cluster_types = ['Icosahedra-like', 'Full-icosahedra'] # Need this for final/further analysis
    #clusters = ag.get_connected_clusters_with_neighs(cutoff, 'Icosahedra-like', 'Full-icosahedra')
    #clusters = ag.get_interpenetrating_clusters_with_neighs(cutoff, 'Icosahedra-like', 'Full-icosahedra')
    #clusters = ag.get_interpenetrating_atoms(cutoff, 'Icosahedra-like', 'Full-icosahedra')
    #v,e,f,i = ag.vefi_sharing('Icosahedra-like', 'Full-icosahedra')

    #cluster_prefix = 'fi.golden.t3.'
    #cluster_types = ['Full-icosahedra'] # Need this for final/further analysis
    #clusters = ag.get_connected_clusters_with_neighs(cutoff, 'Full-icosahedra')
    #clusters = ag.get_interpenetrating_clusters_with_neighs(cutoff, 'Full-icosahedra')
    #clusters = ag.get_interpenetrating_atoms(cutoff, 'Full-icosahedra')
    #v,e,f,i = ag.vefi_sharing('Full-icosahedra')

    #cluster_prefix = 'xtal.t3.'
    #cluster_types = ['Crystal-like'] # Need this for final/further analysis
    #clusters = ag.get_connected_clusters_with_neighs(cutoff, 'Crystal-like')
    #clusters = ag.get_interpenetrating_clusters_with_neighs(cutoff, 'Crystal-like')
    #clusters = ag.get_interpenetrating_atoms(cutoff, 'Crystal-like')
    #v,e,f,i = ag.vefi_sharing('Crystal-like')

    #cluster_prefix = 'mix.t3'
    #cluster_types = ['Mixed'] # Need this for final/further analysis
    #clusters = ag.get_connected_clusters_with_neighs(cutoff, 'Mixed')
    #clusters = ag.get_interpenetrating_clusters_with_neighs(cutoff, 'Mixed')
    #clusters = ag.get_interpenetrating_atoms(cutoff, 'Mixed')
    #v,e,f,i = ag.vefi_sharing('Mixed')

    cluster_prefix = 'undef.t3'
    cluster_types = ['Undef'] # Need this for final/further analysis
    #clusters = ag.get_connected_clusters_with_neighs(cutoff, 'Undef')
    clusters = ag.get_interpenetrating_clusters_with_neighs(cutoff, 'Undef')
    #clusters = ag.get_interpenetrating_atoms(cutoff, 'Undef')
    v,e,f,i = ag.vefi_sharing('Undef')

    print("V: {0}  E: {1}  F: {2}  I: {3}".format(v,e,f,i))



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
        # This changes all the atoms in the original cluster
        # (ie of type set above) to have no atomic number, 
        # which displays as Te (gold) in vesta for me
        if(golden):
            for atom in cluster:
                if(atom.vp.type in cluster_types):
                    atom.z = 0
        # Save cluster files
        cluster_model = Model("Orig cluster {0} contains {1} atoms.".format(i,len(cluster)),model.lx, model.ly, model.lz, cluster)
        #for atom in cluster_model.atoms:
        #    if(atom.vp.type in cluster_types): print('  {0}\t{1}'.format(atom,atom.vp.type))
        cluster_model.write_cif('{1}cluster{0}.cif'.format(i,cluster_prefix))
        cluster_model.write_our_xyz('{1}cluster{0}.xyz'.format(i,cluster_prefix))


    # Expand the clusters
    # This overwrites the old ones
    # Include all atoms already in the cluster and all atoms that have at least *5* neighbors in the cluster
    #for i,cluster in enumerate(clusters):
    #    new_clust = [ atom for atom in model.atoms if atom in clusters[i] or ( nnnic(atom,model,clusters[i]) >= 5 ) ]
    #    new_clust = list(set(new_clust))
    #    clusters[i] = new_clust[:]
    #
    #j = 0
    #for i,cluster in enumerate(clusters):
    #    print("New, bigger cluster {0} contains {1} atoms.".format(i,len(cluster)))
    #    #for atom in cluster:
    #    #    if atom not in orig_clusters[i]:
    #    #        #print("{0} {1}".format(atom.vesta(),atom.vp))
    #    #        print(atom.vesta())
    #    #        j += 1
    #    #    else:
    #    #        print(atom.vesta())
    #    # Save cluster files
    #    cluster_model = Model("New, bigger cluster {0} contains {1} atoms.".format(i,len(cluster)),model.lx, model.ly, model.lz, cluster)
    #    cluster_model.write_cif('{1}cluster{0}.big.cif'.format(i,cluster_prefix))
    #    cluster_model.write_our_xyz('{1}cluster{0}.big.xyz'.format(i,cluster_prefix))

    allclusters = []
    for cluster in clusters:
        for atom in cluster:
            if(atom not in allclusters):
                allclusters.append(atom)
                #if(atom.vp.type in cluster_types): print('  {0}\t{1}'.format(atom,atom.vp.type))
    allclusters = Model("All clusters.",model.lx, model.ly, model.lz, allclusters)
    allclusters.write_cif('{0}allclusters.cif'.format(cluster_prefix))
    allclusters.write_our_xyz('{0}allclusters.xyz'.format(cluster_prefix))
    print("{0}allclusters.cif and {0}allclusters.xyz contain {1} atoms.".format(cluster_prefix, allclusters.natoms))

    if(not golden):
        x_cluster = []
        for i,atom in enumerate(model.atoms):
            if atom not in allclusters.atoms:
                x_cluster.append(atom)
        x_cluster = Model("Opposite cluster of {0}".format(cluster_prefix),model.lx, model.ly, model.lz, x_cluster)
        x_cluster.write_cif('{0}opposite.cif'.format(cluster_prefix))
        x_cluster.write_our_xyz('{0}opposite.xyz'.format(cluster_prefix))
        print('{0}opposite.cif and {0}opposite.xyz contain {1} atoms.'.format(cluster_prefix, x_cluster.natoms))
    
    # Further analysis
    cn = 0.0
    for atom in model.atoms:
        cn += atom.cn
    cn /= model.natoms

    vpcn = 0.0
    count = 0
    for atom in ag.model.atoms:
        if( atom.vp.type in cluster_types ):
            vpcn += atom.cn
            count += 1
    vpcn /= count

    natomsinVPclusters = allclusters.natoms # Number of atoms in VP clusters
    nVPatoms = count # Number of VP atoms
    numsepVPatoms = nVPatoms * vpcn # Number of atoms in VP clusters if all clusters were separated
    maxnumatoms = model.natoms # Max number of atoms in VP clusters if all clusters were separated but still within the model size

    print('Average CN is {0}'.format(cn))
    print('Average CN of VP atoms is {0}'.format(vpcn))
    print('# atoms in all clusters: {0}. # VP atoms * vpcn: {1}. # VP atoms: {2}'.format(natomsinVPclusters,numsepVPatoms,nVPatoms))
    print('~ Number of VP that can fit in the model: {0}'.format(maxnumatoms/vpcn))
    print('Ratio of: (# atoms involved in VP clusters)/(# atoms involved in VP clusters if all clusters were completely separated):                          {0}%  <--- Therefore {1}% sharing.'.format(round(float(natomsinVPclusters)/(numsepVPatoms)*100.0,3),100.0-round(float(natomsinVPclusters)/(numsepVPatoms)*100.0,3)))
    print('Ratio of: (# atoms involved in VP clusters)/(# atoms involved in VP clusters if all clusters were separated as much as possible within the model): {0}%  <--- Therefore {1}% sharing.'.format(round(float(natomsinVPclusters)/min(numsepVPatoms,maxnumatoms)*100.0,3),100.0-round(float(natomsinVPclusters)/min(numsepVPatoms,maxnumatoms)*100.0,3) if numsepVPatoms < maxnumatoms else round(float(natomsinVPclusters)/min(numsepVPatoms,maxnumatoms)*100.0,3)))

    vor_instance = Vor()
    vor_instance.runall(modelfile,cutoff)
    index = vor_instance.get_indexes()
    vor_instance.set_atom_vp_indexes(model)
    vp_dict = categorize_vor.load_param_file('/home/jjmaldonis/model_analysis/scripts/categorize_parameters_iso.txt')
    atom_dict = categorize_vor.generate_atom_dict(index,vp_dict)
    categorize_vor.set_atom_vp_types(model,vp_dict)
    # Count the number of common neighbors in each of the VP
    vp_atoms = []
    for atom in model.atoms:
        if(atom.vp.type in cluster_types):
            vp_atoms.append(atom)
    common_neighs = 0.0
    atom_pairs = []
    for atomi in vp_atoms:
        for atomj in vp_atoms:
            if(atomi != atomj):
                if(atomi in atomj.neighs and [atomi,atomj] not in atom_pairs and [atomj,atomi] not in atom_pairs):
                    common_neighs += 1
                    atom_pairs.append([atomi,atomj])
                #if(atomj in atomi.neighs): common_neighs += 0.5
                for n in atomi.neighs:
                    if(n in atomj.neighs and [n,atomj] not in atom_pairs and [atomj,n] not in atom_pairs):
                        common_neighs += 1
                        atom_pairs.append([n,atomj])
                #for n in atomj.neighs:
                #    if(n in atomi.neighs): common_neighs += 0.5
    # Now common_neighs is the number of shared atoms
    #print(common_neighs)
    print('Percent shared based on common neighsbors: {0}'.format(100.0*common_neighs/natomsinVPclusters))
    
    # Debugging stuff, mainly for get_connected_clusters_with_neighs function
    #print(clusters[0][57], clusters[0][57].vesta())
    #for atom in ag.model.atoms:
    #cm = Model('temp',model.lx,model.ly,model.lz,clusters[0])
    #cm.generate_neighbors(3.5)
    #for i,atom in enumerate(cm.atoms):
    #    found = False
    #    for n in atom.neighs:
    #        if(n.vp.type in cluster_types):
    #            found = True
    #    if(not found and atom.vp.type not in cluster_types):
    #        #print('Found an atom that isnt connected to a VP type! {0} {1}'.format(allclusters.atoms.index(atom),atom.neighs))
    #        print('Found an atom that isnt connected to a VP type! {0} {1} {2} {3}'.format(i+1,atom,atom.neighs,atom.vp.type))
    #    #if(atom.neighs == []):
    #    #    for atom2 in ag.model.atoms:
    #    #        if(atom in atom2.neighs):
    #    #            print('Found an atom with empty neighbors as a neighbor of another! {0} {1} {2}'.format(atom,atom2,atom2.neighs))
    #    #if(clusters[0][57] in atom.neighs):
    #        #print(atom,atom.neighs)


if __name__ == "__main__":
    main()
