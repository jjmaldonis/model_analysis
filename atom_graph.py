import copy
import warnings
from pprint import pprint
import os
import time
import random
import string
from vor import Vor
from voronoi_3d import voronoi_3d
import categorize_vor
from hutch import Hutch
from model import Model
import sys
import numpy as np

class AtomGraph(object):
    """ Implements a graph made up of atoms from a model.
        A cutoff is necessary to determine atom neighbors. """

    def __init__(self,modelfile,cutoff):
        """ Constructor
        @param cutoff is the cutoff we will use to determine neighbors """

        super(AtomGraph,self).__init__()
        
        self.model = Model(modelfile)
        #self.model.generate_neighbors(cutoff)

        #self.model.generate_coord_numbers()
        #print('Coordination numbers:')
        #pprint(self.model.coord_numbers)
        #self.model.print_bond_stats()

        # Generate CNs for different cutoffs. I can plot this and find
        # where it changes the least (ie deriv=0); this is a good spot
        # to set the cutoff distances because then the neighbors are
        # the least dependent on the cutoff distance.
        # I should make this into a function in model.py TODO
        #for cut in np.arange(2.0,4.6,0.1):
        #    self.model.generate_neighbors(cut)
        #    self.model.generate_coord_numbers()
        #    print("Cutoff: {0}".format(cut))
        #    for key in self.model.coord_numbers:
        #        if(len(key) < 4):
        #            print('  {0}: {1}'.format(key,self.model.coord_numbers[key]))

        #vor_instance = Vor()
        #vor_instance.runall(modelfile,cutoff)
        #index = vor_instance.get_indexes()
        #vor_instance.set_atom_vp_indexes(self.model)
        voronoi_3d(self.model,cutoff)
        
        #vorcats = VorCats('/home/jjmaldonis/OdieCode/vor/scripts/categorize_parameters.txt')
        self.vp_dict = categorize_vor.load_param_file('/home/jjmaldonis/model_analysis/scripts/categorize_parameters_iso.txt')
        #Eself.vp_dict = categorize_vor.load_param_file('/home/maldonis/model_analysis/scripts/categorize_parameters_iso.txt')
        #self.atom_dict = categorize_vor.generate_atom_dict(index,self.vp_dict)
        #vorcats.save(index)
        categorize_vor.set_atom_vp_types(self.model,self.vp_dict)

        #self.atom_dict = vorcats.get_atom_dict()
        ##for key in self.atom_dict:
        ##    print("{0} {1}".format(key,self.atom_dict[key]))
        #for key in self.atom_dict:
        #    for i in range(0,len(self.atom_dict[key])):
        #        self.atom_dict[key][i] = self.atom_dict[key][i][0]
        ##for key in self.atom_dict:
        ##    print("{0} {1}".format(key,len(self.atom_dict[key])))
        ##    #print("{0} {1}".format(key,self.atom_dict[key]))

        ##self.values = {}
        ##for atom in self.model.atoms:
        ##    for key in self.atom_dict:
        ##        if atom.id in self.atom_dict[key]:
        ##            self.values[atom] = key[:-1]
        ##    if atom not in self.values:
        ##        print("CAREFUL! SOMETHING WENT WRONG!")
        ##        #self.values[atom] = "Undef"
        ##for key in self.values:
        ##    print("{0}: {1}".format(key,self.values[key]))

        # Let's also run our VP algorithm to generate all that info.
        #voronoi_3d.voronoi_3d(self.model,cutoff)


    def get_neighs(self,atom):
        return atom.neighs

    def get_common_neighs(self,atom,*types):
        neighs = self.get_neighs(atom)
        list = [ atom for atom in neighs if atom.compute_vp_type(self.vp_dict) in types ]
        #list = [ atom for atom in neighs if atom.vp.type in types ]
        #list = []
        #for atom in neighs:
        #    if atom.get_vp_type(self.model.vp_dict) in types:
        #        list.append(atom)
        return list

    def get_unvisited_common_neighs(self,atom,visited,*types):
        comm_neighs = self.get_common_neighs(atom,*types)
        for atom in visited:
            if visited[atom] and atom in comm_neighs:
                comm_neighs.remove(atom)
        return comm_neighs

    def get_clusters(self,*cluster_types):
        """ This searches for interpenetrating clusters of cluster_types
        see http://arxiv.org/pdf/1302.1895.pdf """
        warnings.warn("Deprecated function, use one of the better ones!", DeprecationWarning)
        connections = {}
        #print("Connections:")
        # Generate a list for every atom that contains its neighbors of the type(s) we desire
        # This is what we use to generate our clusters
        for atom in self.model.atoms:
            # Does not include itself in get_common_neighs
            connections[atom] = self.get_common_neighs(atom,*cluster_types)
            #print("Atom {0}: {1}".format(atom,connections[atom]))

        clusters = []

        # Find an atom of type cluster_types
        atoms = self.model.atoms
        # Get all the atom INDEXES that have the same type as a type in cluster_types
        temp_atoms = [x[0] for ctype in cluster_types for x in self.atom_dict[ctype+':'] ]
        #print(len(temp_atoms))

        # This prints out a histogram of the number of interpenetrating connections
        # The x-axis printed out isn't quite right, you have to pick the integer that's
        # inside the bin and plot the Y axis vs that.
        lenconnections = [len(connections[x]) for x in connections if x.id in temp_atoms]
        #print(np.histogram(lenconnections,max(lenconnections)+1))

        # Go thru each atom I just found and append to a cluster all atoms ...
        for atom in temp_atoms:
            start_atom = atoms[atom] # atom represents the id in for this format - it's an int! not an atom
        
            visited = {start_atom:True}

            #neighs = self.get_unvisited_common_neighs(start_atom,visited,cluster_types)
            # neighs now contains all the neighbors of start_atom that have the same vp type as it
            # and that we have not visited already.

            paths = [] # this will contain all possible paths we find!
            path = [start_atom] # this will contain our current path
            #print("Path: {0}".format(path))
            queue = connections[start_atom]
            #print("Queue: {0}".format(queue))
            here = True
            prepop = -1
            while(len(path)):
                for atom in visited:
                    if atom in queue and visited[atom]:
                        queue.remove(atom)
                #print(len(queue))
                if(len(queue)): #if there is something in the queue
                    # move to the first atom in the queue
                    path.append(queue[0])
                    visited[queue[0]] = True
                    queue = connections[queue[0]]
                    here = True
                else: #if the queue is empty (ie we reached the end of a path)
                    if(here):
                        #print("Appending: {0}".format(path))
                        paths.append(path[:])
                    #else:
                    #    print("Not including: {0}".format(path))
                    here = False
                    # Go back to the last atom in 'path' and look thru the rest of its connections
                    if(type(prepop) != type(-1)):
                        visited[prepop] = False
                    #
                    prepop = path.pop()
                    if not len(path):
                        # if we have explored everything, stop
                        break
                    queue = connections[path[-1]]
            cluster = [ atom  for pathi in paths for atom in pathi ]
            cluster = list(set(cluster)) # remove duplicates
            cluster.sort()
            app = True
            if cluster != [] and cluster not in clusters:
                # Was having a problem before including clusters of size '1', fixed it now.
                if(len(cluster) == 1):
                    for c in clusters:
                        if(cluster[0] in c):
                            app = False
                if(app): clusters.append(cluster)
        #for cluster in clusters:
        #    print(cluster)
        #print(sum([len(l) for l in clusters]))
        return clusters


    def get_connected_clusters_with_neighs(self,cutoff,*cluster_types):
        """ Connected cluster finding via vertex sharing.
            Finds O -- O bonds and O -- X -- O bonds, where O
            represents an atom of the VP type(s) """
        # This code currently gives me a first nearest neighbor search. (Vertex sharing)
        m = Model(self.model.comment, self.model.lx, self.model.ly, self.model.lz, self.model.atoms[:])
        m.generate_neighbors(cutoff)
        count = 0
        for atom in m.atoms:
            keep = False
            if( atom.vp.type in cluster_types):
                keep = True
                #print('Keeping due to atom')
            ncount = 0
            if(not keep):
                temp = [n for n in atom.neighs if n.vp.type in cluster_types]
                if(len(temp) >= 1):
                    #keep = True
                    atom.neighs = [n for n in atom.neighs if n.vp.type in cluster_types]
            if(not keep):
                atom.neighs = [n for n in atom.neighs if n.vp.type in cluster_types]
                #print('Removing neighbors')
                #print(self.model.atoms[atom.id].neighs)
            else:
                count += 1
                if(atom.vp.type not in cluster_types): print(len(temp),ncount,atom)
        print('Total number of {0} atoms: {1}'.format(cluster_types,count))
        # Now I should be able to go through the graph/model's neighbors.
        clusters = []
        for atom in m.atoms:
            already_found = False
            for cluster in clusters:
                if atom in cluster:
                    already_found = True
            # If the VP atom isn't already in a cluster:
            if(not already_found and atom.vp.type in cluster_types):
                # Breadth first search
                queue = []
                visited = {}
                queue.append(atom)
                visited[atom] = True
                while( len(queue) ):
                    t = queue.pop()
                    for n in t.neighs:
                        if( not visited.get(n,False) ):
                            queue.append(m.atoms[m.atoms.index(n)])
                            visited[n] = True
                clusters.append(list(visited))
        for i,atom in enumerate(clusters[0]):
            found = atom.vp.type in cluster_types
            for n in atom.neighs:
                if(n.vp.type in cluster_types):
                    found = True
            if(not found):
                print('AG found an atom that isnt connected to a VP type! {0} {1} {2} {3}'.format(i+1,atom,atom.neighs,atom.vp.type))
                for atom2 in m.atoms:
                    if atom in atom2.neighs:
                        print('  It is connected to {0} {1} {2}'.format(atom2,atom2.neighs,atom2.vp.type))
                        for n in atom2.neighs:
                            print('   Dist from {0} to neighbor {1}: {2}. n.vp.type={3}'.format(atom2,n,m.dist(atom2,n),n.vp.type))
        for cluster in clusters:
            for atom in cluster:
                if(cluster.count(atom) > 1):
                    print('     ERROR!!!!')
                    #cluster.remove(atom)
        return clusters

    def get_interpenetrating_clusters_with_neighs(self,cutoff,*cluster_types):
        """ Interpenetrating cluster finding.
            O -- O bonds only, where O represents 
            an atom of the VP type(s) """
        clusters = self.get_interpenetrating_atoms(cutoff,*cluster_types)
        # Add neighbors on
        for cluster in clusters:
            neighs = []
            for atom in cluster:
                for n in self.model.atoms[self.model.atoms.index(atom)].neighs:
                    if n not in neighs and n not in cluster:
                        neighs.append(n)
            cluster += neighs
        for cluster in clusters:
            for atom in cluster:
                if(cluster.count(atom) > 1):
                    print('     ERROR!!!!')
                    #cluster.remove(atom)
        return clusters

    #def get_interpenetrating_atoms(self,cutoff,*cluster_types):
    #    """ Interpenetrating atom finding.
    #        O -- O bonds only, where O represents 
    #        an atom of the VP type(s).
    #        Neighbors are not included. """
    #    if(type(cluster_types[0]) == type(())):
    #        cluster_types = cluster_types[0]
    #    m = Model(self.model.comment, self.model.lx, self.model.ly, self.model.lz, self.model.atoms[:])
    #    count = 0
    #    for atom in m.atoms:
    #        if( atom.vp.type in cluster_types):
    #            atom.neighs = [ n for n in atom.neighs if n.vp.type in cluster_types ]
    #            count += 1
    #        else:
    #            atom.neighs = []
    #    print('Total number of {0} atoms: {1}'.format(cluster_types,count))
    #    # Now I should be able to go through the graph/model's neighbors.
    #    clusters = []
    #    for atom in m.atoms:
    #        already_found = False
    #        for cluster in clusters:
    #            if atom in cluster:
    #                already_found = True
    #        if( not already_found and atom.vp.type in cluster_types):
    #            # Breadth first search
    #            queue = []
    #            visited = {}
    #            queue.append(atom)
    #            visited[atom] = True
    #            while( len(queue) ):
    #                t = queue.pop()
    #                for n in t.neighs:
    #                    if( not visited.get(n,False) ):
    #                        queue.append(m.atoms[m.atoms.index(n)])
    #                        visited[n] = True
    #            clusters.append(list(visited))
    #    for cluster in clusters:
    #        for atom in cluster:
    #            if(cluster.count(atom) > 1):
    #                cluster.remove(atom)
    #    ## Add neighbors on too
    #    #for cluster in clusters:
    #    #    neighs = []
    #    #    for atom in cluster:
    #    #        for n in self.model.atoms[self.model.atoms.index(atom)].neighs:
    #    #            if n not in neighs and n not in cluster:
    #    #                neighs.append(n)
    #    #    cluster += neighs
    #    return clusters


    def vefi_sharing(self,*cluster_types):
        """ This function returns the number of vertex sharing,
        edge sharing, face sharing, and interpenetrating
        atoms in the clusters of atoms of type
        cluster_types + their nearest neighbors
        that are found in the model.
        The numbers are returned via a tuple in the order
        written in the function name: v, e, f, i. """
        # Vertex shared atoms will have will have 1 and
        # only 1 shared atom between the two VP.
        # Edge shared atoms will have 2 and only 2
        # shared atoms between the two VP.
        # Face shared atoms will have 3+ atoms shared
        # between the two VP.
        # Interpenetrating atoms will have the VPs as
        # neighbors of each other.
        # I need to go through each pair of VP and 
        # check for 1) interpenetrating, 2) face,
        # 3) edge, and then 4) vertex shared atoms.
        # When found, I should increment each counter.
        # Return the counters at the end in a tuple.
        if(type(cluster_types[0]) == type(())):
            cluster_types = cluster_types[0]
        vp_atoms = []
        for atom in self.model.atoms:
            if(atom.vp.type in cluster_types):
                vp_atoms.append(atom)
        vertex = 0.0
        edge = 0.0
        face = 0.0
        interpenetrating = 0.0
        atom_pairs = []
        for atomi in vp_atoms:
            for atomj in vp_atoms:
                common_neighs = 0.0
                if(atomi != atomj):
                    if(atomi in atomj.neighs and [atomi,atomj] not in atom_pairs and [atomj,atomi] not in atom_pairs):
                        interpenetrating += 1.0
                        atom_pairs.append([atomi,atomj])
                    for n in atomi.neighs:
                        #if(n in atomj.neighs and [n,atomj] not in atom_pairs and [atomj,n] not in atom_pairs):
                        if(n in atomj.neighs):
                            common_neighs += 1.0
                            #atom_pairs.append([n,atomj])
                    if(common_neighs == 1):
                        vertex += 0.5
                    elif(common_neighs == 2):
                        edge += 0.5
                    elif(common_neighs >= 3):
                        face += 0.5
        return(vertex,edge,face,interpenetrating)
    
    def get_interpenetrating_atoms(self,cutoff,*cluster_types):
        return self.get_sharing_clusters(cutoff,0,*cluster_types)
    def get_vertex_sharing_clusters(self,cutoff,*cluster_types):
        return self.get_sharing_clusters(cutoff,1,*cluster_types)
    def get_edge_sharing_clusters(self,cutoff,*cluster_types):
        return self.get_sharing_clusters(cutoff,2,*cluster_types)
    def get_face_sharing_clusters(self,cutoff,*cluster_types):
        return self.get_sharing_clusters(cutoff,3,*cluster_types)

    def get_sharing_clusters(self,cutoff,numneighs,*cluster_types):
        """ Connected cluster finding via vertex/edge/face sharing. 
            The last argument (1,2, or 3) specifies which. """
        if(type(cluster_types[0]) == type(())):
            cluster_types = cluster_types[0]
        if(numneighs == 0):
            temp = 'interepenetrating'
        elif(numneighs == 1):
            temp = 'interepenetrating and vertex'
        elif(numneighs ==2):
            temp = 'interepenetrating, vertex, and edge'
        elif(numneighs == 3):
            temp = 'interpenetrating, vertex, edge, and face'
        else:
            raise Exception("Wrong argument passsed to vertex/edge/face sharing cluster finding!")
        m = Model(self.model.comment, self.model.lx, self.model.ly, self.model.lz, self.model.atoms[:])
        m.generate_neighbors(cutoff)
        vp_atoms = []
        neighs = [[]]*m.natoms
        vp_atoms = [atom.copy() for atom in m.atoms if atom.vp.type in cluster_types]
        neighs = [[n for n in atom.neighs if n.vp.type in cluster_types] for atom in m.atoms]
        numfound = 0
        if(numneighs > 0): # Look for vertex, edge, or face sharing
            for i,atomi in enumerate(vp_atoms):
                print(i)
                # Interpenetrating
                ind = m.atoms.index(atomi)
                for j,atomj in enumerate(vp_atoms[vp_atoms.index(atomi)+1:]):
                    # Get all the neighbors they have in common
                    common_neighs = [n for n in atomi.neighs if n in atomj.neighs]
                    if(len(common_neighs) and (len(common_neighs) <= numneighs or numneighs == 3) ):
                        ind = m.atoms.index(atomi)
                        neighs[ind] = neighs[ind] + copy.copy([x for x in common_neighs if x not in neighs[ind]])
                        ind = m.atoms.index(atomj)
                        neighs[ind] = neighs[ind] + copy.copy([x for x in common_neighs if x not in neighs[ind]])
                        numfound += 1
        else:
            interpenetrating = sum(1 for atomi in vp_atoms for atomj in vp_atoms[vp_atoms.index(atomi)+1:] if atomi in atomj.neighs)
            numfound = interpenetrating
        for i,tf in enumerate(neighs):
            m.atoms[i].neighs = tf
        print('Total number of {0} atoms: {1}'.format(cluster_types,len(vp_atoms),temp))
        print('Total number of {2} sharing {0} atoms: {1}'.format(cluster_types,numfound,temp))
        # Now I should be able to go through the graph/model's neighbors.
        return self.search(m,cluster_types)

    def get_clusters_with_n_numneighs(self,cutoff,numneighs,cluster_types):
        m = Model(self.model.comment, self.model.lx, self.model.ly, self.model.lz, self.model.atoms[:])
        m.generate_neighbors(cutoff)
        vp_atoms = []
        #print(cluster_types)
        neighs = [[]]*m.natoms
        for i,atom in enumerate(m.atoms):
            if(atom.vp.type in cluster_types):
                vp_atoms.append(atom.copy())
        numfound = 0
        for i,atomi in enumerate(vp_atoms):
            for j,atomj in enumerate(vp_atoms[vp_atoms.index(atomi)+1:]):
                # Get all the neighbors they have in common
                #common_neighs = [n for n in atomi.neighs if n in atomj.neighs if n.vp.type not in cluster_types]
                common_neighs = [n for n in atomi.neighs if n in atomj.neighs]
                if(len(common_neighs) >= numneighs):
                    ind = m.atoms.index(atomi)
                    neighs[ind] = neighs[ind] + copy.copy([x for x in common_neighs if x not in neighs[ind]])
                    ind = m.atoms.index(atomj)
                    neighs[ind] = neighs[ind] + copy.copy([x for x in common_neighs if x not in neighs[ind]])
                    for n in common_neighs:
                        ind = m.atoms.index(n)
                        neighs[ind] = neighs[ind] + [x for x in [atomi,atomj] if x not in neighs[ind]]
                    numfound += 1
        for i,tf in enumerate(neighs):
            m.atoms[i].neighs = tf
        m.check_neighbors()
        print('Total number of {0} atoms: {1}'.format(cluster_types,len(vp_atoms)))
        print('Total number of {2}-sharing {0} atoms: {1}'.format(cluster_types,numfound,numneighs))
        # Now I should be able to go through the graph/model's neighbors.
        return self.search(m,cluster_types)


    def search(self,m,cluster_types):
        clusters = []
        for atom in m.atoms:
            #print(atom)
            already_found = False
            for cluster in clusters:
                if atom in cluster:
                    already_found = True
            # If the VP atom isn't already in a cluster:
            if(not already_found and atom.vp.type in cluster_types):
                # Breadth first search
                queue = []
                #visited = {}
                visited = [False]*m.natoms
                queue.append(m.atoms[m.atoms.index(atom)])
                #visited[m.atoms[m.atoms.index(atom)]] = True
                visited[atom.id] = True
                while( len(queue) ):
                    t = queue.pop()
                    #if(len(clusters) == 0):
                    #    print(list(visited))
                    #    print(t)
                    #    print(t.neighs)
                    #print(t.neighs)
                    for n in t.neighs:
                        #if( not visited.get(n,False) ):
                        if( not visited[n.id] ):
                            #print("going to", n)
                            queue.append(m.atoms[m.atoms.index(n)])
                            #visited[m.atoms[m.atoms.index(n)]] = True
                            visited[n.id] = True
                #clusters.append(list(visited))
                for i in visited:
                    if(i != m.atoms[i].id):
                        raise Exception("Calculating 'visited' this way will not work!")
                visited = [m.atoms[i] for i,x in enumerate(visited) if x]
                if(len(visited) > 1):
                    clusters.append(copy.copy(visited))
        return clusters


 


def main():
    model = Model(sys.argv[1])
    cluster_prefix = 'jason'
    #cluster_types = 'Crystal-like'
    #cluster_types = ['Icosahedra-like', 'Full-icosahedra']
    cluster_types = 'Icosahedra-like','Full-icosahedra'
    #cluster_types = 'Full-icosahedra'
    #cluster_types = 'Icosahedra-like'

    ag = AtomGraph(sys.argv[1],3.5)
    #clusters = ag.get_connected_clusters_with_neighs(3.5,cluster_types)
    #clusters = ag.get_connected_clusters_with_neighs(3.5,'Icosahedra-like','Full-icosahedra')
    #clusters = ag.get_interpenetrating_clusters_with_neighs(3.5,cluster_types)
    #clusters = ag.get_clusters(cluster_types)

    #clusters = ag.get_vertex_sharing_clusters(3.5,cluster_types)
    #clusters = ag.get_edge_sharing_clusters(3.5,cluster_types)
    #clusters = ag.get_face_sharing_clusters(3.5,cluster_types)
    #clusters = ag.get_interpenetrating_clusters_with_neighs(3.5,cluster_types)
    clusters = ag.get_interpenetrating_atoms(3.5,cluster_types)

    orig_clusters = clusters[:]
    # Print orig clusters
    j = 0
    for i,cluster in enumerate(clusters):
        print("Orig cluster {0} contains {1} atoms.".format(i,len(cluster)))
        # This changes all the atoms in the original cluster
        # (ie of type set above) to have no atomic number, 
        # which displays as Te (gold) in vesta for me
        #for atom in cluster:
        #    if(i==0): print(atom.vp.type)
        #    if(atom.vp.type in cluster_types):
        #        atom.z = 0
        # Save cluster files
        cluster_model = Model("Orig cluster {0} contains {1} atoms.".format(i,len(cluster)),model.lx, model.ly, model.lz, cluster)
        #for atom in cluster_model.atoms:
        #    if(atom.vp.type in cluster_types): print('  {0}\t{1}'.format(atom,atom.vp.type))
        cluster_model.write('{1}cluster{0}.cif'.format(i,cluster_prefix))
        cluster_model.write('{1}cluster{0}.xyz'.format(i,cluster_prefix))
        #for atom in cluster:
        #    print(atom)
   

if __name__ == "__main__":
    main()
        
