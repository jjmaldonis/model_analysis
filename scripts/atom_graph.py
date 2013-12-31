import os
import time
import random
import string
import gen_paramfile
from vor import Vor
from categorize_vor import VorCats
from hutch import Hutch
import sys

class AtomGraph(object):
    """ Implements a graph made up of atoms from a model.
        A cutoff is necessary to determine atom neighbors. """

    def __init__(self,modelfile,cutoff):
        """ Constructor
        @param hd is the hutch dictionary implemented in hutch.py
        @param cutoff is the cutoff we will use to determine neighbors """

        super(AtomGraph,self).__init__()
        
        hd = Hutch(modelfile)
        self.atoms = hd.get_all_atoms()

        self.neighs = {}
        for atom in self.atoms:
            self.neighs[atom] = hd.get_atoms_in_cutoff(atom,cutoff)
            for i in range(0,len(self.neighs[atom])):
                # Neighbor syntax: (neighbor atom, dist to main atom)
                self.neighs[atom][i] = (self.neighs[atom][i],hd.dist(atom,self.neighs[atom][i]))

        vor_instance = Vor()
        vor_instance.runall(cutoff,modelfile)
        index = vor_instance.get_index()
        
        vorcats = VorCats('categorize_parameters.txt')
        vorcats.save(index)
        
        self.vp_dict = vorcats.get_vp_dict()

        sorted_atoms = vorcats.get_atom_dict()
        self.atom_dict = sorted_atoms
        for key in self.atom_dict:
            for i in range(0,len(self.atom_dict[key])):
                self.atom_dict[key][i] = self.atom_dict[key][i][0]

        self.values = {}
        for atom in self.atoms:
            for key in self.atom_dict:
                if atom.id in self.atom_dict[key]:
                    self.values[atom] = key[:-1]
            if atom not in self.values:
                self.values[atom] = "Undef"

    def get_vp_type(self,atom):
        return self.values[atom]
    def get_neighs(self,atom):
        return self.neighs[atom]

    def get_common_neighs(self,atom,type):
        neighs = self.neighs[atom]
        list = []
        for atom,dist in neighs:
            if self.get_vp_type(atom) == type:
                list.append(atom)
        return list

    def get_unvisited_common_neighs(self,atom,type,visited):
        comm_neighs = self.get_common_neighs(atom,type)
        for atom in visited:
            if visited[atom] and atom in comm_neighs:
                comm_neighs.remove(atom)
        return comm_neighs

    def get_all_atoms(self):
        return self.atoms


def main():
    cluster_type = 'Crystal-like'
    #cluster_type = 'Icosahedra-like'
    ag = AtomGraph(sys.argv[1],3.5)

    connections = {}
    #print("Connections:")
    for atom in ag.get_all_atoms():
        connections[atom] = ag.get_common_neighs(atom,cluster_type)
        #print("Atom {0}: {1}".format(atom,connections[atom]))

    clusters = []

    # Find an atom of type cluster_type
    atoms = ag.get_all_atoms()
    for atom in ag.atom_dict[cluster_type+':']:
        start_atom = atoms[atom] # atom represents the id in for this format - it's an int! not an atom
    
        visited = {start_atom:True}

        neighs = ag.get_unvisited_common_neighs(start_atom,cluster_type,visited)
        # neighs now contains all the neighbors of start_atom that have the same vp type as it
        # and that we have not visited already.

        paths = [] # this will contain all possible paths we find!
        path = [start_atom] # this will contain our current path
        queue = connections[start_atom]
        here = False
        prepop = -1
        while(len(path)):
            for atom in visited:
                if atom in queue and visited[atom]:
                    queue.remove(atom)
            if(len(queue)):
                path.append(queue[0])
                visited[queue[0]] = True
                queue = connections[queue[0]]
                here = True
            else:
                if(here):
                    #print(path)
                    paths.append(path[:])
                #else:
                #    print("Not including: {0}".format(path))
                here = False
                if(prepop != -1):
                    visited[prepop] = False
                #
                prepop = path.pop()
                if not len(path):
                    break
                queue = connections[path[-1]]
        cluster = []
        for pathi in paths:
            #print(pathi)
            for atom in pathi:
                if atom not in cluster:
                    cluster.append(atom)
        cluster.sort()
        if cluster != [] and cluster not in clusters:
            clusters.append(cluster)
    
    i = 1
    for cluster in clusters:
        print("Unique cluster {0}:".format(i))
        i += 1
        for atom in cluster:
            print(atom.convert_to_sym())

        



    

if __name__ == "__main__":
    main()
        
