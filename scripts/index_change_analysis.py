import sys
import copy
from model import Model
from vor import fortran_voronoi_3d
from voronoi_3d import voronoi_3d
#from categorize_vor import categorize_index
import networkx as nx

def path_length(path,d):
    length = 0.0
    for i,elem in enumerate(path):
        if(i == 0):
            continue
        length += d[(path[i-1],path[i])]
    return length

def create_graph_from_dict(d):
    G = nx.DiGraph()
    for indpair,val in d.iteritems():
        if(val != 0.0):
            #G.add_edge(indpair[0],indpair[1],{'weight':1.0/val})
            G.add_edge(indpair[0],indpair[1],{'weight':1.0-val})
    return G


def create_matrix_from_dict(d):
    import numpy as np
    keys = d.keys()
    indexes = []
    for key in keys:
        if(key[0] not in indexes):
            indexes.append(key[0])
        if(key[1] not in indexes):
            indexes.append(key[1])
    nelems = len(indexes)

    print("LENGTH OF INDEXES: {0}".format(nelems))
    # Create dictionary for fast index index lookup
    lookup = {}
    for i,ind in enumerate(indexes):
        lookup[ind] = i

    mat = np.zeros((nelems,nelems),dtype=np.float)
    for indpair in keys:
        mat[ lookup[indpair[0]] ][ lookup[indpair[1]] ] = d[indpair]

    print(lookup)
    #print_matrix(mat)
    #print_row_of_matrix_as_column(mat,lookup[(0, 3, 6, 4, 0, 0, 0, 0)])
    #print_row_of_matrix_as_column(mat,lookup[(0, 2, 8, 2, 0, 0, 0, 0)])
    #print_transitions(mat,lookup,(0, 2, 8, 2, 0, 0, 0, 0))
    print_all_transitions(mat,lookup)
    return lookup,mat

def print_matrix(mat):
    print('\t|\t'.join([str(x) for x in range(0,len(mat))]))
    print('')
    for i,line in enumerate(mat):
        print(str(i) + '\t' + '\t|\t'.join([str(round(x*100,1)) for x in line]))

#def print_column_of_matrix(mat):

def print_row_of_matrix_as_column(mat,row):
    for i,x in enumerate(mat[row]):
        print(str(i) + '\t' + str(round(x*100,1)))

def print_transitions(mat,lookup,index):
    print("Transitions of {0}:".format(index))
    row = lookup[index]
    for i,x in enumerate(mat[row]):
        if(round(x*100,2) >= 1):
            key = (key for key,value in lookup.items() if value==i).next()
            print('  --> ' + str(key) + ':  \t' + str(round(x*100,2))+'%')

def print_all_transitions(mat,lookup):
    for index,row in lookup.items():
        print("Transitions of {0}:".format(index))
        for i,x in enumerate(mat[row]):
            if(round(x*100,2) >= 1):
                key = (key for key,value in lookup.items() if value==i).next()
                print('  --> ' + str(key) + ':  \t' + str(round(x*100,2))+'%')

def sum_dicts(dict_list):
    dict = {}
    for d in dict_list:
        for key in d:
            dict[key] = dict.get(key,0) + d[key]
    return dict

def drange(start, stop, step):
    r = start
    while r < stop:
        yield r
        r += step

def compare_models_vp_indexes(m1,m2):
    """ Compares each atom.id's atom.vp.index in m1 and m2
    and return a dictionary of the index changes.
    Each key in the dictionary is a 2-tuple of VP indexes,
    the first of which is the original index and the second
    is the index the original changed to. They value of this
    2-tuple key is the number of times that change occured. """
    td = count_number_of_indexes_in_model(m1)
    dict = {}
    for atom in m1.atoms:
        dict[(tuple(m1.atoms[atom.id].vp.index),tuple(m2.atoms[atom.id].vp.index))] = dict.get((tuple(m1.atoms[atom.id].vp.index),tuple(m2.atoms[atom.id].vp.index)),0) + 1.0#/td[tuple(m1.atoms[atom.id].vp.index)]
    return dict

def count_number_of_indexes_in_model(m):
    dict = {}
    for atom in m.atoms:
        dict[tuple(atom.vp.index)] = dict.get(tuple(atom.vp.index),0.0) + 1.0
    return dict

def main():
    modelfile = sys.argv[1]
    keys = [(40,40),(29,29),(13,13),(40,29),(40,13),(29,13)]
    cuttol = 0.2
    cutdelta = 0.1
    cutoff = {}
    cutoff[(40,40)] = 3.8
    cutoff[(13,29)] = 3.5
    cutoff[(29,13)] = 3.5
    cutoff[(40,13)] = 3.7
    cutoff[(13,40)] = 3.7
    cutoff[(29,40)] = 3.5
    cutoff[(40,29)] = 3.5
    cutoff[(13,13)] = 3.5
    cutoff[(29,29)] = 3.5

    models = []
    num = 1
    #print("Number of runs: {0}".format((cuttol/cutdelta*2) * len(keys)+1))
    #models.append(fortran_voronoi_3d(modelfile,cutoff))
    print("Running {0} of {1}".format(num,(cuttol/cutdelta*2) * len(keys)+1))
    m = Model(modelfile)
    voronoi_3d(m,cutoff)
    models.append(m)
    for key in keys:
        for deltacut in drange(-cuttol,cuttol+cutdelta,cutdelta):
        #for deltacut in drange(0,cuttol+cutdelta,cutdelta):
            if(deltacut == 0):
                continue
            num += 1
            print("Running {0} of {1}".format(num,(cuttol/cutdelta*2) * len(keys)+1))
            #if(num > 2): break
            #print(key,deltacut)
            cutoff2 = copy.deepcopy(cutoff)
            cutoff2[key] += deltacut
            if(key != key[::-1]):
                cutoff2[key[::-1]] += deltacut
            #print(cutoff2)
            #models.append(fortran_voronoi_3d(modelfile,cutoff2))
            m = Model(modelfile)
            voronoi_3d(m,cutoff2)
            models.append(m)

    print("Got {0} models.".format(len(models)))
    index_in_models = [count_number_of_indexes_in_model(m) for m in models]
    #for item,key in index_in_models[1].items():
    #    print("{0}: {1}".format(item,key))
    models_containing_index = {}
    for d in index_in_models: # Go thru each model
        for index in d: # if the index was in the model, increment
            models_containing_index[index] = models_containing_index.get(index,0) + 1
    #index_in_models = sum_dicts(index_in_models)
    dict = {}
    dict2 = {}
    for i,modeli in enumerate(models):
        #for j,modelj in enumerate(models[i+1:]):
        #for j,modelj in enumerate(models[1:]):
        for j,modelj in enumerate(models):
            #realj = i+1+j
            if(i != j):
                td = compare_models_vp_indexes(modeli,modelj)
                for key in td:
                    #print(key, td[key], index_in_models[i][key[0]])
                    #if(key[0] == (1, 2, 6, 5, 1, 0, 0, 0)):
                        #print("The transition {0} to {1} took place {2} times. There were {3} {0} indexes in the model and {4} {1} indexes.".format(key[0],key[1],td[key],index_in_models[i][key[0]], index_in_models[i].get(key[1],0)))
                        #print("Added {1} to index {0}".format(key, td[key]/index_in_models[i][key[0]]))
                        #print(dict.get(key,0.0),td[key]/index_in_models[i][key[0]])
                    try:
                        #dict[key] = dict.get(key,0.0) + td[key]/(index_in_models[i][key[0]]-td.get((key[0],key[0]),0))
                        dict[key] = dict.get(key,0.0) + td[key]/index_in_models[i][key[0]]
                    except ZeroDivisionError:
                        dict[key] = dict.get(key,0.0) + 0
                    dict2[key] = dict2.get(key,0.0) + td[key]
    
    #for key in index_in_models:
    #    print("{0}: {1}".format(index_in_models[key],key))
    #for key in dict:
        #print(key)
        #print(key[0])
        #dict[key] /= float(index_in_models[key[0]])

    for key in dict:
        dict[key] /= (models_containing_index[key[0]]*(len(models)-1))

    transitions_of_type = {}
    for key in dict2:
        transitions_of_type[key[0]] = transitions_of_type.get(key[0],0) + dict2[key]

    #double = True
    #while(double == True):
    #    for key in dict:
    #        if(key[::-1] in dict and key != key[::-1]):
    #            dict[key] += dict[key[::-1]]
    #            del dict[key[::-1]]
    #            break
    #    else:
    #        double = False

    dict3 = {}
    for key in dict:
        #print("There were {0} {1} -> {2} transitions which composes {3}% of the non-identity {1} transitions.".format(dict2[key]-dict2.get((key[0],key[0]),0),key[0],key[1],round(100.0*dict2[key]/indexes_in_dict[key[0]],1)))
        if(key[0] != key[1]):
            dict3[key] = dict2[key]/( transitions_of_type[key[0]]-dict2.get((key[0],key[0]),0))
            #if( 100.0*dict2[key]/( transitions_of_type[key[0]]-dict2.get((key[0],key[0]),0)) > 10 and dict2[key] > 50):
                #print("There were {0} {1} -> {2} transitions which composes {3}% of the non-identity {1} transitions.".format(dict2[key],key[0],key[1], 100.0*dict2[key]/( transitions_of_type[key[0]]-dict2.get((key[0],key[0]),0) )))

                #print("{3}%: {1} -> {2} totaled {0}".format(dict2[key],key[0],key[1], round(100.0*dict2[key]/( transitions_of_type[key[0]]-dict2.get((key[0],key[0]),0)),2) ))

        else:
            #dict3[key] = dict[key]
            dict3[key] = 0
            #print("There were {0} {1} -> {2} transitions which composes {3}% of the {1} transitions.".format(dict2[key],key[0],key[1],round(dict[key]*100.0)))

    lookup,mat = create_matrix_from_dict(dict3)
    G = create_graph_from_dict(dict3)

    keys = dict3.keys()
    indexes = []
    for key in keys:
        if(key[0] not in indexes):
            indexes.append(key[0])
        if(key[1] not in indexes):
            indexes.append(key[1])
    nelems = len(indexes)
    for ind in indexes:
        try:
            path = nx.dijkstra_path(G,ind,(0, 0, 12, 0, 0, 0, 0, 0))
            print("Shortest path from {0} to FI:".format(ind))
            print(path)
            print("Length = {0}".format(path_length(path,dict3)))
        except:
            pass
        try:
            path = nx.dijkstra_path(G,(0, 0, 12, 0, 0, 0, 0, 0),ind)
            print("Shortest path from FI to {0}:".format(ind))
            print(path)
            print("Length = {0}".format(path_length(path,dict3)))
        except:
            pass
        try:
            path = nx.dijkstra_path(G,ind,(0, 6, 0, 8, 0, 0, 0, 0))
            print("Shortest path from {0} to BCC:".format(ind))
            print(path)
            print("Length = {0}".format(path_length(path,dict3)))
        except:
            pass
        try:
            path = nx.dijkstra_path(G,(0, 6, 0, 8, 0, 0, 0, 0),ind)
            print("Shortest path from BCC to {0}:".format(ind))
            print(path)
            print("Length = {0}".format(path_length(path,dict3)))
        except:
            pass
        print('')

    print("Shortest path from BCC to FI:")
    path = nx.dijkstra_path(G,(0, 6, 0, 8, 0, 0, 0, 0), (0, 0, 12, 0, 0, 0, 0, 0))
    print(path)
    print("Length = {0}".format(path_length(path,dict3)))

    print("Shortest path from HCP to FI:")
    print(path)
    path = nx.dijkstra_path(G,(0, 6, 0, 2, 0, 0, 0, 0), (0, 0, 12, 0, 0, 0, 0, 0))
    print("Length = {0}".format(path_length(path,dict3)))

    #l = [(dict[key],key) for key in dict]
    #l.sort()
    #for item in l:
    #    #if(item[1][0] != item[1][1]):
    #    if(dict2[item[1]] > (len(models)-1)*2):
    #    #and item[1][0] != item[1][1]):
    #        #print("{0}%:\t\t{1}. There were {2} total.".format(int(item[0]*100),item[1],dict2[item[1]]))
    #        print("{0}%:\t\t{1}.".format(int(item[0]*100),item[1]))
    #for item in l:
    #    sum = 0.0
    #    for item2 in l:
    #        if(item[1][0] == item2[1][0]):
    #            sum += item2[0]
    #    print("Sum of {0} = {1}     expected to be {2}".format(item[1][0],sum,models_containing_index[item[1][0]]))
    #    print("Sum of {0} = {1} correct??? expected to be {2}".format(item[1][0],sum-dict[(item[1][0],item[1][0])],models_containing_index[item[1][0]]))
        


if __name__ == "__main__":
    main()
