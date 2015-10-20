import sys, os, operator
import natsort, copy
from categorize_vor import voronoi_3d, index_stats
from model import Model
import numpy as np
from tools import save_obj, load_obj
import networkx as nx

def main():
    load_data_from_pickles = True
    if(not load_data_from_pickles): # If false, will skip to just load the matrix from file
        modelfiles = sys.argv[1:]
        modelfiles = [file for file in modelfiles if '.xyz' in file]
        modelfiles = natsort.natsorted(modelfiles)
        modelfiles = modelfiles[0:5000]
        new_modelfiles = [file for file in modelfiles if(not os.path.isfile(file[:-4]+'_VP.txt'))]

        cutoff = 3.6
        #print(modelfiles)
        print("Calculating VP for all models...")
        calculate_and_save_vp(new_modelfiles,cutoff)

        print("Finding all transformed VPs...")
        model_indexes = []
        for model in modelfiles:
            indexes = open(model[:-4]+'_VP.txt').readlines()
            indexes = [eval(ind.strip()) for ind in indexes]
            model_indexes.append(copy.copy(indexes))
        d_vpIndex_to_matIndex, d_matIndex_to_vpIndex, matrix = analyze(model_indexes)
        save_obj(d_vpIndex_to_matIndex, 'd_vpIndex_to_matIndex.pkl')
        save_obj(d_matIndex_to_vpIndex, 'd_matIndex_to_vpIndex.pkl')
        save_obj(matrix, 'matrix.pkl')
        f = open('output.txt','w')
        for row in matrix:
            row = ','.join([str(x) for x in row]) + '\n'
            f.write(row)
        f.close()

    # Post analysis
    if(load_data_from_pickles):
        d_vpIndex_to_matIndex = load_obj('d_vpIndex_to_matIndex.pkl')
        d_matIndex_to_vpIndex = load_obj('d_matIndex_to_vpIndex.pkl')
        matrix = load_obj('matrix.pkl')

    normalized_matrix = copy.copy(matrix)
    for i,row in enumerate(normalized_matrix):
        normalized_matrix[i][i] = 0
    row_total = [ 1.0*sum(row) for row in normalized_matrix ]
    #print(row_total)
    normalized_matrix = [ [x/row_total[i] for x in row] for i,row in enumerate(normalized_matrix) ]

    if(True):
        count = 0
        to_print = [[] for row in normalized_matrix]
        for i,row in enumerate(normalized_matrix):
            for j,x in enumerate(row):
                #if(i==j): continue
                if(x > .01):
                    line = "{0}:\t{1} -> {2}".format(round(100.0*x,3), d_matIndex_to_vpIndex[i], d_matIndex_to_vpIndex[j])
                    to_print[i].append(tuple([row_total[i],line]))
                    count += x*row_total[i]/100.0
        to_print = natsort.natsorted(to_print)
        to_print = [natsort.natsorted(row) for row in to_print]
        for row in to_print:
            for x,line in row:
                print(line + '\t' + str(100.0*x))
            print('')
        print("Total transformations: {0}".format(count))


    # Find shortest paths to (0, 0, 12, 0) and (0, 6, 0, 8)
    ico_index = (0, 0, 12, 0, 0, 0, 0, 0)
    bcc_index = (0, 6,  0, 8, 0, 0, 0, 0)
    graph = nx.Graph()
    for i,row in enumerate(normalized_matrix):
        for j,x in enumerate(row):
            if(i==j): continue
            if(x > 0.00):
                #graph.add_edge( d_matIndex_to_vpIndex[i], d_matIndex_to_vpIndex[j] )
                #graph[d_matIndex_to_vpIndex[i]][d_matIndex_to_vpIndex[j]]['weight'] = x
                graph.add_edge( d_matIndex_to_vpIndex[j], d_matIndex_to_vpIndex[i] )
                graph[d_matIndex_to_vpIndex[j]][d_matIndex_to_vpIndex[i]]['weight'] = 1-x
    #test = []
    bcc_dist = {}
    ico_dist = {}
    for ind in d_vpIndex_to_matIndex.keys():
        try:
            path = nx.shortest_path(graph, source=ind, target=ico_index, weight='weight')
            dist = 1.0
            for i in range(len(path)-1):
                dist = dist * (1-graph[path[i]][path[i+1]]['weight'])
            ico_dist[ind] = dist
        except nx.exception.NetworkXNoPath:
            ico_dist[ind] = 0.0
            #test.append(tuple([ dist*100,ind, len(path) ]))

        try:
            path = nx.shortest_path(graph, source=ind, target=bcc_index, weight='weight')
            dist = 1.0
            for i in range(len(path)-1):
                dist = dist * (1-graph[path[i]][path[i+1]]['weight'])
            bcc_dist[ind] = dist
        except nx.exception.NetworkXNoPath:
            bcc_dist[ind] = 0.0
    #test.sort()
    #for t in test:
    #    print(t)
    test = []
    for key in ico_dist:
        #print(key, ico_dist[key], bcc_dist[key], ico_dist[key]/bcc_dist[key], sum(matrix[d_vpIndex_to_matIndex[key]))
        test.append([key, ico_dist[key], bcc_dist[key], ico_dist[key]/bcc_dist[key], sum(matrix[d_vpIndex_to_matIndex[key]])])
    test.sort(key=operator.itemgetter(3))
    test.reverse()
    for t in test:
        t = [str(x) for x in t]
        t = '$'.join(t)
        print(t)
    


def calculate_and_save_vp(modelfiles,cutoff):
    #models = []
    for model in modelfiles:
        print(model)
        m = Model(model)
        voronoi_3d(m,cutoff)
        #models.append(m)
        f = open(model[:-4]+'_VP.txt','w')
        for atom in m.atoms:
            f.write(str(atom.vp.index)+'\n')
        f.close()

def analyze(model_indexes):
    all_indexes = set()
    for index_list in model_indexes:
        for index in index_list:
            all_indexes.add(index)
    all_indexes = list(all_indexes)
    numIndexes = len(all_indexes)
    print("Found {0} unique indexes.".format(numIndexes))

    d_vpIndex_to_matIndex = {}
    for i,index in enumerate(all_indexes):
        d_vpIndex_to_matIndex[tuple(index)] = i
    d_matIndex_to_vpIndex = dict((v, k) for k, v in d_vpIndex_to_matIndex.iteritems())

    matrix = np.zeros((numIndexes,numIndexes), dtype=np.int)
    count = 0
    for atomId in range(len(model_indexes[0])): # Each list in model_index contains the vp indexes for all the atoms. We want to go one atom at a time.
        #vps = [ model.atoms[atomId].vp.index for model in models ]
        vps_for_atom = [tuple(index_list[atomId]) for index_list in model_indexes]
        for i in range(0,len(vps_for_atom)-1):
            matrix[d_vpIndex_to_matIndex[vps_for_atom[i]]][d_vpIndex_to_matIndex[vps_for_atom[i+1]]] += 1
            matrix[d_vpIndex_to_matIndex[vps_for_atom[i+1]]][d_vpIndex_to_matIndex[vps_for_atom[i]]] += 1 # Symmetric
            
        #vps_set = set(vps)
        #if(len(vps_set) != 1):
        #    print(atomId, vps)
        #    count += 1
    #print(count)
    return d_vpIndex_to_matIndex, d_matIndex_to_vpIndex, matrix





if __name__ == '__main__':
    main()
