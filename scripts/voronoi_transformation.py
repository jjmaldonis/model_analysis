import sys, os
import natsort, copy
from categorize_vor import voronoi_3d, index_stats
from model import Model
import numpy as np
from tools import save_obj

def main():
    modelfiles = sys.argv[1:]
    modelfiles = [file for file in modelfiles if '.xyz' in file]
    modelfiles = natsort.natsorted(modelfiles)
    #modelfiles = modelfiles[0:10]
    new_modelfiles = [file for file in modelfiles if(not os.path.isfile(file[:-4]+'_VP.txt'))]

    cutoff = 3.6
    #print(modelfiles)
    print("Calculating VP for all models...")
    calculate_and_save_vp(new_modelfiles,cutoff)

    print("Finding any transformed VP...")
    model_indexes = []
    for model in modelfiles:
        indexes = open(model[:-4]+'_VP.txt').readlines()
        indexes = [eval(ind.strip()) for ind in indexes]
        model_indexes.append(copy.copy(indexes))
    d_vpIndex_to_matIndex, d_matIndex_to_vpIndex, matrix = analyze(model_indexes)

    # Post analysis
    save_obj(d_vpIndex_to_matIndex, 'd_vpIndex_to_matIndex.pkl')
    save_obj(d_matIndex_to_vpIndex, 'd_matIndex_to_vpIndex.pkl')
    save_obj(matrix, 'matrix.pkl')
    f = open('output.txt','w')
    for row in matrix:
        row = ','.join([str(x) for x in row]) + '\n'
        f.write(row)
    f.close()
    count = 0
    for i,row in enumerate(matrix):
        for j,x in enumerate(row):
            if(i==j): continue
            if(x > 0):
                print("{0}:   {1} -> {2}".format(x, d_matIndex_to_vpIndex[i], d_matIndex_to_vpIndex[j]))
                count += x
    print("Total transformations: {0}".format(count))


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
            
        #vps_set = set(vps)
        #if(len(vps_set) != 1):
        #    print(atomId, vps)
        #    count += 1
    #print(count)
    return d_vpIndex_to_matIndex, d_matIndex_to_vpIndex, matrix





if __name__ == '__main__':
    main()
