import sys
import copy
from model import Model
from vor import fortran_voronoi_3d
from voronoi_3d import voronoi_3d
#from categorize_vor import categorize_index

def drange(start, stop, step):
    r = start
    while r < stop:
        yield r
        r += step

def compare_models_vp_indexes(m1,m2):
    """ Compares each atom.id's atom.vp.index in m1 and m2
    and return a dictionary of the index changes. """
    dict = {}
    for atom in m1.atoms:
        dict[(tuple(m1.atoms[atom.id].vp.index),tuple(m2.atoms[atom.id].vp.index))] = dict.get((tuple(m1.atoms[atom.id].vp.index),tuple(m2.atoms[atom.id].vp.index)),0) + 1
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
            if(deltacut == 0):
                continue
            num += 1
            print("Running {0} of {1}".format(num,(cuttol/cutdelta*2) * len(keys)+1))
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

    print(models)
    dict = {}
    for i,modeli in enumerate(models[:-1]):
        for modelj in models[i+1:]:
            td = compare_models_vp_indexes(modeli,modelj)
            for key in td:
                dict[key] = dict.get(key,0) + td[key]
    
    double = True
    while(double == True):
        for key in dict:
            if(key[::-1] in dict):
                dict[key] += dict[key[::-1]]
                del dict[key[::-1]]
                break
        else:
            double = False

    l = [(dict[key],key) for key in dict]
    l.sort()
    for item in l:
        if(item[1][0] != item[1][1]):
            print("{0}:\t{1}".format(item[0],item[1]))


if __name__ == "__main__":
    main()
