import sys
import numpy as np
import os
from model import Model
import categorize_vor
import vor

def moransI(xx,yy):
    # Assume atom 1 in m1 is the same atom as atom 1 in m2, just after some number of moves.

    # Not what I wanted...

    # have to use spacial weights to do this, which doesnt make sense here at all

    meanx = np.mean(xx)
    meany = np.mean(yy)

    sum2rxi = 0.0
    numx = 0.0
    for i in range(0,len(xx)):
        rxi = xx[i] - meanx
        #print(rxi)
        sum2rxi = sum2rxi + (rxi - meanx)**2
        for j in range(0,len(yy)):
            #rxj = yy[j] - meany
            #print('   ',rxj)
            #numx = numx + rxi*rxj
            numx = numx + (xx[i] - meanx)*(yy[j] - meany)
    #print(numx,sum2rxi)
    numx = round(numx,7)
    sum2rxi = round(sum2rxi,7)
    print(numx,sum2rxi)
    if(numx == sum2rxi == 0): return 0
    numx = numx / sum2rxi

    return np.array(numx)

def main():
    paramfile = sys.argv[1]
    models_prefix = sys.argv[2]

    vp_dict = categorize_vor.load_param_file(paramfile)

    print("Collecting paths and models...")
    path = models_prefix[:models_prefix.rfind('/')+1]
    files = os.listdir(path)
    modelfiles = []
    for file in files:
        if models_prefix[models_prefix.rfind('/')+1:] in file:
            modelfiles.append(path+file)
    print("Sorting by filename...")
    modelfiles.sort()

    print('Last filename: {0}'.format(modelfiles[-1]))

    print("Running voronoi analysis on all models.")
    vorrun = vor.Vor()
    i = 1
    for model in modelfiles:
        m = Model(model)
        good = True
        try:
            vorrun.runall(model,3.5)
        except:
            good = False
            print("Model {0} failed!".format(model))
        if(good):
            print("Step: {0}; Model: {1}".format(i,model))
            vorrun.set_atom_vp_indexes(m)
            categorize_vor.set_atom_vp_types(m,vp_dict)
            categorize_vor.vor_stats(m) # Prints what you probably want
        i += 1



if __name__ == "__main__":
    main()
