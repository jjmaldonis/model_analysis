import sys
import math
import numpy as np
import itertools
from model import Model
import znum2sym

def gr(m, nbins=0, types=None):
    if(m.lx != m.ly != m.lz): raise Exception("Can only do cubic boxes currently.")
    if(types == None):
        types = list(m.atomtypes) # These are ints (atomic numbers)

    count = 0
    for x in types:
        count += m.atomtypes[x]

    rho = m.natoms/(m.lx*m.ly*m.lz)
    if(nbins == 0):
        nbins = int(m.lx/rho)
    delg = m.lx*math.sqrt(3)/2.0/(nbins)
    #delg = 20.0/(nbins)

    g = np.zeros(nbins)
    for i in xrange(0,m.natoms):
        if(m.atoms[i].z in types):
            for j in xrange(i+1,m.natoms):
                if(m.atoms[j].z in types):
                    r = m.dist(m.atoms[i],m.atoms[j])
                    if(r < 20.0):
                        ig = int(r/delg)
                        g[ig] += 2 # contribution for particles i and j

    # Normalize gr
    for i in xrange(0,nbins):
        r = delg*(i+0.5)
        vb = ((i+1)**3-i**3)*(delg**3)
        nid = (4.0/3.0)*np.pi*vb*rho
        g[i] = g[i] / (nid*m.natoms)
        
    xaxis = np.arange(0.0,delg*nbins,delg)
    return (xaxis,g)

def generate_cutoffs(m):
    outcontent = []

    cutoff_types = []
    types = list(m.atomtypes)
    types.sort()
    for L in range(1,len(types)):
    #for L in range(1,len(types)+1):
        for subset in itertools.combinations(types, L):
            #print(subset)
            s = list(subset)
            xx,g = gr(m,types=s)
            g = list(g)
            s = [znum2sym.z2sym(y) for y in s]
            #g.insert(0,'gr_{0}'.format('_'.join(s)))
            outcontent.append(g)
            cutoff_types.append(s)
    cutoff_types = [ [x[0],x[0]] if len(x) == 1 else x for x in cutoff_types ]
    cutoff_types = [ [znum2sym.sym2z(y) for y in x] for x in cutoff_types ]
    cutoff_types = [ [types.index(y) for y in x] for x in cutoff_types ]
    #print(cutoff_types)

    cutoffs = np.zeros((len(m.atomtypes), len(m.atomtypes)))
    
    xx = list(xx)
    #xx.insert(0,'dr')
    outcontent.insert(0,xx)
    for i,arr in enumerate(outcontent[1:]):
        peak1 = arr.index(max(arr))
        #print('peak1',peak1,xx[peak1],arr[peak1])
        #print(arr[peak1:int(round(peak1*3.5))])
        min1 = peak1 + arr[peak1:int(round(peak1*3.5))].index(min(arr[peak1:int(round(peak1*3.5))]))
        #print('min1',min1,xx[min1],arr[min1])
        peak2 = min1 + arr[min1:].index(max(arr[min1:]))
        #print('peak2',peak2,xx[peak2],arr[peak2])
        #print('')
        cutoffs[cutoff_types[i][0]][cutoff_types[i][1]] = xx[min1]
        cutoffs[cutoff_types[i][1]][cutoff_types[i][0]] = xx[min1]
    #print(cutoffs)
    return cutoffs

def main():
    print("WARNING! The partial g(r)'s in this function might be a bit off (compared to RINGS calculations). Also, the total is a bit off too. However, if you add up all the partials you DO get the correct total (assuming RINGS is right).")
    modelfile = sys.argv[1]
    m = Model(modelfile)

    # This block of code prints the cutoffs
    cutoffs = generate_cutoffs(m)
    print(cutoffs)

    ## This block of code saves all the gr's to gr.txt
    #outcontent = []
    #types = list(m.atomtypes)
    #for L in range(1,len(types)+1):
    #    for subset in itertools.combinations(types, L):
    #        print(subset)
    #        s = list(subset)
    #        x,g = gr(m,types=s)
    #        g = list(g)
    #        s = [znum2sym.z2sym(y) for y in s]
    #        g.insert(0,'gr_{0}'.format('_'.join(s)))
    #        outcontent.append(g)
    #print("Writing output to gr.txt.")
    #x = list(x)
    #x.insert(0,'dr')
    #outcontent.insert(0,x)
    #outcontent = [ [str(y) for y in x] for x in zip(*outcontent)]
    #of = open('gr.txt','w')
    #for i in xrange(0,len(outcontent)):
    #    of.write(' '.join(outcontent[i])+'\n')

    ##outcontent = [ [ float(x) for x in l] for l in outcontent[1:]]
    #outcontent = [ [float(y) for y in x] for x in zip(*outcontent[1:])]
    ##print(outcontent)

    #x,g = gr(m)
    #x,g = gr(m,types=[13,40])
    #for i in xrange(0,len(g)):
    #    print("{0} {1}".format(x[i],g[i]))
    


if __name__ == "__main__":
    main()
