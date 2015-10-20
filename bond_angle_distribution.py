import sys
from model import Model
from pprint import pprint
from math import acos
import numpy as np

def bad(m,nbins):
    dtheta = np.pi/nbins
    divisor = 0
    g = np.zeros(nbins+1, np.int)
    for atomi in m.atoms:
        for j in range(0,atomi.cn):
            rij = m.dist(atomi,atomi.neighs[j])
            for k in range(j+1,atomi.cn):
                rik = m.dist(atomi,atomi.neighs[k])
                rjk = m.dist(atomi.neighs[j],atomi.neighs[k])
                temp = round( (rij**2 + rik**2 - rjk**2)/(2.0*rij*rik) , 10 )
                tijk = acos(temp)
                #try:
                #    tijk = acos(temp)
                #except ValueError:
                #    tijk = acos(round(temp))
                bin = int(tijk/dtheta)
                g[bin] += 1
        divisor += atomi.cn*(atomi.cn-1)
    #if(divisor > 0):
    #    g = [float(x)/float(divisor)*200 for x in g]

    dtheta = np.arange(0.0,180.0+dtheta*180.0/np.pi,dtheta*180.0/np.pi)
    return (dtheta,g)


def main():
    modelfile = sys.argv[1]
    m = Model(modelfile)
    m.generate_neighbors(3.50)

    dtheta = 0.5*np.pi/180.0
    nbins = np.pi/dtheta
    g = bad(m,nbins)
    for i in range(0,len(g[0])):
        print('{0}\t{1}'.format(g[0][i],g[1][i]))
    if(round(sum(g[1]),10) != 100.0): raise Exception("ERROR! You need to increase the neighbor distance. Only got {0}% of the distribution.".format(sum(g[1])))
    #print(sum(g[1]))


if __name__ == "__main__":
    main()
