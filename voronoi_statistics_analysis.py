from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import scipy, scipy.stats, scipy.integrate
from math import sqrt, log, exp, pi
import itertools

def normal_distribution(mu, sigma):
    return lambda x: 1.0/(sigma*sqrt(2*pi))*exp(-(x-mu)*(x-mu)/(2*sigma*sigma))

def drange(start, stop, step):
    r = start
    while r < stop:
        yield r
        r += step

binning = 500

mean1 = 504
mean2 = 594

ico = [555, 504, 594, 546, 625, 547, 592, 551, 540, 595]
mixed = [447, 474, 432, 474, 436, 463, 444, 485, 496, 485]

my_ico = [587, 542, 549]
my_xtal = [183, 190, 208]
my_mixed = [357, 381, 341]
my_natoms = 1250

jwh_ico = [544, 584, 600]
jwh_xtal = [193, 190, 194]
jwh_mixed = [424, 385, 373]
jwh_natoms = 1425

#my_ico = [0.4696, 0.4336, 0.4392]
#my_xtal = [0.1464, 0.152, 0.1664]
#my_mixed = [0.2856, 0.3048, 0.2728]

#jwh_ico = [0.3818, 0.4098, 0.4211]
#jwh_xtal = [0.1354, 0.1333, 0.1361]
#jwh_mixed = [0.2975, 0.2702, 0.2618]

ts = ['t1', 't2', 't3']
f = open('values.txt', 'w')
for i,mean1 in enumerate(my_xtal):
    for j,mean2 in enumerate(jwh_xtal):
        #if i>=j: continue
        xx1 = np.array(list(drange(-mean1/my_natoms/2,mean1/my_natoms/2,mean1/my_natoms/binning))[1:]) + mean1/my_natoms
        pdf1_func = normal_distribution(mean1/my_natoms, sqrt(mean1)/my_natoms)

        xx2 = np.array(list(drange(-mean2/jwh_natoms/2,mean2/jwh_natoms/2,mean2/jwh_natoms/binning))[1:]) + mean2/jwh_natoms
        pdf2_func = normal_distribution(mean2/jwh_natoms, sqrt(mean2)/jwh_natoms)

        newx = sorted(set(xx1).union(set(xx2)))
        minpdf = [min(pdf1_func(x), pdf2_func(x)) for k,x in enumerate(newx)]

        overlap = scipy.integrate.simps(minpdf, newx)

        print('my_'+ts[i], 'jwh_'+ts[j], mean1/my_natoms, mean2/jwh_natoms, round(overlap*100,1))
        f.write('{0}\n'.format(overlap))

        if True: # Plot the curves
            pdf1 = [pdf1_func(x) for x in xx1]
            pdf2 = [pdf2_func(x) for x in xx2]
            plt.plot(xx1,pdf1)
            plt.plot(xx2,pdf2)
            plt.plot(newx,minpdf)
            plt.show()
        break
    break
f.close()

