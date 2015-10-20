import sys, copy
from rot_3d import calc_rot_array, rot as rotate
from tools import printTable
from model import Model
from math import sqrt, sin
from scipy.optimize import curve_fit
import numpy as np


def rotate_and_save(prototype,alpha,beta,gamma):
    rot_arr = calc_rot_array(alpha,beta,gamma)
    rotated = copy.deepcopy(prototype)
    rotate(rotated,rot_arr)
    rotated.write_real_xyz(rotated.filename[:-4] + '_rotated_{0}.{1}.{2}'.format(alpha,beta,gamma)
        + rotated.filename[-4:])
    return rotated

def generate_arrays(alpha,beta,gamma,start=0,stop=360,step=15):
    a = alpha; b = beta; c = gamma;
    arrays = []
    array_code = []
    print('set rotations /1*{0}/;'.format((stop/step)**3))
    for a in range(start,stop,step):
        for b in range(start,stop,step):
            for c in range(start,stop,step):
                rot = calc_rot_array(a,b,c)
                arrays.append(rot)
                gamsrot = [ [0 for i in range(len(rot[0])+1)] for j in range(len(rot)+1)]
                for i in range(1,len(gamsrot)):
                    for j in range(1,len(gamsrot[i])):
                        gamsrot[i][j] = rot[i-1][j-1]
                for i in range(len(gamsrot)):
                    gamsrot[i][0] = i
                for j in range(len(gamsrot[0])):
                    gamsrot[0][j] = j
                gamsrot[0][0] = '_'
                ptable = printTable(gamsrot)
                ptable = ptable.replace('_',' ')
                array_code.append( ['* Angle = {0}'.format(a), 'table rotation{0}(d1,d2)'.format(a/step+1), ptable+';'] )
                print('* Angles = {0},{1},{2}'.format(a,b,c))
                print('table rotation{0}(d1,d2)'.format(len(array_code)))
                print(ptable+';')
    print('parameter RotationMats(rotations,d1,d2);')
    i = 1
    for a in range(start,stop,step):
        for b in range(start,stop,step):
            for c in range(start,stop,step):
                print("RotationMats('{0}',d1,d2) = rotation{0}(d1,d2);".format(i))
                i += 1
    return arrays
    
def main():
    prototype = Model(sys.argv[1])
    #polygon   = Model(sys.argv[2])
    alpha = 15
    beta = 45
    gamma = 135
    # Assume no permutation of the atoms
    polygon = rotate_and_save(prototype,alpha,beta,gamma)

    # Take the prototype and the polygon and rotate the polygon to the prototype
    # for every rotation matrix, and calculate the value of epsilon for that
    # rotation. Store the angle,epsilon values in two arrays for later analysis.
    start = 0
    stop = 360
    step = 15
    xdata = []
    ydata = []
    m = (100000,(1,2,3,),)
    #m = 5.9520869802
    count = 0
    for a in xrange(start,stop,step):
        for b in xrange(start,stop,step):
            for c in xrange(start,stop,step):
                rot = calc_rot_array(-a,-b,-c)
                rotpolygon   = copy.deepcopy(polygon)
                rotate(rotpolygon, rot)
                epsilon = sum( dist(prototype.atoms[i].coord,rotpolygon.atoms[i].coord) for i in range(prototype.natoms) )
                #if(epsilon < min[0]): m = (epsilon, (-a,-b,-c,),)
                #if(epsilon == m): print( count, (epsilon, (-a,-b,-c,),) )
                #print(round(epsilon,2))
                xdata.append(a)
                ydata.append(epsilon)
                count += 1
    print(min)

    # Do a fit to abs(sin(...)) based on the data recorded above
    xdata = np.array(xdata)
    ydata = np.array(ydata)
    popt, pcov = curve_fit(sin2, xdata, ydata,p0=[83.5,700,10,0])
    a,b,c,d = tuple(popt)
    #print(a,b,c,d)
    #print(pcov)
    for x,y in zip(xdata,ydata):
        print("{0}\t{1}".format(x,y))

    ## Reorder the data to make it all concave down
    #m = list(ydata).index(min(ydata))
    ##array_code = array_code[m:] + array_code[:m]
    #z = zip(xdata,ydata)
    #z = z[m:] + [(x + stop,y) for x,y in z[:m]]
    #xdata,ydata = zip(*z)
    #grad = np.gradient(ydata)
    #grad = np.gradient(grad) # The second derivative of the function should be negative everywhere
    #i = 0
    #for x,y in zip(xdata,ydata):
    #    print("{0}\t{1}\t{2}\t{3}".format(x,y,sin2(x,a,b,c,d),grad[i]))
    #    i += 1

    # Print the rotation code to paste into gams now that the rotations have been reordered
    # This is what gams is supposed to do itself though, this is all very costly.
    #for line in array_code:
    #    print('\n'.join(line))

def dist(x,y):
    return round( sqrt( (x[0]-y[0])**2 + (x[1]-y[1])**2 + (x[2]-y[2])**2 ), 10)

def sin2(x,a,b,c,d):
    #print(x,a,b,c)
    #return a * np.exp(-b * x) + c
    #return (a*np.sin(b*x+c))**2 + d
    return a*abs(np.sin(2*np.pi/b*x+c)) + d

def func(x, a, b, c):
    return a * np.exp(-b * x) + c


if __name__ == '__main__':
    main()
