import sys
import numpy as np
from math import cos,sin,sqrt
from pprint import pprint
import scipy.spatial.distance

from fastModel import fastModel
from model import Model


def rot_normal_model(model,arr):
    """ arr should be a 9 element rotation numpy array, which we will reshape here """
    if(type(arr) == list): arr = np.array(arr)
    if(arr.shape == (9,)): arr = arr.reshape((3,3))
    arr = np.linalg.inv(arr)

    for i in range(0,model.natoms):
        old_coord = [model.atoms[i].coord[0], model.atoms[i].coord[1], model.atoms[i].coord[2]]
        new_coord = np.dot(np.asarray(old_coord), arr)
        model.atoms[i].coord = (new_coord[0], new_coord[1], 0.0)


def rot(model,arr):
    """ arr should be a 9 element rotation numpy array, which we will reshape here """
    if(type(arr) == list): arr = np.array(arr)
    if(arr.shape == (9,)): arr = arr.reshape((3,3))
    arr = np.linalg.inv(arr)

    for i in range(0,model.natoms):
        old_coord = [model.atoms[i][0], model.atoms[i][1], model.atoms[i][2]]
        new_coord = np.dot(np.asarray(old_coord), arr)
        model.atoms[i][0] = new_coord[0]
        model.atoms[i][1] = new_coord[1]
        model.atoms[i][2] = 0.0#new_coord[2]


def calc_rot_array(t1,t2,t3):
    """ We construct the rotation matrix based on t1,t2,t3
        NOTE! Order matters! """
    # http://mathworld.wolfram.com/RotationMatrix.html
    rx = np.array( [ [1,0,0], [0,cos(t1),sin(t1)], [0,-sin(t1),cos(t1)] ] )
    ry = np.array( [ [cos(t2),0,-sin(t2)], [0,1,0], [sin(t2),0,cos(t2)] ] )
    rz = np.array( [ [cos(t3),sin(t3),0], [-sin(t3),cos(t3),0], [0,0,1] ] )
    arr = np.dot( np.dot(ry,rx) ,rz)
    return np.linalg.inv(arr)


def dist_2d(atomi,atomj,m):
    x = (atomj[0] - atomi[0])
    y = (atomj[1] - atomi[1])
    while(x > m.lx/2): x = -m.lx + x
    while(y > m.ly/2): y = -m.ly + y
    while(x < -m.lx/2): x = m.lx + x
    while(y < -m.ly/2): y = m.ly + y
    return sqrt(x**2+y**2)

def dist_2d_arry(x0,x1,dimensions):
    delta = np.abs(x0-x1)
    delta = np.where(delta > 0.5*dimensions, dimensions-delta,delta)
    return np.sqrt((delta**2).sum(axis=-1))


def main():
    modelfile = sys.argv[1]
    m = Model(modelfile)
    #(93, 440, 80, 17, 0.1121997376282069, 2.9919930034188509, 0.63579851322650582)
    #ra = calc_rot_array(0.1121997376282069, 2.9919930034188509, 0.63579851322650582)
    #(87, 440, 80, 28, 0.1121997376282069, 2.9919930034188509, 1.0471975511965979)
    ra = calc_rot_array(0.1121997376282069, 2.9919930034188509, 1.0471975511965979)
    rot_normal_model(m,ra)
    m.write_cif('found_rot_3.cif')
    #m.write_our_xyz('found_rot_1.xyz')
    #(93, 440, 80, 59, 0.1121997376282069, 2.9919930034188509, 2.2065948400214026)
    #m = Model(modelfile)
    #ra = calc_rot_array(0.1121997376282069, 2.9919930034188509, 2.2065948400214026)
    #rot_normal_model(m,ra)
    #m.write_cif('found_rot_2.xyz')
    #m.write_our_xyz('found_rot_2.xyz')


    return
    modelfile = sys.argv[1]
    m = fastModel(modelfile)
    for i in range(0,m.natoms):
        m.atoms[i][2] = 0.0

    deltaDist = 0.3
    numangles = 84

    alpha = np.arange(0.0,np.pi,np.pi/numangles)
    beta  = np.arange(0.0,np.pi,np.pi/numangles)
    gamma = np.arange(0.0,np.pi,np.pi/numangles)

    numdone = 0
    results = []
    for i,a in enumerate(alpha):
        print("{0}% done".format(100*i/float(numangles)))
        for j,b in enumerate(beta):
            for k,g in enumerate(gamma):
                mrot = fastModel(modelfile)
                ra = calc_rot_array(a,b,g)
                rot(mrot,ra)
                #mrot.write_cif("{0}{1}{2}.cif".format(i,j,k))
                count = 0
                for i in range(0,mrot.natoms):
                    dists = dist_2d_arry(mrot.atoms[i+1:],[mrot.atoms[i]],np.array([mrot.lx,mrot.ly,mrot.lz]))
                    for d in dists:
                        if(d < deltaDist): count += 1
                results.append( (count,i,j,k,a,b,g) )
                numdone += 1
                print("finished a model... {0} {1}".format(numdone,count))

    results.sort(key=lambda x:x[0])
    pprint(results)



if __name__ == '__main__':
    main()
