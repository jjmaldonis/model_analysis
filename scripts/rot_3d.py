import sys
import numpy as np
from model import Model
from math import cos, sin, atan2, sqrt

def rotate(m, a, axis):
    # This is rot from http://www.ks.uiuc.edu/Research/vmd/doxygen/Matrix4_8C-source.html
    # a is the angle in degrees
    mat = [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]
    angle = a*np.pi/180.0 # in radians
    if(axis == 'x'):
        mat = [[1,0,0],[0,cos(angle),sin(angle)],[0,-sin(angle),cos(angle)]]
    elif(axis == 'y'):
        mat = [[cos(angle),0,-sin(angle)],[0,1,0],[sin(angle),0,cos(angle)]]
    elif(axis == 'z'):
        mat = [[cos(angle),sin(angle),0],[-sin(angle),cos(angle),0],[0,0,1]]
    m = np.array(m)
    mat = np.array(mat)
    #print(np.dot(m,mat), np.dot(mat,m))
    return np.dot(mat,m)

def calc_rot_array_from_hkl(h,k,l):
    # This is transvecinv from http://www.ks.uiuc.edu/Research/vmd/doxygen/Measure_8C-source.html
    mat = [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]
    x = h
    y = k
    z = l
    if(x == 0.0 and y == 0.0):
        if(z == 0.0):
            return -1
        if( z > 0 ):
            mat = rotate(mat,-90,'y')
        else:
            mat = rotate(mat,90,'y')
        return 0
    theta = atan2(y,x)
    length = sqrt(x*x + y*y)
    phi = atan2(float(z), length)
    #mat = rotate(mat,phi*180.0/np.pi,'y')
    #mat = rotate(mat,-theta*180.0/np.pi,'z')
    mat = rotate(mat,theta*180.0/np.pi,'z')
    mat = rotate(mat,-phi*180.0/np.pi,'y')
    return mat


def rot(model, arr):
    """ arr should be a 9 element rotation numpy array, which we will reshape here
        NOTE! This changes the model atom positions in the model """
    if(type(arr) == list): arr = np.array(arr)
    if(arr.shape == (9,)): arr = arr.reshape((3,3))
    #print(arr)
    arr = np.linalg.inv(arr)
    
    for i,atom in enumerate(model.atoms):
        old_coord = [atom.coord[0], atom.coord[1], atom.coord[2]]
        new_coord = np.dot(np.asarray(old_coord), arr)
        atom.set_coord(new_coord[0],new_coord[1],new_coord[2])


def calc_rot_array(t1,t2,t3,deg=True):
    """ We construct the rotation matrix based on t1,t2,t3
        NOTE! Order matters! """
    if(deg):
        t1 = t1*np.pi/180.0 # in radians
        t2 = t2*np.pi/180.0 # in radians
        t3 = t3*np.pi/180.0 # in radians
    rx = np.array( [ [1,0,0], [0,cos(t1),-sin(t1)], [0,sin(t1),cos(t1)] ] )
    ry = np.linalg.inv(np.array( [ [cos(t2),0,sin(t2)], [0,1,0], [-sin(t2),0,cos(t2)] ] ))
    rz = np.linalg.inv(np.array( [ [cos(t3),-sin(t3),0], [sin(t3),cos(t3),0], [0,0,1] ] ))
    arr = np.dot( np.dot(ry,rx) ,rz)
    return np.linalg.inv(arr)


def main():
    #print(rotate([[2,0,0],[0,1,0],[0,0,1]], 90, 'y'))
    #print(calc_rot_array_from_hkl(19,-24,28))
    modelfile = sys.argv[1]
    m = Model(modelfile)
    rot_arr = calc_rot_array_from_hkl(41,60,-6)
    rot(m,rot_arr)
    m.write_real_xyz('temp.real.xyz')
    return

    # Below is a (the?) rotation matrix of Pei's t1 that gives some planes. Oriented for a specific plane ~.
    #rot_arr = [ -0.977103, -0.123352, -0.173361, -0.130450, 0.990997, 0.030118, 0.168085, 0.052043, -0.984398 ]
    #rot(m,rot_arr)

    # Angles in radians
    # Note that these are semi difficult to figure out from the vesta rotation matrix,
    # partly because there are negative angles, so you may need to do 2pi - angle you found.
    #t1 = np.pi*2 - 0.0371505
    #t2 = 0.162790
    #t3 = 0
    #rot_arr = calc_rot_array(m,t1,t2,t3)
    #rot(m,rot_arr)

    kx = -0.76094085
    ky = 0.028182994
    kz = -0.648208872
    t2 = np.arctan(-ky/kx)
    t3 = 0.0
    t1 = np.arctan( (kx*np.cos(t2)-ky*cos(t2))/kz )
    t1 = 0.0
    print(t1,t2,t3)
    rot_arr = calc_rot_array(t1,t2,t3)
    rot(m,rot_arr)

    # Write cif file to screen
    #m.write_cif()
    #m.write_our_xyz()
    #m.write_real_xyz()

if __name__ == "__main__":
    main()
