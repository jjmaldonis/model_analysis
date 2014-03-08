import sys
import numpy as np
from model import Model
from math import cos, sin

def rot(model, arr):
    """ arr should be a 9 element rotation numpy array, which we will reshape here
        NOTE! This changes the model atom positions in the model """
    if(type(arr) == list): arr = np.array(arr)
    if(arr.shape == (9,)): arr = arr.reshape((3,3))
    arr = np.linalg.inv(arr)
    
    for i,atom in enumerate(model.atoms):
        old_coord = [atom.coord[0], atom.coord[1], atom.coord[2]]
        new_coord = np.dot(np.asarray(old_coord), arr)
        atom.set_coord(new_coord[0],new_coord[1],0)#new_coord[2])


def calc_rot_array(t1,t2,t3):
    """ We construct the rotation matrix based on t1,t2,t3
        NOTE! Order matters! """
    rx = np.array( [ [1,0,0], [0,cos(t1),-sin(t1)], [0,sin(t1),cos(t1)] ] )
    ry = np.array( [ [cos(t2),0,sin(t2)], [0,1,0], [-sin(t2),0,cos(t2)] ] )
    rz = np.array( [ [cos(t3),-sin(t3),0], [sin(t3),cos(t3),0], [0,0,1] ] )
    arr = np.dot( np.dot(ry,rx) ,rz)
    return np.linalg.inv(arr)


def main():
    modelfile = sys.argv[1]
    m = Model(modelfile)

    # Below is a (the?) rotation matrix of Pei's t1 that gives some planes. Oriented for a specific plane ~.
    rot_arr = [ -0.977103, -0.123352, -0.173361, -0.130450, 0.990997, 0.030118, 0.168085, 0.052043, -0.984398 ]
    rot(m,rot_arr)

    # Angles in radians
    # Note that these are semi difficult to figure out from the vesta rotation matrix,
    # partly because there are negative angles, so you may need to do 2pi - angle you found.
    #t1 = np.pi*2 - 0.0371505
    #t2 = 0.162790
    #t3 = 0
    #rot_arr = calc_rot_array(m,t1,t2,t3)
    #rot(m,rot_arr)

    # Write cif file to screen
    #m.write_cif()
    m.write_our_xyz()

if __name__ == "__main__":
    main()
