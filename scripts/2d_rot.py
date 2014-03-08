import numpy as np
import sys

from model import Model


def rot(model, arr):
    """ arr should be a 9 element rotation numpy array, which we will reshape here """
    arr = arr.reshape((3,3))
    arr = np.linalg.inv(arr)
    
    for i,atom in enumerate(model.atoms):
        old_coord = [atom.coord[0], atom.coord[1], atom.coord[2]]
        new_coord = np.dot(np.asarray(old_coord), arr)
        atom.set_coord(new_coord[0],new_coord[1],0)#new_coord[2])

def main():
    modelfile = sys.argv[1]
    m = Model(modelfile)

    # Below is the identity rotation array.
    #rot_arr = [1,0,0,0,1,0,0,0,1]

    # Below is the correct rotation matrix for JWH t3 icofrac to get it into the orientation in his PRL paper.
    #rot_arr = [-0.031777, 0.998843, 0.036102, 0.986602, 0.025563, 0.161133, 0.160023, 0.040739, -0.986272]

    # Below is a (the?) rotation matrix for JWH t1 icofrac.
    #rot_arr = [ 0.954646, -0.233932, 0.184194, 0.280650, 0.913581, -0.294287, -0.099433, 0.332633, 0.937800 ]

    # Below is a (the?) rotation matrix of Pei's t1 that gives some planes. Oriented for a specific plane ~.
    #rot_arr = [ -0.977103, -0.123352, -0.173361, -0.130450, 0.990997, 0.030118, 0.168085, 0.052043, -0.984398 ]

    # Below is a (the?) rotation matrix of Pei's t2 that gives some planes. Oriented for a specific plane ~.
    rot_arr = [ 0.985478, -0.010230, -0.169493, 0.009247, 0.999936, -0.006586, 0.169549, 0.004923, 0.985509]

    # Below is a (the?) rotation matrix of Pei's t3 that gives some planes. Oriented for a specific plane ~.
    #rot_arr = [0.981624,-0.002765,-0.190808, -0.003436,0.999477,-0.032163, 0.190797,0.032228,0.981100]

    npra = np.asarray(rot_arr)
    rot(m,npra)

    # Write cif file to screen
    #m.write_cif()
    m.write_our_xyz()

if __name__ == "__main__":
    main()
