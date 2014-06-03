import numpy as np
import sys

from model import Model


def rot(model, arr):
    """ arr should be a 9 element rotation numpy array, which we will reshape here """
    if( all( i == 0 for i in arr )):
        for i,atom in enumerate(model.atoms):
            atom.set_coord(atom.coord[0],atom.coord[1],0)
    else:
        arr = arr.reshape((3,3))
        arr = np.linalg.inv(arr)
        
        for i,atom in enumerate(model.atoms):
            old_coord = [atom.coord[0], atom.coord[1], atom.coord[2]]
            new_coord = np.dot(np.asarray(old_coord), arr)
            atom.set_coord(new_coord[0],new_coord[1],0)#new_coord[2])

def main():
    modelfile = sys.argv[1]
    outfilebase = sys.argv[2]
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
    #rot_arr = [ 0.985478, -0.010230, -0.169493, 0.009247, 0.999936, -0.006586, 0.169549, 0.004923, 0.985509]

    # Below is a (the?) rotation matrix of Pei's t3 that gives some planes. Oriented for a specific plane ~.
    #rot_arr = [0.981624,-0.002765,-0.190808, -0.003436,0.999477,-0.032163, 0.190797,0.032228,0.981100]

    # Below is a rotation matrix of Pei's Zr50 t3 that gives some planes. Oriented for a specific plane ~.
    #rot_arr = [0.983910,-0.035278,-0.175147,0.048246,0.996356,0.070340,0.172028,-0.077659,0.982026]

    # Below is a rotation matrix of Pei's Zr50 t3, cluster 0 from 'Crystal-like' with an alternate plane pattern.
    #rot_arr = [-0.078966,-0.539181,0.838479,-0.486443,0.755011,0.439695,-0.870137,-0.373151,-0.321901]
    #rot_arr = [-0.550168, -0.592703, -0.588233, -0.675963, 0.729699, -0.103023, 0.490295, 0.340944, -0.802102] # a second attempt
    #rot_arr = [ -0.102056, 0.795154, -0.597758, 0.867018, 0.365702, 0.338439, 0.487712, -0.483728, -0.726735] # eh
    #rot_arr = [ -0.485998, 0.815787, 0.313525, -0.151211, -0.431820, 0.889195, 0.860779, 0.384739, 0.333220]
    #rot_arr = [ -0.508856, -0.631880, -0.584631, -0.710228, 0.691926, -0.129672, 0.486458, 0.349237, -0.800870]

    # Below is a rotation matrix of Pei's Zr50 t3, all clusters, with one promising spot. Mostly noise probably.
    #rot_arr = [0.563829, 0.290821, -0.772995, -0.641563, -0.435169, -0.631684, -0.520090, 0.852086, -0.058781]
    #rot_arr = [-0.832982, -0.489634, -0.257683, -0.212007, -0.147731, 0.966038, -0.511072, 0.859322, 0.019252]
    #rot_arr = [ 0.494247,  0.236924,  -0.836413, 0.704623,  0.454325,  0.545064, 0.509142,  -0.858752,  0.057607]

    #rot_arr = [0.235445, 0.089792, 0.967731, -0.967993, 0.110713, 0.225236, -0.086917, -0.989788, 0.112985]
    #rot_arr = [0.216234, -0.047608, 0.975180,-0.955161, 0.196601, 0.221394,-0.202262, -0.979327, -0.002961]
    # BELOW: t3 - 0 rot
    #rot_arr = [ 0.984104, 0.007169, -0.177450, 0.008311, 0.996231, 0.086340, 0.177400, -0.086442, 0.980335]
    #rot_arr = [ 0.983102, -0.003388, -0.183028, 0.003365, 0.999994, -0.000437, 0.183028, -0.000187, 0.983108]
    
    # Pei's t2
    #rot_arr = [ 0.985568, 0.001345, -0.169276, 0.000677, 0.999929, 0.011883, 0.169280, -0.011826, 0.985497] # 0rot
    #rot_arr = [ -0.406881, -0.461953, -0.788065, -0.548766, 0.813295, -0.193412, 0.730277, 0.353768, -0.584418] #1rot
    #rot_arr =[ -0.466130, 0.820964, 0.329760, -0.171769, -0.449615, 0.876551, 0.867881, 0.351944, 0.350596] # 2rot
    #rot_arr = [ -0.303913, -0.572611, -0.761415, -0.608626, 0.731562, -0.307232, 0.732947, 0.370045, -0.570838] # 3rot
    #rot_arr = [ -0.653632, -0.441386, -0.614771, -0.749000, 0.260854, 0.609061, -0.108465,0.858565, -0.501100] # 4rot

    # Pei's t1
    #rot_arr =  [0.985435, 0.000143, -0.170053, 0.003167, 0.999811, 0.019194, 0.170023, -0.019453, 0.985248] #0rot
    #rot_arr = [ -0.181186, -0.653769, -0.734682, -0.632103, 0.649717, -0.422273, 0.753404, 0.387885, -0.530968] #1rot
    #rot_arr = [ -0.408942, 0.873456, 0.264275, -0.252622, -0.386635, 0.886959, 0.876897, 0.295953, 0.378765] #2rot
    #rot_arr = [ -0.243749, -0.741828, -0.624722, -0.637231, 0.608094, -0.473453, 0.731111, 0.282689, -0.620938] #3rot

    # JWH's t1
    #rot_arr = [ 0.915013, -0.004459, -0.403401, -0.010527, 0.999335, -0.034924, 0.403288, 0.036202, 0.914357]
    # JWH's t2
    #rot_arr = [ 0.983576, 0.009003, -0.180268, 0.010862, 0.993992, 0.108912, 0.180165, -0.109081, 0.977569]
    #JWH's t3
    #rot_arr = [ 0.985667, 0.014122, 0.168112, -0.020755, 0.999071, 0.037763, -0.167422, -0.040711, 0.985044]
    #JWH's t3 all.xtal cluster
    #rot_arr = [ 0.985237, 0.004144, 0.171146, -0.008192, 0.999703, 0.022952, -0.171000, -0.024015, 0.984978] 

    rot_arr = [0,0,0,0,0,0,0,0,0]

    npra = np.asarray(rot_arr)
    rot(m,npra)

    # Write cif file to screen
    m.write_cif(outfilebase+'.cif')
    m.write_our_xyz(outfilebase+'.xyz')

if __name__ == "__main__":
    main()
