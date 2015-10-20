import sys
import copy
import random
import math
import numpy as np
from math import cos,sin


def dist2(atom1,atom2):
    return (atom1.coord[0]-atom2.coord[0])**2 + (atom1.coord[1]-atom2.coord[1])**2 + (atom1.coord[2]-atom2.coord[2])**2


class MonteCarlo(object):
    def __init__(self, m1, m2, map, temperature):
        self.m1 = m1
        self.m2 = m2
        self.alpha = 0.
        self.beta = 0.
        self.gamma = 0.
        self.value = 1000000000
        self.temperature = temperature
        self.numsteps = 0

        self.map = map

        self._m1 = copy.deepcopy(m1)
        self._m2 = copy.deepcopy(m2)
        self._alpha = self.alpha
        self._beta = self.beta
        self._gamma = self.gamma

    def rot(self, model, arr):
        """ arr should be a 9 element rotation numpy array, which we will reshape here
            NOTE! This changes the model atom positions in the model """
        if(type(arr) == list): arr = np.array(arr)
        if(arr.shape == (9,)): arr = arr.reshape((3,3))
        #print(arr)
        arr = np.linalg.inv(arr)
        
        for i,atom in enumerate(model.atoms):
            old_coord = [atom.coord[0], atom.coord[1], atom.coord[2]]
            new_coord = np.dot(np.asarray(old_coord), arr)
            atom.coord = (new_coord[0],new_coord[1],new_coord[2])


    def calc_rot_array(self, t1, t2, t3, deg=True):
        """ We construct the rotation matrix based on t1,t2,t3
            NOTE! Order matters! """
        t1 = t1*np.pi/180.0 # in radians
        t2 = t2*np.pi/180.0 # in radians
        t3 = t3*np.pi/180.0 # in radians
        rx = np.array( [ [1,0,0], [0,cos(t1),-sin(t1)], [0,sin(t1),cos(t1)] ] )
        ry = np.linalg.inv(np.array( [ [cos(t2),0,sin(t2)], [0,1,0], [-sin(t2),0,cos(t2)] ] ))
        rz = np.linalg.inv(np.array( [ [cos(t3),-sin(t3),0], [sin(t3),cos(t3),0], [0,0,1] ] ))
        arr = np.dot( np.dot(ry,rx) ,rz)
        return np.linalg.inv(arr)


    @property
    def rot_arr(self):
        return self.calc_rot_array(self.alpha, self.beta, self.gamma)

    def cost_func(self):
        #print(dist2(self.m1.atoms[0], self.m2.atoms[0]))
        #print(self.m1.atoms[0].coord)
        #print(self.m2.atoms[0].coord)
        #print(self.rot_arr)
        return max( dist2(self.m1.atoms[a1], self.m2.atoms[a2]) for a1,a2 in self.map.items() )
        #s = 0.
        #for a1,a2 in self.map.items():
        #    s += dist2(self.m1.atoms[a1], self.m2.atoms[a2])
        #return s

    def accept(self):
        self._alpha = self.alpha
        self._beta = self.beta
        self._gamma = self.gamma
        pass

    def reject(self):
        self.alpha = self._alpha
        self.beta = self._beta
        self.gamma = self._gamma 
        pass

    def step_forward(self):
        s = 1.
        astep = random.uniform(-s, s)
        bstep = random.uniform(-s, s)
        gstep = random.uniform(-s, s)
        self.alpha += astep
        self.beta  += bstep
        self.gamma += gstep
        if(self.alpha > 360): self.alpha -= 360
        if(self.alpha < 0):   self.alpha += 360
        if(self.beta  > 360): self.beta -= 360
        if(self.beta  < 0):   self.beta += 360
        if(self.gamma > 360): self.gamma -= 360
        if(self.gamma < 0):   self.gamma += 360
        self.m2 = copy.deepcopy(self._m2)
        self.rot(self.m2, self.rot_arr)
        self.numsteps += 1

    def run(self):
        self.step_forward()
        val = self.cost_func()
        delta = val - self.value
        #if delta < 0 or math.log(1-random.random()) < -delta/self.temperature:
        if delta < 0:
            self.accept()
            self.value = val
        else:
            self.reject()
            #self.accept()
            #self.value = val






