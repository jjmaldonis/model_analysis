import sys
import copy
import random
import math
import numpy as np
from math import cos,sin,sqrt,pi


def dist2(atom1,atom2):
    return (atom1.coord[0]-atom2.coord[0])**2 + (atom1.coord[1]-atom2.coord[1])**2 + (atom1.coord[2]-atom2.coord[2])**2


class MonteCarlo2(object):
    def __init__(self, m1, m2, temperature, unitvec):
        self.m1 = m1
        self.m2 = m2
        self.theta = 180.
        self.unitvec = unitvec
        self.value = 1000000000
        self.temperature = temperature
        self.numsteps = 0

        self._m1 = copy.deepcopy(m1)
        self._m2 = copy.deepcopy(m2)
        self._theta = 0.

    def rot_point(self, point, linepoint, unitvec, theta):
        # See last result:
        # http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/ArbitraryAxisRotation.pdf
        # theta must be in radians
        # unitvec is the line direction, it will be normalized to a unit vector no matter what
        # linepoint is a point on the line through which to rotate
        # point is the x,y,z point you want to rotate
        x,y,z = point
        if x==y==z==0.: return x,y,z
        a,b,c = linepoint
        u,v,w = unitvec
        L = sqrt( sum(unitvec[i]*unitvec[i] for i in range(len(unitvec))) )
        u /= L
        v /= L
        w /= L

        nx = (a*(v**2 + w**2) - u*(b*v + c*w - u*x - v*y - w*z)) * (1 - cos(theta)) + x*cos(theta) + (-c*v + b*w - w*y + v*z)*sin(theta)

        ny = (b*(u**2 + w**2) - v*(a*u + c*w - u*x - v*y - w*z)) * (1 - cos(theta)) + y*cos(theta) + ( c*u - a*w + w*x - u*z)*sin(theta)

        nz = (c*(u**2 + v**2) - w*(a*u + b*v - u*x - v*y - w*z)) * (1 - cos(theta)) + z*cos(theta) + (-b*u + a*v - v*x + u*y)*sin(theta)
        return nx,ny,nz



    def rot(self, model, theta):
        for atom in model.atoms:
            atom.coord = self.rot_point(atom.coord, (0,0,0), self.unitvec, theta)

    def cost_func(self):
        #print(dist2(self.m1.atoms[0], self.m2.atoms[0]))
        #print(self.m1.atoms[0].coord)
        #print(self.m2.atoms[0].coord)
        #print(self.rot_arr)
        #return max( dist2(self.m1.atoms[a1], self.m2.atoms[a2]) for a1 in range(self.m1.natoms) for a2 in range(self.m2.natoms) )
        return max( dist2(a1,a2) for a1 in self.m1.atoms for a2 in self.m2.atoms )
        #s = 0.
        #for a1,a2 in self.map.items():
        #    s += dist2(self.m1.atoms[a1], self.m2.atoms[a2])
        #return s

    def accept(self):
        self._theta = self.theta
        pass

    def reject(self):
        self.theta = self._theta
        pass

    def step_forward(self):
        s = 1.
        tstep = random.uniform(-s, s)
        #tstep = 1.e-3
        self.theta += tstep
        if(self.theta > 360): self.theta -= 360
        if(self.theta < 0):   self.theta += 360
        self.m2 = copy.deepcopy(self._m2)
        self.rot(self.m2, self.theta)
        self.numsteps += 1

    def run(self):
        self.step_forward()
        val = self.cost_func()
        delta = val - self.value
        rand = random.random()
        if delta < 0 or math.log(1-rand) < -delta/self.temperature:
            #if delta >= 0:
            #    print("{0}\tAccepted! Delta = {1} and {2} < {3}".format(self.numsteps, delta, math.log(1-rand), -delta/self.temperature))
            self.accept()
            self.value = val
        else:
            self.reject()
            #print("{0}\tRejected! Delta = {1} and {2} > {3}".format(self.numsteps, delta, math.log(1-rand), -delta/self.temperature))
            #self.accept()
            #self.value = val






