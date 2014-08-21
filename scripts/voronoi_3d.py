import sys
import os
import random
import math,time
import numpy as np
from model import Model


class MyError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


def frac(atom,model):
    # Returns the coordinate of atom in fractional coordinates between 0 and 1
    return (atom.coord[0]/model.lx+0.5, atom.coord[1]/model.ly+0.5, atom.coord[2]/model.lz+0.5)

def voronoi_3d(model,cutoff):
    n0 = 50000
    maxcan = 75
    maxver = 100
    maxepf = 20
    atol = 0.03
    tol = 0.03
    tltol = 0.03

    nnab = []
    vvol = np.zeros(model.natoms)
    nbad = 0
    ibad = []
    nablst = []
    nedges = []

    nedges,nnab,nablst,vvol = vp_analysis(model,cutoff,nedges,nnab,nablst,vvol,ibad,nbad,tol,atol,tltol)
    save_vp_atom_data(model,nedges,nnab,nablst,vvol)


def vp_analysis(model,cutoff,nedges,nnab,nablst,vvol,ibad,nbad,tol,atol,tltol):
    sumvol = 0.0
    volratio = 0.0
    vol = model.lx*model.ly*model.lz
    # model.atomtypes has integer keys, not stings
    weight_sp = model.atomtypes.copy()
    for key in weight_sp:
        weight_sp[key] = 1.0
    print_percent = 0.0
    for i,atomi in enumerate(model.atoms):
        p = []
        mtag = []
        #print("Calculating VP for atom {0}".format(i))
        if(100.0*i/model.natoms > print_percent):
            print("{0}% done...".format(print_percent))
            print_percent += 5.0
        noi = i # Number of i, just make a copy of it
        #print("  Selecting candidates")
        for j,atomj in enumerate(model.atoms):
            if(i != j):
                rxij = atomj.coord[0]/model.lx - atomi.coord[0]/model.lx
                ryij = atomj.coord[1]/model.ly - atomi.coord[1]/model.ly
                rzij = atomj.coord[2]/model.lz - atomi.coord[2]/model.lz
                rxij = rxij - round(rxij) #PBCs
                ryij = ryij - round(ryij)
                rzij = rzij - round(rzij)
                # These four lines implement weighted voronoi analysis
                ratio = 2*weight_sp[atomi.z]/(weight_sp[atomi.z]+weight_sp[atomj.z])
                rxij=rxij*ratio*model.lx
                ryij=ryij*ratio*model.ly
                rzij=rzij*ratio*model.lz
                rijsq = rxij**2 + ryij**2 + rzij**2
                # Select all atoms within cutoff of atomi
                # Can set a cutoff for each species, but I dont implement it.
                try:
                    thiscut = cutoff[(atomi.z,atomj.z)]
                    if(rijsq < thiscut**2):
                        p.append([rxij, ryij, rzij, rijsq]) # I dont know why we save this information, but it is the distance information
                        mtag.append(j) #???
                except TypeError:
                    if(rijsq < cutoff**2):
                        p.append([rxij, ryij, rzij, rijsq]) # I dont know why we save this information, but it is the distance information
                        mtag.append(j) #???

        # Candidates have been selected
        nc = len(p)
        #print("Sorting mtag and p")
        # Sort mtag and p
        mtag = [x for (y,x) in sorted(zip(p,mtag),key=lambda x: x[0][3])]
        p.sort(key=lambda x: x[3])
        #print("Calling work")
        #nv,nf,ne,nbad,ibad,nepf,nloop,mvijk,v = work(nc,tol,nbad,ibad,p)
        try:
            nv,nf,ne,nepf,nloop,mvijk,v = work(nc,tol,nbad,ibad,p)
            good = True
        except MyError as inst:
            good = False
            print(inst)
            nv = 0
            nf = 0
            ne = 0
            nepf = []
            nloop = []
            mvijk = []
            nc = 0
            p = []
            mtag = []

        if(good):
            # Do the area calculation, which currently doesnt work.
            # This could be a problem TODO
            # Sort vertices ring in faces and calculate edge length and face are (???)
            tarea = 0.0
            avglen = 0.0
            avgarea = 0.0
            sleng = []
            tleng = []
            area = []

            #for ic in range(0,nc):
                #print(nepf[ic])
                #for ie in range(0,nepf[ic]+1):
                #    print("{0} {1} {2}".format(ie,ic,nloop[ie][ic]))

            #print(nloop[4][2])
            #print("  Connects modify nloop")
            # For each connection
            for ic in range(0,nc):
                if(nepf[ic] >= 0): #This is unnecessary
                    #print(ic,nepf[ic])
                    for ie in range(0,nepf[ic]+1):
                        #print(ic,ie)
                        #if(ic == 4 and ie == 2):
                        #    iv = nloop[ie][ic]
                        #    iv1 = nloop[ie+1][ic]
                        #    print("{0} {1} {2} {3}".format(iv,iv1,mvijk[iv],mvijk[iv1]))
                        if(ie == nepf[ic]):
                            #print(ic,ie)
                            #print("{0} {1} DEBUG 1".format(ic,ie))
                            iv = nloop[ie][ic]
                            i1 = nloop[1][ic]
                            #print(ie+1,ic)
                            #print(nloop[ie][ic]) # TODO
                            #print("{0} {1} {2}".format(ie,ic,nloop[ie][ic]))
                            #print(nloop[ie])
                        elif(ie == nepf[ic]-1):
                            #print("{0} {1} DEBUG 2".format(ic,ie))
                            iv = nloop[ie][ic]
                            iv1 = nloop[ie+1][ic]
                        else:
                            #print("{0} {1} DEBUG 3".format(ic,ie))
                            iv = nloop[ie][ic]
                            iv1 = nloop[ie+1][ic]
                            #print("{0} {1} {2} {3}".format(iv,iv1,mvijk[iv],mvijk[iv1]))
                            if( not connect(iv,iv1,mvijk)):
                                #print("{2} {3} Not connect\t{0}\t{1}".format(iv,iv1,ic,ie))
                                for je in range(ie+2,nepf[ic]+1):
                                    #print("je=",je)
                                    jv = nloop[je][ic]
                                    if(connect(iv,jv,mvijk)):
                                        nloop[ie+1][ic] = jv
                                        nloop[je][ic] = iv1
                                        #print("Connected",iv,jv)
                                        #print("Set nloop[{0}][{1}]={2}".format(ie+1,ic,nloop[ie+1][ic]))
                                        #print("Set nloop[{0}][{1}]={2}".format(je,ic,nloop[je][ic]))
                                        break
                                    #else:
                                        #print("iv,jv not connected",iv,jv)
            
            #for ic in range(0,nc):
            #    #print(nepf[ic])
            #    for ie in range(0,nepf[ic]+1):
            #        print("{0} {1} {2}".format(ie,ic,nloop[ie][ic]))
                
            #print("  Volume and area calculations")
            for ic in range(0,nc):
                sleng.append([])
                tleng.append(0.0)
                area.append(0.0)
                if(nepf[ic] != 0):
                    for j in range(0,nepf[ic]):
                        #print('DEBUG-3',ic,j)
                        ivs = nloop[j][ic]
                        if(j == nepf[ic]):
                            ive = nloop[1][ic]
                        else:
                            ive = nloop[j+1][ic]
                        sleng[ic].append( math.sqrt((v[ivs][0]-v[ive][0])**2+(v[ivs][1]-v[ive][1])**2+(v[ivs][2]-v[ive][2])**2) )
                        tleng[ic] += sleng[ic][j]
                        x1 = v[ivs][0] - v[nloop[0][ic]][0]
                        y1 = v[ivs][1] - v[nloop[0][ic]][1]
                        z1 = v[ivs][2] - v[nloop[0][ic]][2]
                        x2 = v[ive][0] - v[nloop[0][ic]][0]
                        y2 = v[ive][1] - v[nloop[0][ic]][1]
                        z2 = v[ive][2] - v[nloop[0][ic]][2]
                        #print('DEBUG0',ic,nloop[0][ic],v[ivs][0],v[nloop[0][ic]][0])
                        #print('DEBUG1',ic,nloop[0][ic],v[ivs][1],v[nloop[0][ic]][1])
                        #print('DEBUG2',ic,nloop[0][ic],v[ivs][2],v[nloop[0][ic]][2])
                        #print('DEBUG-2',x1,y1,z1,x2,y2,z2)
                        #print('DEBUG-1',(y1*z2-z1*y2)**2,(z1*x2-z2*x1)**2,(x1*y2-x2*y1)**2)
                        #print(ic,j,ivs,ive)
                        #####print(ic,j,nloop[1][ic]) # THIS LINE IS THE BEST CURRENTLY FOR COMPARING TO vor_v4.f90
                        #print('DEBG0',ic,j,0.5*math.sqrt( (y1*z2-z1*y2)**2 + (z1*x2-z2*x1)**2 + (x1*y2-x2*y1)**2 ))
                        area[ic] += 0.5*math.sqrt( (y1*z2-z1*y2)**2 + (z1*x2-z2*x1)**2 + (x1*y2-x2*y1)**2 )
                        #print('DEBUG 0',area[ic])
                    #print(round(area[ic],6))
                    tarea += area[ic]
                    vvol[i] += area[ic]*math.sqrt(p[ic][3])/6.0
            sumvol += vvol[i]
            #print(vvol[i])

            # drop small faces / edges, optional
            avgarea = tarea/nf
            for ic in range(0,nc):
                if(nepf[ic] != 0):
                    if( (area[ic] != 0) and (area[ic] < atol*tarea) ):
                    #    for j in range(0,len(sleng)):
                    #        try:
                    #            sleng[j][ic] = 0
                    #        except:
                    #            pass
                    #    print("Dropped a face!")
                        break
                    avglen = tleng[ic] / float(nepf[ic])
                    #print(nepf[ic], sleng[j])
                    #for j in range(0,nepf[ic]):
                    #    try:
                    #        if((sleng[j][ic] != 0.0) and (sleng[j][ic] < tltol*avglen)):
                    #            sleng[j][ic] = 0
                    #            print("Dropped an edge!")
                    #    except:
                    #        pass
            
        #print("  Generating nablst")
        # nedges will create the vp indexes
        nedges.append([]) # Create the ith position in nedges
        nnab.append(0)
        nablst.append([])
        for ic in range(0,nc):
            nedges[i].append(0)
            if(nepf[ic] != 0):
                for j in range(0,nepf[ic]):
                    if(sleng[ic][j] != 0):
                        nedges[i][ic] += 1 #???
                if(nedges[i][ic] != 0):
                    nnab[i] += 1 #???
                    nablst[i].append(0)
                    nablst[i][nnab[i]-1] = mtag[ic]
        nedges[i] = [x for x in nedges[i] if x != 0]

    volratio = sumvol/vol
    print("percentages of volume counted: {0}".format(volratio))
        
    return nedges,nnab,nablst,vvol
    #End of vp_analysis


def work(nc,tol,nbad,ibad,p):
    """ returns (nv,nf,ne,nbad,ibad,nepf,nloop,mvijk,v) """
    if(nc < 4): raise MyError("Less than 4 points given to work: {0}".format(nc))
    mvijk = [] # This will be a 2D matrix of size nv x 3
    v = []
    iv = 0
    for i in range(0,nc-2):
        ai = p[i][0] #rx for connection i
        bi = p[i][1] #ry for connection i
        ci = p[i][2] #rz for connection i
        di = -p[i][3] #dist for connection i
        for j in range(i+1,nc-1):
            aj = p[j][0]
            bj = p[j][1]
            cj = p[j][2]
            dj = -p[j][3]
            ab = ai * bj - aj * bi
            bc = bi * cj - bj * ci
            ca = ci * aj - cj * ai
            da = di * aj - dj * ai
            db = di * bj - dj * bi
            dc = di * cj - dj * ci
            for k in range(j+1,nc):
                ak = p[k][0]
                bk = p[k][1]
                ck = p[k][2]
                dk = -p[k][3]
                det = ak * bc + bk * ca + ck * ab
                if ( abs ( det ) > tol ):
                    detinv = 1.0 / det
                    vxijk = ( - dk * bc + bk * dc - ck * db ) * detinv
                    vyijk = ( - ak * dc - dk * ca + ck * da ) * detinv
                    vzijk = (   ak * db - bk * da - dk * ab ) * detinv
                    ok = True
                    l = 0
                    # Go thru all the other connections besides i,j,k
                    while( ok and (l < nc) ):
                        if( (l!=i) and (l!=j) and (l!=k) ):
                            ok = ( (p[l][0]*vxijk+p[l][1]*vyijk+p[l][2]*vzijk) <= p[l][3] )
                        l += 1
                    # If all of those connetions are within dist (p[l][3]) then append a vertex
                    if(ok):
                        mvijk.append([i,j,k])
                        v.append([0.5*vxijk, 0.5*vyijk, 0.5*vzijk])
                        #print('DEBUG2',iv,vxijk,vyijk,vzijk)
                        iv += 1
    # Set the number of vertices found
    nv = iv
    if(nv < 4): raise MyError("Less than 4 vertices found in work: {0}".format(nv))

    nepf = [0]*nc # Number of edges per face
    nloop = np.zeros((75,nv),dtype=int).tolist() # ??? Change 75 to something more concrete, or use append somehow
    # This seems strange but I'm okay with it I think.
    for iv in range(0,nv):
        # This is done in fortran
        nepf[mvijk[iv][0]] += 1
        nepf[mvijk[iv][1]] += 1
        nepf[mvijk[iv][2]] += 1
        nloop[nepf[mvijk[iv][0]]-1][mvijk[iv][0]] = iv
        nloop[nepf[mvijk[iv][1]]-1][mvijk[iv][1]] = iv
        nloop[nepf[mvijk[iv][2]]-1][mvijk[iv][2]] = iv

    nf = sum([1 for x in nepf if x > 0])
    ne = sum(nepf)
    if( ne%2 != 0): raise MyError("Something got screwed up in work! {0}".format(ne))
    ne = ne/2

    # Need to do this.
    nepf = [x-1 for x in nepf]
    
    if(nv-ne+nf != 2):
        raise MyError("Bad atom!")
    #return (nv,nf,ne,nbad,ibad,nepf,nloop,mvijk,v)
    return (nv,nf,ne,nepf,nloop,mvijk,v)
    #End of work()


def save_vp_atom_data(model,nedges,nnab,nablst,vvol):
    nnabsp = []
    for i,atomi in enumerate(model.atoms):
        nnabsp.append({})
        for key in model.atomtypes:
            nnabsp[i][key] = 0
        atomi.vp.index = [0,0,0,0,0,0,0,0]

        for j in range(0,len(nedges[i])):
            if(nedges[i][j] > 0):
                nnabsp[i][model.atoms[nablst[i][j]].z] += 1
                # nedges is one off, thats why we start at 2 and not 3
                if(nedges[i][j] == 2):
                    atomi.vp.index[0] += 1
                elif(nedges[i][j] == 3):
                    atomi.vp.index[1] += 1
                elif(nedges[i][j] == 4):
                    atomi.vp.index[2] += 1
                elif(nedges[i][j] == 5):
                    atomi.vp.index[3] += 1
                elif(nedges[i][j] == 6):
                    atomi.vp.index[4] += 1
                elif(nedges[i][j] == 7):
                    atomi.vp.index[5] += 1
                elif(nedges[i][j] == 8):
                    atomi.vp.index[6] += 1
                elif(nedges[i][j] == 9):
                    atomi.vp.index[7] += 1
                elif(nedges[i][j] == 10):
                    atomi.vp.index[8] += 1

        nablst[i] = [model.atoms[x] for x in nablst[i]]
        atomi.vp.nnabsp = nnabsp[i]
        atomi.vp.neighs = nablst[i]
        atomi.vp.vol = vvol[i]
        if(nnab[i] != sum(atomi.vp.index)):
            print("ERROR!!! nnab is not right for atom {2}! {0} {1}".format(nnab[i],atomi.vp.index,i))


def print_data(model):
    for i,atomi in enumerate(model.atoms):
        keys = atomi.vp.nnabsp.keys()
        s = str( atomi.vp.nnabsp[keys[0]] )
        for j in range(1,len(keys)):
            s += '\t' + str( atomi.vp.nnabsp[keys[j]] )
        #print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(i,atomi.z,nnab[i],s,atomi.vp.index,atomi.vp.vol))
        print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(i,atomi.z,sum(atomi.vp.index),s,atomi.vp.index,atomi.vp.vol))



def connect(i,j,mvijk):
    np = 0
    for ii in range(0,3):
        for jj in range(0,3):
            if(mvijk[i][ii] == mvijk[j][jj]): np += 1

    if(np == 0): raise Exception("Not connected, same ic?")
    if(np == 1): return False
    if(np == 2): return True
    if(np == 3): raise Exception("Wrong connection")


def main():
    # NOTE: Cutoff can either be a single integer or it
    # can be a dictionary where the keys are two-tuples
    # of atomic numbers (e.g. (40,13)=3.5 for Zr,Al).
    modelfile = sys.argv[1]
    m = Model(modelfile)
    try:
        cut = float(sys.argv[2])
        cutoff = {}
        for z1 in m.atomtypes:
            for z2 in m.atomtypes:
                cutoff[(z1,z2)] = cut
                cutoff[(z2,z1)] = cut
    except:
        print("You didn't input a cutoff so you much define it in the code.")
    voronoi_3d(m,cutoff)
    print_data(m)

if __name__ == '__main__':
    main()
