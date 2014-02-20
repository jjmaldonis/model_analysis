import sys
import os
import random
import math
import numpy as np
from model import Model




# Everything seems to work execpt the volume calculation. nloop is wrong, and I cant figure out how.
# This could also pose a problem for the 'connects' stuff ~ line 90.





def frac(atom,model):
    return (atom.coord[0]/model.lx+0.5, atom.coord[1]/model.ly+0.5, atom.coord[2]/model.lz+0.5)

def voronoi_3d(model,cutoff):
    n0 = 50000
    maxcan = 75
    maxver = 100
    maxepf = 20
    atol = 0.01
    tol = 0.01

    nnab = []
    vvol = []
    nbad = 0
    ibad = []
    nablst = []
    nedges = []

    nedges,nnab,nablst,vvol = vtanal(model,cutoff,nedges,nnab,nablst,vvol,ibad,nbad,tol,atol)
    outvt(model,nedges,nnab,nablst,vvol)


def vtanal(model,cutoff,nedges,nnab,nablst,vvol,ibad,nbad,tol,atol):
    sumvol = 0.0
    volratio = 0.0
    vol = model.lx*model.ly*model.lz
    # model.atomtypes has integer keys, not stings
    radii_sp = model.atomtypes.copy()
    for key in radii_sp:
        radii_sp[key] = 1.0
    for i,atomi in enumerate(model.atoms):
        p = []
        mtag = []
        print("Calculating VP for atom {0}".format(i))
        noi = i
        vvol.append(0.0)
        #print("  Selecting candidates")
        for j,atomj in enumerate(model.atoms):
            if(i != j):
                rxij = frac(atomj,model)[0] - frac(atomi,model)[0]
                ryij = frac(atomj,model)[1] - frac(atomi,model)[1]
                rzij = frac(atomj,model)[2] - frac(atomi,model)[2]
                rxij = rxij - round(rxij)
                ryij = ryij - round(ryij)
                rzij = rzij - round(rzij)
                # This line implements weighted voronoi analysis
                ratio = 2*radii_sp[atomi.z]/(radii_sp[atomi.z]+radii_sp[atomj.z])
                rxij=rxij*ratio*model.lx
                ryij=ryij*ratio*model.ly
                rzij=rzij*ratio*model.lz
                rijsq = rxij**2 + ryij**2 + rzij**2
                # Can set a cutoff for each species, but I dont implement it.
                if(rijsq < cutoff**2):
                    p.append([rxij, ryij, rzij, rijsq]) #???
                    mtag.append(j) #???

        # Candidates have been selected
        nc = len(p)
        #print("Sorting mtag and p")
        # Sort mtag and p
        mtag = [x for (y,x) in sorted(zip(p,mtag),key=lambda x: x[0][3])]
        p.sort(key=lambda x: x[3])
        #print("Calling work")
        nv,nf,ne,nbad,ibad,nepf,nloop,mvijk,v = work(nc,tol,nbad,ibad,p)

        # Sort vertices ring in faces and calculate edge length and face are (???)
        tarea = 0.0
        avglen = 0.0
        avgarea = 0.0
        sleng = []
        tleng = []
        area = []

        #print("  Connects modify nloop")
        for ic in range(0,nc):
            if(nepf[ic] != 0):
                for ie in range(0,nepf[ic]):
                    if(ie == nepf[ic]):
                        iv = nloop[ie][ic]
                        i1 = nloop[1][ic]
                    elif(ie == nepf[ic]-1):
                        iv = nloop[ie][ic]
                        iv1 = nloop[ie+1][ic]
                    else:
                        iv = nloop[ie][ic]
                        iv1 = nloop[ie+1][ic]
                        if( not connect(iv,iv1,mvijk)):
                            for je in range(ie+1,nepf[ic]):
                                jv = nloop[je][ic]
                                if(connect(iv,jv,mvijk)):
                                    nloop[ie+1][ic] = jv
                                    nloop[je][ic] = iv1
                                    break
        
        #print("  Volume and area calculations")
        for ic in range(0,nc):
            sleng.append([])
            tleng.append(0.0)
            area.append(0.0)
            if(nepf[ic] != 0):
                for j in range(0,nepf[ic]+1):
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
                #print('DEBUG1',area[ic])
                tarea += area[ic]
                vvol[i] += area[ic]*math.sqrt(p[ic][3])/6.0
        sumvol += vvol[i]

        # drop small faces / edges, optional
        avgarea = tarea/nf
        for ic in range(0,nc):
            if(nepf[ic] != 0):
                if( (area[ic] != 0) and (area[ic] < atol*tarea) ):
                    break
                avglen = tleng[ic] / float(nepf[ic])
        
        #print("  Generating nablst")
        nedges.append([])
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

    volratio = sumvol/vol
    print("percentages of volume counted: {0}".format(volratio))
        
    return nedges,nnab,nablst,vvol
    #End of vtanal


def work(nc,tol,nbad,ibad,p):
    """ returns (nv,nf,ne,nbad,ibad,nepf,nloop,mvijk,v) """
    if(nc < 4): raise Exception("Less than 4 points given to work: {0}".format(nc))
    mvijk = []
    v = []
    iv = 0
    for i in range(0,nc-2):
        ai = p[i][0]
        bi = p[i][1]
        ci = p[i][2]
        di = -p[i][3]
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
                    while( ok and (l < nc) ):
                        if( (l!=i) and (l!=j) and (l!=k) ):
                            p[l]
                            ok = ( (p[l][0]*vxijk+p[l][1]*vyijk+p[l][2]*vzijk) <= p[l][3] )
                        l += 1
                    if(ok):
                        mvijk.append([i,j,k])
                        v.append([0.5*vxijk, 0.5*vyijk, 0.5*vzijk])
                        #print('DEBUG2',iv,vxijk,vyijk,vzijk)
                        iv += 1
    nv = iv
    if(nv < 4): raise Exception("Less than 4 vertices found in work: {0}".format(nv))

    nepf = [0]*nc
    nloop = np.zeros((75,nv)).tolist() # ??? Change 75 to something more concrete, or use append somehow
    for ii in range(0,len(nloop)):
        for jj in range(0,len(nloop[ii])):
            nloop[ii][jj] = int(nloop[ii][jj])
    for iv in range(0,nv):
        nepf[mvijk[iv][0]] += 1
        nepf[mvijk[iv][1]] += 1
        nepf[mvijk[iv][2]] += 1
        nloop[nepf[mvijk[iv][0]]-1][mvijk[iv][0]] = iv
        nloop[nepf[mvijk[iv][1]]-1][mvijk[iv][1]] = iv
        nloop[nepf[mvijk[iv][2]]-1][mvijk[iv][2]] = iv

    nf = 0
    ne = 0
    for i in range(0,nc):
        if(nepf[i] > 0): nf += 1
        ne = ne + nepf[i]
    if( ne%2 != 0): raise Exception("Something got screwed up in work!")
    ne = ne/2

    nepf = [x-1 for x in nepf]
    
    if(nv-ne+nf != 2):
        nbad += 1
        ibad.append[noi]
    return (nv,nf,ne,nbad,ibad,nepf,nloop,mvijk,v)
    #End of work()


def outvt(model,nedges,nnab,nablst,vvol):
    nnabsp = []
    for i,atomi in enumerate(model.atoms):
        nnabsp.append({})
        for key in model.atomtypes:
            nnabsp[i][key] = 0
        atomi.vp.index = [0,0,0,0,0,0,0,0]

        for j in range(0,nnab[i]):
            nnabsp[i][model.atoms[nablst[i][j]].z] += 1
            # nedges is one off, thats why we start at 3 and not 3
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

        atomi.vp.nnabsp = nnabsp[i]
        atomi.vp.neighs = nablst[i]

        keys = nnabsp[i].keys()
        #print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}".format(i,atomi.z,nnab[i],nnabsp[i][keys[0]],nnabsp[i][keys[1]],nnabsp[i][keys[2]],atomi.vp.index,vvol[i]))


def connect(i,j,mvijk):
    np = 0
    for ii in range(0,3):
        for jj in range(0,3):
            if(mvijk[i][ii] == mvijk[j][jj]): np += 1

    if(np == 0): raise Exception("Not connected, same ic?")
    if(np == 1): return True
    if(np == 2): return False
    if(np == 3): raise Exception("Wrong connection")


def main():
    modelfile = sys.argv[1]
    cutoff = 3.5
    m = Model(modelfile)
    voronoi_3d(m,cutoff)

if __name__ == '__main__':
    main()
