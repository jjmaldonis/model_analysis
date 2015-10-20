import sys
from model import Model
from tools import printTable
from rot_3d import rotate
from categorize_vor import voronoi_3d,load_param_file,set_atom_vp_types,vor_stats,index_stats
from recenter_model import recenter_model
from math import sqrt


def convert(m,tablename,filename):
    of = open(filename,'w')
    of.write('*set Dim /1*3/;\n')
    of.write('*set Coords /1*{0}/;\n'.format(m.natoms+1))
    of.write('table {0}(Coords,Dim)\n'.format(tablename))
    table = [ [0 for i in range(0,3+1)] for j in range(0,m.natoms+1)]
    #of.write('\t'.join(range(1,m.natoms+1) + '\n')
    for i in range(0,m.natoms+1):
        table[i][0] = i
    for i in range(1,4):
        table[0][i] = i
    for i,atom in enumerate(m.atoms):
        for j in range(0,3):
            table[i+1][j+1] = atom.coord[j]
    table[0][0] = '_'
    ptable = printTable(table,True)
    ptable = ptable.replace('_',' ')
    of.write(ptable)
    of.write(';')
    of.close()


def main():
    modelfile = sys.argv[1]
    paramfile = sys.argv[2]

    m = Model(modelfile)
    m.generate_neighbors(3.45)

    cutoff = {}
    cutoff[(40,40)] = 3.6
    cutoff[(13,29)] = 3.6
    cutoff[(29,13)] = 3.6
    cutoff[(40,13)] = 3.6
    cutoff[(13,40)] = 3.6
    cutoff[(29,40)] = 3.6
    cutoff[(40,29)] = 3.6
    cutoff[(13,13)] = 3.6
    cutoff[(29,29)] = 3.6

    cutoff[(41,41)] = 3.7
    cutoff[(28,28)] = 3.7
    cutoff[(41,28)] = 3.7
    cutoff[(28,41)] = 3.7

    cutoff[(46,46)] = 3.45
    cutoff[(14,14)] = 3.45
    cutoff[(46,14)] = 3.45
    cutoff[(14,46)] = 3.45

    m.generate_neighbors(cutoff)

    voronoi_3d(m,cutoff)

    vp_dict = load_param_file(paramfile)
    set_atom_vp_types(m,vp_dict)

    #vor_stats(m) # Prints what you probably want
    #index_stats(m)
    count = 0
    for atom in m.atoms:
        #print(atom.vp.index[:4])
        if(atom.vp.index[:4] == (0,1,10,2)):
            atoms = [a for a in atom.neighs]+[atom]
            for atom in atoms[1:]:
                #if( abs(atom.coord[0]-atoms[0].coord[0]) > 10 ):
                fix_atom(m, atoms[0], atom)
            vp = Model('0,1,10,2 vp',100,100,100,atoms)
            recenter_model(vp)
            vp.write('vp_data/vp{0}.xyz'.format(count))
            convert(vp,'polyhedron','vp_data/vp{0}.txt'.format(count))
            count += 1
            #for a in vp.atoms:
            #    print(a.coord)
            #break


    #convert(m,'y',modelfile[:-4]+'.txt')
    #angle = 15
    #for atom in m.atoms:
    #    atom.coord = rotate(atom.coord,angle,'x')
    #m.write_real_xyz(modelfile[:-4]+'_rotated_'+str(angle)+'.xyz')
    #convert(m,'x',modelfile[:-4]+'_rotated_{0}.txt'.format(angle))


def dist(a1,a2):
    x1,y1,z1 = a1.coord
    x2,y2,z2 = a2.coord
    return sqrt( (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2 )

def fix_atom(model, ref, tomove):
    eps = 0.0000001
    d2 = model.dist2(ref,tomove)
    if(d2-eps < (ref.coord[0]-tomove.coord[0])**2 + (ref.coord[1]-tomove.coord[1])**2 + (ref.coord[2]-tomove.coord[2])**2 < d2+eps):
        return tomove
    else: # move tomove over the pbcs
        #print("NEED TO FIX SOMETHING!")
        if(ref.coord[0]-tomove.coord[0] > 10):
            #print(dist(ref,tomove))
            tomove.coord = (tomove.coord[0]+model.lx, tomove.coord[1], tomove.coord[2],)
            #print("Fixed x. New dist = {0}".format(dist(ref,tomove)))
        if(ref.coord[0]-tomove.coord[0] < -10):
            #print(dist(ref,tomove))
            tomove.coord = (tomove.coord[0]-model.lx, tomove.coord[1], tomove.coord[2],)
            #print("Fixed x. New dist = {0}".format(dist(ref,tomove)))
        if(ref.coord[1]-tomove.coord[1] > 10):
            #print(dist(ref,tomove))
            tomove.coord = (tomove.coord[0], tomove.coord[1]+model.ly, tomove.coord[2],)
            #print("Fixed y. New dist = {0}".format(dist(ref,tomove)))
        if(ref.coord[1]-tomove.coord[1] < -10):
            #print(dist(ref,tomove))
            tomove.coord = (tomove.coord[0], tomove.coord[1]-model.ly, tomove.coord[2],)
            #print("Fixed y. New dist = {0}".format(dist(ref,tomove)))
        if(ref.coord[2]-tomove.coord[2] > 10):
            #print(dist(ref,tomove))
            tomove.coord = (tomove.coord[0], tomove.coord[1], tomove.coord[2]+model.lz,)
            #print("Fixed z. New dist = {0}".format(dist(ref,tomove)))
        if(ref.coord[2]-tomove.coord[2] < -10):
            #print(dist(ref,tomove))
            tomove.coord = (tomove.coord[0], tomove.coord[1], tomove.coord[2]-model.lz,)
            #print("Fixed z. New dist = {0}".format(dist(ref,tomove)))
        return tomove


if __name__ == '__main__':
    main()
