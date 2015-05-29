import sys
from model import Model
from tools import printTable
from rot_3d import rotate


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
    m = Model(modelfile)
    convert(m,'y',modelfile[:-4]+'.txt')
    angle = 15
    for atom in m.atoms:
        atom.coord = rotate(atom.coord,angle,'x')
    m.write_real_xyz(modelfile[:-4]+'_rotated_'+str(angle)+'.xyz')
    convert(m,'x',modelfile[:-4]+'_rotated_{0}.txt'.format(angle))



if __name__ == '__main__':
    main()
