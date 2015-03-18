import sys
from atom import Atom
from model import Model
import operator

def frac_to_cart(x,y,z, lx,ly,lz):
    """ Assuming x,y,z, are initially between 0 and 1 """
    x = x*lx - lx/2.0
    y = y*ly - ly/2.0
    z = z*lz - lz/2.0
    return x,y,z


def dump_to_modelfiles(dump,masses,base_modelname):
    flag = None
    for line in open(dump):
        if('ITEM: TIMESTEP' in line):
            if(flag):
                # Save the model
                comment = 'timestep {0}'.format(timestep)
                lx = bounds[0]
                ly = bounds[1]
                lz = bounds[2]
                m = Model(comment,lx,ly,lz, atoms)
                m.atoms = sorted(m.atoms, key=operator.itemgetter('id'))
                filename = base_modelname + '_' + str(timestep) + '.xyz'
                m.write_real_xyz(filename)
            flag = 'TIMESTEP'
            continue
        elif('ITEM: NUMBER OF ATOMS' in line):
            flag = 'NUMBER OF ATOMS'
            continue
        elif('ITEM: BOX BOUNDS pp pp pp' in line):
            flag = 'BOX BOUNDS'
            bounds = []
            continue
        elif('ITEM: ATOMS id type xs ys zs' in line):
            flag = 'ATOMS'
            atoms = []
            continue

        if(flag == 'TIMESTEP'):
            timestep = int(line.strip())
        elif(flag == 'NUMBER OF ATOMS'):
            natoms = int(line.strip())
        elif(flag == 'BOX BOUNDS'):
            bounds.append( [float(x) for x in line.strip().split()] )
            bounds[-1] = sum([abs(x) for x in bounds[-1]])
        elif(flag == 'ATOMS'):
            line = line.strip().split()
            id = int(line[0])
            znum = masses[ int(line[1]) ]
            x,y,z = frac_to_cart(float(line[2]), float(line[3]), float(line[4]), bounds[0],bounds[1],bounds[2])
            atoms.append( Atom(id,znum,x,y,z) )
    # Save the last model
    comment = 'timestep {0}'.format(timestep)
    lx = bounds[0]
    ly = bounds[1]
    lz = bounds[2]
    m = Model(comment,lx,ly,lz, atoms)
    m.atoms = sorted(m.atoms, key=operator.itemgetter('id'))
    filename = base_modelname + '_' + str(timestep) + '.xyz'
    m.write_real_xyz(filename)




def main():
    dumpfile = sys.argv[1]
    masses = {
      1:46,
      2:14 }
    dump_to_modelfiles(dumpfile, masses, 'modelfiles/md_model')



if __name__ == '__main__':
    main()
