import random
import math

def main():
    filename = 'Zr50Cu45Al5_156atoms.xyz'
    natoms = 156
    composition = {'Zr': 0.50,  'Cu': 0.45, 'Al':0.05}
    cutoff = 2.3
    xsize = 10*math.sqrt(2)
    ysize = 10*math.sqrt(2)
    zsize = 10*math.sqrt(2)
    m = generate_model(filename, natoms, composition, xsize, ysize, zsize, cutoff)


def generate_model(filename, natoms, composition, xsize, ysize, zsize, cutoff):
    cutoff2 = cutoff * cutoff

    natoms_per_elem = {k: int(round(v*natoms)) for k,v in composition.items()}
    natoms_per_elem = sorted([(v,k) for k,v in natoms_per_elem.items()])  # Sort composition from lowest percentage to highest
    for v,k in natoms_per_elem:
        print('{0}: {1} atoms'.format(k,v))
    assert sum(v for v,_ in natoms_per_elem) == natoms

    atoms = []
    content = []
    content.append('{0}\n'.format(natoms))
    content.append( str(xsize) + ' ' + str(ysize) + ' ' + str(zsize) + "\n")

    for v,k in natoms_per_elem:
        i = 0
        while i < v:
            good = True
            x = random.uniform(-xsize/2,xsize/2)
            y = random.uniform(-ysize/2,ysize/2)
            z = random.uniform(-zsize/2,zsize/2)
            for atom in atoms:
                xx = atom[1]-x
                yy = atom[2]-y
                zz = atom[3]-z
                xx = xx-xsize*int(round(xx/xsize))
                yy = yy-ysize*int(round(yy/ysize))
                zz = zz-zsize*int(round(zz/zsize))
                r2 = xx**2 + yy**2 + zz**2
                if( r2 < cutoff2 ):
                    good = False
                    break
            if(good):
                atoms.append( (k, x,y,z) )
                i += 1
            else:
                print("Stuck at "+ str(len(atoms)))

    for atom in atoms:
        content.append( str(atom[0]) + " " + str(atom[1]) + " " + str(atom[2]) + " " + str(atom[3]) + "\n" )
    
    with open(filename, 'w') as f:
        for line in content:
            f.write(line)

    

if __name__ == "__main__":
    main();
