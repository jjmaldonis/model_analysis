import sys
from model import Model
from recenter_model import recenter_model
from math import sqrt
from random import shuffle

def dist2(atom1,atom2):
    return (atom1.coord[0]-atom2.coord[0])**2 + (atom1.coord[1]-atom2.coord[1])**2 + (atom1.coord[2]-atom2.coord[2])**2

def find_center_atom(m):
    return sorted([(atom.coord[0]**2 + atom.coord[1]**2 + atom.coord[2]**2, atom) for atom in m.atoms])[0][1]

def find_closest(m,atom):
    return sorted([((atom.coord[0]-atom2.coord[0])**2 + (atom.coord[1]-atom2.coord[1])**2 + (atom.coord[2]-atom2.coord[2])**2, atom2) for atom2 in m.atoms])[0][1]
    
def closest_list(m,atom):
    return sorted([((atom.coord[0]-atom2.coord[0])**2 + (atom.coord[1]-atom2.coord[1])**2 + (atom.coord[2]-atom2.coord[2])**2, atom2) for atom2 in m.atoms])

def main():
    m1 = sys.argv[1]
    m2 = sys.argv[2]
    m1 = Model(m1)
    m2 = Model(m2)
    c1 = find_center_atom(m1)
    c2 = find_center_atom(m2)

    f = open('temp.txt','w')
    g = open('good.txt','w')
    bestscore = 99999999
    bestmap = {}
    for i in range(10000):
        map = {}
        map[c1.id] = c2.id
        temp = [atom.id for atom in m1.atoms]
        shuffle(temp)
        #for atom in m1.atoms:
        for id in temp:
            atom = [x for x in m1.atoms if x.id==id][0]
            #if(atom.id == c1.id): continue
            if(atom.id == c1.id): continue
            #closest = find_closest(m2,atom)
            closest = closest_list(m2,atom)
            #print(atom.id, closest)
            for d,a in closest:
                if(a.id not in map.values()):
                    map[atom.id] = a.id
                    #print("Added {0}={1} to mapping".format(atom.id, a.id))
                    break
            #map[atom.id] = closest.id

        score = sum(dist2(m1.atoms[key],m2.atoms[val]) for key,val in map.items())
        score = score/m1.natoms
        if(score < bestscore):
            bestscore = score
            bestmap = {k:v for k,v in map.items()}
        f.write('{0}\n'.format(score))
        #if(score < 7.562):
        #if(score < 24):
        #    g.write('{0}'.format([(k,v) for k,v in map.items()]))
        #    print(score)
        #    print(map)
    f.close()
    g.close()
    print(bestscore)

    print('')
    map = bestmap
    score = sum(dist2(m1.atoms[key],m2.atoms[val]) for key,val in map.items())/m1.natoms
    #for key,val in bestmap.items():
    #    print(dist2(m1.atoms[key],m2.atoms[val]))
    print(score)
    print(map)

    print('')
    print('Permutation.fx(i,j) = 0;') # Over written in the next loop
    for k,v in map.items():
        print("Permutation.fx('{0}','{1}') = 1;".format(v+1,k+1))



if __name__ == '__main__':
    main()
