""" This file generates a perfect <0,0,12,0> icosahedron.
    It is actually called a dodecahedron (platonic solid).
    You can specify the bond length between nearest neighbor atoms in the code """


def main():
    p = 1.61803398875
    b = 2. # Set interatomic bond distance here
    b = b * 0.5 * p
    coords = [[0,0,0]]

    coords.append([p,0,1./p])
    coords.append([-p,0,1./p])
    coords.append([-p,0,-1./p])
    coords.append([p,0,-1./p])

    coords.append([1./p, p, 0])
    coords.append([1./p, -p, 0])
    coords.append([-1./p, -p, 0])
    coords.append([-1./p, p, 0])

    coords.append([0, 1./p, p])
    coords.append([0, 1./p, -p])
    coords.append([0, -1./p, -p])
    coords.append([0, -1./p, p])

    coords.append([1,1,1])
    coords.append([1,-1,1])
    coords.append([-1,-1,1])
    coords.append([-1,1,1])

    coords.append([-1,1,-1])
    coords.append([1,1,-1])
    coords.append([1,-1,-1])
    coords.append([-1,-1,-1])

    coords = [ [b*x for x in c] for c in coords]

    f = open('icosahedron.xyz','w')
    f.write(str(len(coords))+'\n')
    f.write('comment\n')
    for c in coords:
        f.write('Si ' + ' '.join([str(x) for x in c]) + '\n')
    f.close()

if __name__ == '__main__':
    main()
