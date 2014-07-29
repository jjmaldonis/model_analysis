import sys


def simple_cubic():
    a = 2.50 # Lattice parameter in A
    world_min = -a*2
    world_max = a*2
    points = []
    for x in range(int(world_min/a),int(world_max/a)):
        for y in range(int(world_min/a),int(world_max/a)):
            for z in range(int(world_min/a),int(world_max/a)):
                points.append( (x*a, y*a, z*a) )
    print(points)

    f = open('sc.xyz', 'w')
    f.write('Comment: Simple Cubic Al Lattice\n')
    f.write(str(world_max*2) + '\t' + str(world_max*2) + '\t' + str(world_max*2) + '\n')
    for i in range(0,len(points)):
        line = ""
        for item in points[i]:
            line = line + '\t' + str(item)
        line = '13\t' + line.strip() + '\n' # Get rid of leading tab and add endline and call the atom Al
        f.write(line)
    f.write('-1')
    f.close()

def fcc():
    a = 2.5 # Lattice parameter in A
    world_min = -2*a
    world_max = 2*a
    points = []
    for x in range(int(round(world_min/a)),int(round(world_max/a))):
        for y in range(int(round(world_min/a)),int(round(world_max/a))):
            for z in range(int(round(world_min/a)),int(round(world_max/a))):
                points.append( (x*a, y*a, z*a) )
                points.append( ((0.5+x)*a, (0.5+y)*a, z*a) )
                points.append( ((0.5+x)*a, y*a, (0.5+z)*a) )
                points.append( (x*a, (0.5+y)*a, (0.5+z)*a) )

    print(len(points))
    print(points)

    f = open('fcc.xyz', 'w')
    f.write('Comment: FCC Lattice\n')
    f.write(str(world_max) + '\t' + str(world_max) + '\t' + str(world_max) + '\n')
    for i in range(0,len(points)):
        line = ""
        for item in points[i]:
            line = line + '\t' + str(item)
        line = '13\t' + line.strip() + '\n' # Get rid of leading tab and add endline and call the atom Al
        f.write(line)
    f.write('-1')
    f.close()


def main():
    #simple_cubic()
    fcc()

if __name__ == '__main__':
    main()
