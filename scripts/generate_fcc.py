

def main():
    world_min = -10
    world_max = 10
    a = 2.5 # Lattice parameter in A
    points = []
    for x in range(world_min/a,world_max/a):
        for y in range(world_min/a,world_max/a):
            for z in range(world_min/a,world_max/a):
                points.append( (x*a, y*a, z*a) )
                points.append( ((0.5+x)*a, (0.5+y)*a, z*a) )
                points.append( ((0.5+x)*a, y*a, (0.5+z)*a) )
                points.append( (x*a, (0.5+y)*a, (0.5+z)*a) )

    print(len(points))
    print(points)

    #i = 0
    #while( i < len(points)):
    #    if(max(points[i]) > world_max or min(points[i]) < world_min):
    #        points.pop(i)
    #        i = i -1
    #    i = i + 1

    #print(len(points))
    #print(points)


    f = open('fcc.xyz', 'w')
    f.write('Comment: FCC Lattice\n')
    f.write(str(world_max*2) + '\t' + str(world_max*2) + '\t' + str(world_max*2) + '\n')
    for i in range(0,len(points)):
        line = ""
        for item in points[i]:
            line = line + '\t' + str(item)
        line = '13\t' + line.strip() + '\n' # Get rid of leading tab and add endline and call the atom Al
        f.write(line)
    f.write('-1')
    f.close()


if __name__ == "__main__":
    main()
