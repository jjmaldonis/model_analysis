
# Sorts a _stat.out file and prints it to screen

import sys
from operator import itemgetter

def main():
    if(len(sys.argv) < 2): sys.exit("ERROR! You must pass a _stat.out file to be sorted.")
    with open( sys.argv[1]) as f:
        lines = [line.split() for line in f]
    for line in lines:
        for i in range(0,len(line)):
            try:
                line[i] = int(line[i])
            except:
                pass
    lines.sort(key=itemgetter(1))
    lines.reverse()

    #f = open(sys.argv[1]+'.new','w')
    for line in lines:
        for i in range(0,len(line)):
            if(type(line[i]) != type('hi')):
                line[i] = str(line[i])
        line = "\t".join(line)
        print(line)

    #f = open(sys.argv[1]+'.new','w')
    #for line in lines:
    #    for i in range(0,len(line)):
    #        if(type(line[i]) != type('hi')):
    #            line[i] = str(line[i])
    #    line = "\t".join(line)
    #    f.write(line+'\n')

if __name__ == "__main__":
    main()
