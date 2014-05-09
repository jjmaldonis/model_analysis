import sys
import numpy as np
from collections import defaultdict

def main():
    with open(sys.argv[1]) as f:
        content = f.readlines()

    content.pop(0)
    content.pop(0)
    content.pop(0)
    content.pop(0)
    content.pop(0)

    steps = []
    d = defaultdict(dict)
    for line in content:
        #print(line,len(line))
        if('.xyz' in line):
            sline = line.strip().split()
            step = int(sline[-1][:-4])
            steps.append(step)
        elif(len(line) != 2):
            sline = line.strip().split()
            count = int(sline.pop(-1))
            cna = eval(''.join(sline))
            d[cna][step] = count

    nsteps = len(steps)
    out = np.zeros((nsteps,len(d)),dtype=int)
    ld = list(d)
    for key in ld:
        for step in steps:
            #out[ld.index(key)][step] = d[key].get(step,0)
            out[steps.index(step)][ld.index(key)] = d[key].get(step,0)
    
    #print(' '.join(str(x).replace(" ","") for x in ld))
    print(' '.join(''.join(str(y) for y in x) for x in ld))
    print('')
    for i,step in enumerate(steps):
        print(str(step) + ' '.join(str(x) for x in out[i]))


if __name__ == '__main__':
    main()

