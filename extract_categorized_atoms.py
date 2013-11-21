import sys

def main():
    f = open(sys.argv[1],'r')
    line = f.readline().strip()
    while(line != 'Outputs:'):
        line = f.readline().strip()
    line = f.readline().strip() # Total number of atoms categorized: ###


    line = f.readline().strip() # Prep line
    cats = {}
    while True:
        if not line: break
        if( ':' in line ):
            cats[line] = []
            current_entry = line
        else:
            cats[current_entry].append(int(line.split("\t")[0]))
        line = f.readline().strip()
    f.close()


    atoms = {}
    for key in cats:
        atoms[key] = []

    mf = open(sys.argv[2],'r')
    ln = 0
    for line in mf:
        ln = ln + 1
        for key in cats:
            if(ln in cats[key]):
                atoms[key].append(line)

    outf = open('temp.out', 'w')
    for key in atoms:
        outf.write(key+'\n')
        for line in atoms[key]:
            outf.write(line.strip() + '\n')

if __name__ == "__main__":
    main()
