
# First arg is file from categorize_vor.py
# Second arg is model file (normal xyz format)
# Output is printed to screen
# Output is a sorted list of which atom is which VP type

import sys

def main():
    print("Arg 1: output from categorize_vor.py\nArg 2: model file in normal xyz format\n\nOutput printed to screen.\n")
    if(len(sys.argv) <= 1): sys.exit("ERROR! Fix your inputs!")

    f = open(sys.argv[1],'r')
    line = f.readline().strip()
    while(line != 'Outputs:'):
        line = f.readline().strip()
    line = f.readline().strip() # Total number of atoms categorized: ###


    # Creates a dictionary where each key contains values
    # that are the atom #s for that category.
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
    ln = -2
    for line in mf:
        ln = ln + 1
        if(ln > 0): # skip top two lines
            for key in cats:
                if(ln in cats[key]):
                    atoms[key].append(line)

    for key in atoms:
        print(key)
        for line in atoms[key]:
            print(line.strip())

    #outf = open('temp.out', 'w')
    #for key in atoms:
    #    outf.write(key+'\n')
    #    for line in atoms[key]:
    #        outf.write(line.strip() + '\n')

if __name__ == "__main__":
    main()
