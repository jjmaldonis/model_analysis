import sys

def main():
    files = sys.argv[1:]
    finalstep = 0
    for file in files:
        with open(file) as f:
            for line in f:
                try:
                    step,chi2,energy = tuple( line.strip().split() )
                    print("{0} {1} {2}".format(int(step)+finalstep,chi2,energy))
                except:
                    pass
        finalstep += int(step)

if __name__ == "__main__":
    main()
