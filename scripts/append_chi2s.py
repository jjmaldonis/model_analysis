import sys

def main():
    file1 = sys.argv[1]
    file2 = sys.argv[2]

    with open(file1) as f1:
        for line in f1:
            try:
                step,chi2,energy = tuple( line.strip().split() )
                print("{0} {1} {2}".format(step,chi2,energy))
            except:
                pass
    finalstep = int(step)
    with open(file2) as f2:
        for line in f2:
            #try:
            step,chi2,energy = tuple( line.strip().split() )
            print("{0} {1} {2}".format(int(step)+finalstep,chi2,energy))
            #except:
            #    pass


if __name__ == "__main__":
    main()
