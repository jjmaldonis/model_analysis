
import sys
from hutch import Hutch


def main():
    # sys.argv[1] should be an xyz cluster
    hutch = Hutch(sys.argv[1])
    dists = hutch.get_all_dists()
    for dist in dists:
        print dist

if __name__ == "__main__":
    main()
