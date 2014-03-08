
import sys
from model import Model

def main():
    modelfile = sys.argv[1]

    m = Model(modelfile)
    m.write_cif()


if __name__ == "__main__":
    main()
