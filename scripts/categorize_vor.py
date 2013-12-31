
# Categorizes the voronoi polyhedra according to the input parameter file, 
# which must be the first arg passed to this script.
# The second arg must be a _index.out file.
# Output is printed to screen.

import sys
class VorCats(object):
    def __init__(self,paramfile):
        super(VorCats,self).__init__()
        self.load_param_file(paramfile)

    def load_index_file(self,indexfile):
        # Open _index.out file.
        with open(indexfile) as f:
            self.index = f.readlines()

    def load_param_file(self,paramfile):
        # Now open the input parameter file.
        # Should be of the form:
        # Crystal:
        #     0,2,8,*
        with open(paramfile) as f:
            vp_params = [line.split(',') for line in f]
        # For each voronoi polyhedra 'structure', change it to an int if it's not a *
        self.vp_dict = {}
        for vps in vp_params:
            for i in range(0,len(vps)):
                vps[i] = vps[i].strip()
            if(len(vps) == 1):
                current_entry = vps[0]
                self.vp_dict[current_entry] = []
            elif(len(vps) == 4):
                for i in range(0,4):
                    if(vps[i] != '*'):
                        vps[i] = int(vps[i])
                self.vp_dict[current_entry].append(vps)
            else:
                sys.exit("ERROR IN INPUT VP PARAMETERS. STOPPING.")

    def save(self,indexes):
        indexes = [line.strip().split() for line in indexes]
        for line in indexes:
            for i in range(0,len(line)):
                try:
                    line[i] = int(line[i])
                except:
                    try:
                        line[i] = float(line[i])
                    except:
                        pass
        # Now all the lines are read in and in integer or float form
        # We care about columns 7-10, which are i3, i4, i5, i6.

        # Create a new dictionary to store all the atoms that are crystalline, icosahedra, etc.
        self.atom_dict = {}
        for key in self.vp_dict:
            self.atom_dict[key] = []

        # Sort atoms
        for line in indexes:
            for key in self.vp_dict:
                for vps in self.vp_dict[key]:
                    found = True
                    for i in range(0,4):
                        if(vps[i] != '*' and vps[i] != line[i+6]):
                            found = False
                    if(found):
                        self.atom_dict[key].append(line)

    def get_atom_dict(self):
        """ In the form:
            { 'Mixed:': [list of mixed atoms], 'Crystal-like:', [list of crystal-like atoms], etc}
            where each atom has the form of the _index.out lines:
            [id, .., .., .., .., .., i3, i4, i5, i6, .....] """
        return self.atom_dict
    def get_vp_dict(self):
        """ In the form:
            {'Crystal-like:': [[0, 4, 4, '*'], [0, 5, 2, '*']], 'Icosahedra-like:': [[0, 2, 8, '*'], [0, 1, 10, '*'], [0, 0, 12, 0], [0, 0, 12, 2]], 'Mixed:': [[0, 3, 6, '*'], [0, 3, 7, 2], [1, 2, 5, 4]]} """
        return self.vp_dict

    def printVorCats(self):
        print("Inputs:")
        for key in self.vp_dict:
            print(key)
            for vps in self.vp_dict[key]:
                print vps
        print("\nOutputs:")
        sum = 0
        for key in self.atom_dict:
            sum = sum + len(self.atom_dict[key])
        print("Total number of atoms categorized: " + str(sum))
        self.atom_dict = self.atom_dict
        for key in self.atom_dict:
            print(key + ' ' + str(len(self.atom_dict[key])))
            for line in self.atom_dict[key]:
                #print(line)
                for i in range(0,len(line)):
                    if(type(line[i]) != type('hi')):
                        line[i] = str(line[i])
                line = "\t".join(line)
                print(line)
            for line in self.atom_dict[key]:
                for i in range(0,len(line)):
                    try:
                        line[i] = int(line[i])
                        #print("Converted "+str(line[i]))
                    except:
                        try:
                            line[i] = float(line[i])
                            #print("Converted "+str(line[i]))
                        except:
                            pass
                            #print("Could not convert: "+str(line[i]))


def main():
    if(len(sys.argv) <= 2): sys.exit("\nERROR! Fix your inputs!\n\nArg 1:  input param file detailing each voronoi 'structure'.\nShould be of the form:\nCrystal:\n    0,2,8,*\n\nArg2: an _index.out file.\n\nOutput is printed to screen.")
    vor_cats = VorCats(sys.argv[1])
    vor_cats.load_index_file(sys.argv[2])
    vor_cats.save(vor_cats.index)
    vor_cats.printVorCats()
    #print(vor_cats.atom_dict)
    #print(vor_cats.vp_dict)

if __name__ == "__main__":
    main()
